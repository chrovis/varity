(ns varity.ref-gene
  "Handles refGene.txt(.gz) and ncbiRefSeq.txt(.gz) content."
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [clj-hgvs.coordinate :as coord]
            [cljam.io.sequence :as cseq]
            [cljam.util.chromosome :refer [normalize-chromosome-key]]
            [cljam.util.sequence :as util-seq]
            [proton.core :refer [as-long as-float]]
            [varity.util :as util]))

;;; Utility

(defn in-cds?
  "Returns true if pos is in the coding region, false otherwise."
  [pos {:keys [cds-start cds-end]}]
  (<= cds-start pos cds-end))

(defn in-exon?
  "Returns true if pos is in the exon region, false otherwise."
  [pos {:keys [exon-ranges]}]
  (->> exon-ranges
       (some (fn [[s e]] (<= s pos e)))
       true?))

;;; Parser

(defn- parse-exon-pos
  [s]
  (map as-long (string/split s #",")))

(defn- exon-ranges
  [starts ends]
  (->> (map vector starts ends)
       (sort-by first)
       vec))

(defn- parse-ref-gene-line
  [s]
  (as-> (zipmap [:bin :name :chr :strand :tx-start :tx-end :cds-start :cds-end
                 :exon-count :exon-start :exon-end :score :name2 :cds-start-stat
                 :cds-end-stat :exon-frames]
                (string/split s #"\t")) m
    (update m :bin as-long)
    (update m :chr normalize-chromosome-key)
    (update m :strand #(case (first %) \+ :forward \- :reverse))
    (update m :tx-start (comp inc as-long))
    (update m :tx-end as-long)
    (update m :cds-start (comp inc as-long))
    (update m :cds-end as-long)
    (update m :exon-count as-long)
    (update m :exon-start #(map inc (parse-exon-pos %)))
    (update m :exon-end parse-exon-pos)
    (update m :score as-long)
    (assoc m :exon-ranges (exon-ranges (:exon-start m) (:exon-end m)))
    (dissoc m :exon-start :exon-end)
    (update m :cds-start-stat keyword)
    (update m :cds-end-stat keyword)
    (update m :exon-frames parse-exon-pos)))

(defn- load-ncbi-file
  [f]
  (with-open [rdr (io/reader (util/compressor-input-stream f))]
    (->> (line-seq rdr)
         (map parse-ref-gene-line)
         (filter #(re-find #"^(NM|NR)_.+$" (:name %)))
         doall)))

(defn load-ref-genes
  {:deprecated "0.8.0"
   :doc "DEPRECATED: Loads f (e.g. refGene.txt(.gz)), returning the all contents as a sequence."}
  [f]
  (load-ncbi-file f))

(defn load-ref-seqs
  "Loads f (e.g. ncbiRefSeq.txt(.gz)), returning the all contents as a sequence."
  [f]
  (load-ncbi-file f))

(defn- ->gencode-attr
  [attr-str kv-sep]
  (->> (string/split attr-str #";")
       (map string/trim)
       (mapcat #(string/split % kv-sep))
       (partition-all 2)
       (map (fn [[k v]]
              [k (string/replace v #"\"" "")]))
       (into {})))

(def ^:private gencode-attr-kv-sep
  {:gtf #"\s"
   :gff3 #"="})

(defn- parse-gencode-line
  [s & {:keys [attr-kv-sep] :or {attr-kv-sep (gencode-attr-kv-sep :gtf)}}]
  (when-not (= \# (first s))
    (let [row (string/split s #"\t")]
      (merge
       {:seqname (nth row 0)
        :source (nth row 1)
        :feature (nth row 2)
        :start (as-long (nth row 3))
        :end (as-long (nth row 4))
        :score (as-float (nth row 5))
        :strand (nth row 6)
        :frame (nth row 7)}
       (when (= (count row) 9)
         {:attribute (->gencode-attr (nth row 8) attr-kv-sep)})))))

(defn- ->feature-map
  ([features]
   (->feature-map features (zipmap [:transcript :exon :cds] (repeat {}))))
  ([gtf-lines data]
   (reduce (fn [data {:keys [seqname feature attribute] :as gtf}]
             (let [base (merge (select-keys gtf [:start :end])
                               {:chr seqname})]
               (if-let [t-id (get attribute "transcript_id")]
                 (let [t-id (re-find #"ENST\d+\.\d+" t-id)]
                   (case feature
                     "transcript"
                     (assoc-in data
                               [:transcript t-id]
                               (merge base
                                      {:name t-id
                                       :name2 (get attribute "gene_name")
                                       :gene-id (when-let [g-id (get attribute "gene_id")]
                                                  (re-find #"ENSG\d+\.\d+" g-id))
                                       :strand (:strand gtf)
                                       :score (:score gtf)}))

                     "exon"
                     (update-in data
                                [:exon t-id]
                                conj
                                (merge base
                                       {:exon-number (as-long (get attribute "exon_number"))
                                        :frame (:frame gtf)
                                        :strand (:strand gtf)}))
                     "CDS"
                     (update-in data [:cds t-id] conj base)

                     "stop_codon"
                     (assoc-in data [:stop-codon t-id] base)

                     data))
                 data)))
           data
           gtf-lines)))

(defn- extend-cds
  "Extend 3'-most cds's `:end` or `:start` depending on the `strand` value
   if the cds doesn't include stop codon in its value"
  [strand {stop-codon-start :start stop-codon-end :end :as stop-codon} cdss]
  (let [{last-cds-start :start last-cds-end :end} (last cdss)]
    (if (and (not= stop-codon-start last-cds-start)
             (not= stop-codon-end last-cds-end))
      (update (vec cdss)
              (dec (count cdss))
              (if (= strand :forward)
                #(merge % (select-keys stop-codon [:end]))
                #(merge % (select-keys stop-codon [:start]))))
      cdss)))

(defn- ->region
  [feature-map [transcript-id transcript]]
  (let [exons (->> (get-in feature-map [:exon transcript-id])
                   (filter :exon-number)
                   (sort-by :exon-number))
        [exons strand] (case (:strand transcript)
                         "+" [exons :forward]
                         "-" [(reverse exons) :reverse])
        cds (extend-cds strand
                        (get-in feature-map [:stop-codon transcript-id])
                        (get-in feature-map [:cds transcript-id]))]
    (cond-> transcript
      (seq exons) (assoc :exon-count (count exons)
                         :exon-start (:start (first exons))
                         :exon-end (:end (last exons))
                         :exon-ranges (mapv #(vector (:start %) (:end %)) exons)
                         :exon-frames (mapv :frame exons))
      (empty? exons) (assoc :exon-count 0
                            :exon-start nil
                            :exon-end nil
                            :exon-ranges []
                            :exon-frames [])
      (seq cds) (assoc :cds-start (apply min (map :start cds))
                       :cds-end (apply max (map :end cds)))
      (empty? cds) (assoc :cds-start (:start transcript)
                          :cds-end (:end transcript))

      true (assoc :tx-start (:start transcript)
                  :tx-end (:end transcript)
                  :strand strand))))

(defn load-gencode
  [f parse-line]
  (with-open [rdr (io/reader (util/compressor-input-stream f))]
    (let [feature-map (->> (line-seq rdr)
                           (keep parse-line)
                           ->feature-map)]
      (->> (:transcript feature-map)
           (map #(->region feature-map %))
           doall))))

(defn load-gtf
  [f]
  (load-gencode f parse-gencode-line))

(defn load-gff3
  [f]
  (load-gencode f #(parse-gencode-line % :attr-kv-sep (:gff3 gencode-attr-kv-sep))))

;; Indexing
;; --------

(defn- round-int
  "e.g. (round-int 138 100) => 100"
  [a b]
  (- a (mod a b)))

(def ^:private pos-index-block 1000000)

(def max-tx-margin 10000)

(defn- locus-index
  [rgs]
  (->> (group-by :chr rgs)
       (map (fn [[chr sub-rgs]]
              (let [fs (round-int (- (apply min (map :tx-start sub-rgs))
                                     max-tx-margin)
                                  pos-index-block)
                    le (round-int (+ (apply max (map :tx-end sub-rgs))
                                     max-tx-margin)
                                  pos-index-block)]
                [chr (loop [s fs, ret {}]
                       (if (<= s le)
                         (let [e (+ s pos-index-block)
                               rgs* (filter (fn [{:keys [tx-start tx-end]}]
                                              (and (<= (- tx-start max-tx-margin) e)
                                                   (<= s (+ tx-end max-tx-margin))))
                                            sub-rgs)]
                           (recur e (assoc ret [s e] rgs*)))
                         ret))])))
       (into {})))

(defn- ref-seq-index
  [rgs]
  (group-by :name rgs))

(defn- gene-index
  [rgs]
  (group-by :name2 rgs))

(defrecord RefGeneIndex [locus ref-seq gene])

(defn index
  "Creates refGene index for search."
  [rgs]
  (RefGeneIndex. (locus-index rgs)
                 (ref-seq-index rgs)
                 (gene-index rgs)))

(defn ref-genes
  "Searches refGene entries with ref-seq, gene or (chr, pos) using index,
  returning results as sequence. See also varity.ref-gene/index."
  ([s rgidx]
   (get-in rgidx (if (re-find #"^ENST|^(NC|LRG|NG|NM|NR|NP)_" s)
                   [:ref-seq s]
                   [:gene s])))
  ([chr pos rgidx] (ref-genes chr pos rgidx 0))
  ([chr pos rgidx tx-margin]
   {:pre [(<= 0 tx-margin max-tx-margin)]}
   (let [pos-r (round-int pos pos-index-block)]
     (->> (get-in rgidx [:locus
                         (normalize-chromosome-key chr)
                         [pos-r (+ pos-r pos-index-block)]])
          (filter (fn [{:keys [tx-start tx-end]}]
                    (<= (- tx-start tx-margin) pos (+ tx-end tx-margin))))))))

(defn in-any-exon?
  "Returns true if chr:pos is located in any ref-gene exon, else false."
  [chr pos rgidx]
  (->> (ref-genes chr pos rgidx)
       (some #(in-exon? pos %))
       (true?)))

(defn tx-region
  "Returns a genomic region of the given gene."
  [{:keys [chr tx-start tx-end strand]}]
  {:chr chr, :start tx-start, :end tx-end, :strand strand})

(defn cds-region
  "Returns a genomic region of a coding sequence of the given gene. Returns nil
  if the gene is a non-coding RNA."
  [{:keys [chr cds-start cds-end strand]}]
  (when (<= cds-start cds-end)
    {:chr chr, :start cds-start, :end cds-end, :strand strand}))

(defn exon-seq
  "Returns a lazy sequence of regions corresponding to each exon in a gene. The
  exons are ordered by their index, thus they're reversed in genomic coordinate
  if the refGene record is on the reverse strand."
  [{:keys [chr strand exon-ranges]}]
  (let [exon-count (count exon-ranges)]
    (->> exon-ranges
         ((if (= :reverse strand) reverse identity))
         (map-indexed
          (fn [i [s e]]
            {:exon-index (inc i), ;; 1-origin
             :exon-count exon-count,
             :chr chr,
             :start s,
             :end e,
             :strand strand})))))

(defn cds-seq
  "Returns a lazy sequence of exons included in a coding region of a
  `ref-gene-record`. Note that exons outside of the CDS are removed and
  partially overlapping ones are cropped in the result. Returns nil if the record
  is a non-coding RNA."
  [{:keys [cds-start cds-end] :as ref-gene-record}]
  (when (<= cds-start cds-end)
    (keep
     (fn [{:keys [start end] :as exon}]
       (when-not (or (< end cds-start) (< cds-end start))
         (-> exon
             (update :start max cds-start)
             (update :end min cds-end))))
     (exon-seq ref-gene-record))))

(defn ^String read-exon-sequence
  "Reads a base sequence of an `exon` from `seq-rdr`."
  [seq-rdr {:keys [strand] :as exon}]
  (cond-> (cseq/read-sequence seq-rdr exon)
    (= :reverse strand) util-seq/revcomp))

(defn ^String read-transcript-sequence
  "Reads a DNA base sequence of a `ref-gene-record` from `seq-rdr`. The sequence
  contains 5'-UTR, CDS and 3'-UTR."
  [seq-rdr ref-gene-record]
  (->> ref-gene-record
       exon-seq
       (map (partial read-exon-sequence seq-rdr))
       string/join))

(defn ^String read-coding-sequence
  "Reads a coding sequence of a ref-gene record `ref-gene-record` from
  `seq-rdr`. Returns nil if the gene is a non-coding RNA."
  [seq-rdr ref-gene-record]
  (some->> ref-gene-record
           cds-seq
           (map (partial read-exon-sequence seq-rdr))
           string/join))

(defn exon-ranges->intron-ranges
  [exon-ranges]
  (mapv (fn [[left right]]
          [(inc (second left)) (dec (first right))])
        (partition 2 1 exon-ranges)))

(defn seek-gene-region
  "Seeks chr:pos through exon entries in refGene and returns those indices"
  ([chr pos rgidx]
   (seek-gene-region chr pos rgidx nil))
  ([chr pos rgidx name]
   (->> (if name
          (ref-genes name rgidx)
          (ref-genes chr pos rgidx))
        (map (fn [rg]
               (let [sgn (case (:strand rg)
                           :forward identity
                           :reverse reverse)
                     exon-ranges (:exon-ranges rg)
                     exon-idx (->> (sgn exon-ranges)
                                   (keep-indexed (fn [i [s e]] (if (<= s pos e) i)))
                                   first)
                     intron-ranges (if-not exon-idx
                                     (->> exon-ranges
                                          (exon-ranges->intron-ranges)
                                          (sgn)))
                     intron-idx (->> intron-ranges
                                     (keep-indexed (fn [i [s e]] (if (<= s pos e) i)))
                                     first)
                     utr (cond
                           (< pos (:cds-start rg)) (first (sgn '({:region "UTR-5"}
                                                                 {:region "UTR-3"})))
                           (> pos (:cds-end rg)) (second (sgn '({:region "UTR-5"}
                                                                {:region "UTR-3"})))
                           :else nil)
                     exon-intron (cond
                                   exon-idx {:region "exon" :index (if exon-idx (inc exon-idx)) :count (count exon-ranges)}
                                   intron-idx {:region "intron" :index (if intron-idx (inc intron-idx)) :count (count intron-ranges)})
                     regions (remove nil? (vector utr exon-intron))]
                 {:regions regions :gene rg}))))))

;;; Calculation of CDS coordinate
;;;
;;; cf. http://varnomen.hgvs.org/bg-material/numbering/

(defn- nearest-edge-and-offset
  [pos {:keys [strand exon-ranges]}]
  (->> exon-ranges
       (mapcat (fn [[s e]]
                 [[s (- pos s) :left]
                  [e (- pos e) :right]]))
       (sort-by (fn [[e ^long o]]
                  [(Math/abs o) (case strand
                                  :forward e
                                  :reverse (- e))]))
       first))

(defn- exon-pos
  [pos strand exon-ranges]
  (->> exon-ranges
       (map (fn [[s e]]
              (case strand
                :forward (max (min (inc (- pos s)) (inc (- e s))) 0)
                :reverse (max (min (inc (- e pos)) (inc (- e s))) 0))))
       (reduce +)))

(defn cds-pos
  [pos {:keys [strand cds-start cds-end exon-ranges]}]
  {:pre [(<= (ffirst exon-ranges) pos (second (last exon-ranges)))
         (<= cds-start cds-end)]}
  (let [pos* (exon-pos pos strand exon-ranges)
        cds-start* (exon-pos cds-start strand exon-ranges)
        cds-end* (exon-pos cds-end strand exon-ranges)
        [start* end*] (cond-> [cds-start* cds-end*]
                        (= strand :reverse) reverse)]
    (cond
      (< pos* start*) [(- start* pos*) :upstream]
      (< end* pos*) [(- pos* end*) :downstream]
      :else [(inc (- pos* start*)) nil])))

(defn cds-coord
  "Converts the genomic position into the coding DNA coordinate. The return
  value is clj-hgvs.coordinate/CodingDNACoordinate record."
  [pos rg]
  (let [[pos* offset] (if (in-exon? pos rg)
                        [pos 0]
                        (nearest-edge-and-offset pos rg))
        tx-edge? (or (= pos* (:tx-start rg)) (= pos* (:tx-end rg)))
        offset (case (:strand rg)
                 :forward offset
                 :reverse (- offset))
        [cds-pos* region] (cds-pos pos* rg)
        [cds-pos* offset] (cond
                            (and tx-edge? (= region :upstream)) [(- cds-pos* offset) 0]
                            (and tx-edge? (= region :downstream)) [(+ cds-pos* offset) 0]
                            :else [cds-pos* offset])]
    (coord/coding-dna-coordinate cds-pos* offset region)))

;;; Calculation of genomic coordinate

(defn cds->genomic-pos
  ([cds-pos rg] (cds->genomic-pos cds-pos nil rg))
  ([cds-pos region {:keys [strand cds-start cds-end exon-ranges]}]
   (let [exon-poss (mapcat (fn [[s e]] (range s (inc e))) exon-ranges)
         cds-poss (cond-> (filter #(<= cds-start % cds-end) exon-poss)
                    (= strand :reverse) reverse)
         utr-up-poss (case strand
                       :forward (reverse (filter #(< % cds-start) exon-poss))
                       :reverse (filter #(< cds-end %) exon-poss))
         utr-down-poss (case strand
                         :forward (filter #(< cds-end %) exon-poss)
                         :reverse (reverse (filter #(< % cds-start) exon-poss)))
         upstream-poss (if (seq utr-up-poss)
                         (concat utr-up-poss
                                 (case strand
                                   :forward (iterate dec (dec (last utr-up-poss)))
                                   :reverse (iterate inc (inc (last utr-up-poss)))))
                         utr-up-poss)
         downstream-poss (if (seq utr-down-poss)
                           (concat utr-down-poss
                                   (case strand
                                     :forward (iterate inc (inc (last utr-down-poss)))
                                     :reverse (iterate dec (dec (last utr-down-poss)))))
                           utr-down-poss)]
     (case region
       nil (nth cds-poss (dec cds-pos) nil)
       :upstream (nth upstream-poss (dec cds-pos) nil)
       :downstream (nth downstream-poss (dec cds-pos) nil)))))

(defn cds-coord->genomic-pos
  "Converts the coding DNA coordinate into the genomic position. coord must be
  clj-hgvs.coordinate/CodingDNACoordinate record."
  [coord {:keys [strand] :as rg}]
  (if-let [base-pos (cds->genomic-pos (:position coord) (:region coord) rg)]
    (let [[_ edge-offset edge-side] (nearest-edge-and-offset base-pos rg)
          offset (cond-> (or (:offset coord) 0)
                   (= strand :reverse) -)]
      (cond
        (zero? offset)
        base-pos

        (and (zero? edge-offset)
             (or (and (= edge-side :left) (neg? offset))
                 (and (= edge-side :right) (pos? offset))))
        (+ base-pos offset)

        :else
        (throw (ex-info "The coordinate is invalid for the refGene."
                        {:type ::invalid-coordinate, :coordinate coord}))))))
