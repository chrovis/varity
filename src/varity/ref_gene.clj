(ns varity.ref-gene
  "Handles refGene.txt(.gz) content."
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [clj-hgvs.coordinate :as coord]
            [cljam.io.sequence :as cseq]
            [cljam.util.chromosome :refer [normalize-chromosome-key]]
            [cljam.util.sequence :as util-seq]
            [proton.core :refer [as-long]]
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

(defn load-ref-genes
  "Loads f (e.g. refGene.txt(.gz)), returning the all contents as a sequence."
  [f]
  (with-open [rdr (io/reader (util/compressor-input-stream f))]
    (->> (line-seq rdr)
         (map parse-ref-gene-line)
         doall)))

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
   (get-in rgidx (if (re-find #"^(NC|LRG|NG|NM|NR|NP)_" s)
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

(defn seek-gene-region
  "Seeks chr:pos through exon entries in refGene and returns those indices"
  ([chr pos rgidx]
   (seek-gene-region chr pos rgidx nil))
  ([chr pos rgidx name]
   ;; TODO seek intron region
   ;; TODO seek UTR-5 or UTR-3 region
   (->> (if name
          (ref-genes name rgidx)
          (ref-genes chr pos rgidx))
        (map (fn [rg]
               (let [sgn (case (:strand rg)
                           :forward identity
                           :reverse reverse)
                     exon-ranges (sgn (:exon-ranges rg))
                     exon-idx (->> exon-ranges
                                   (keep-indexed (fn [i [s e]] (if (<= s pos e) i)))
                                   first)
                     intron-ranges (if exon-idx nil
                                       (->> (:exon-ranges rg)
                                            ((fn [r]
                                               (loop [intron-ranges []
                                                      left (first r)
                                                      right (second r)
                                                      ranges (rest r)]
                                                 (if (empty? (first ranges)) intron-ranges
                                                     (recur (conj intron-ranges [(inc (second left)) (dec (first right))])
                                                            right
                                                            (second ranges)
                                                            (rest ranges))))))
                                            (sgn)))
                     intron-idx (->> intron-ranges
                                     (keep-indexed (fn [i [s e]] (if (<= s pos e) i)))
                                     first)
                     region-type (let [txs (:tx-start rg)
                                       txe (:tx-end rg)]
                                   (cond
                                     (< pos txs) (first (sgn '({:type "UTR-5"} {:type "UTR-3"})))
                                     (> pos txe) (second (sgn '({:type "UTR-5"} {:type "UTR-3"})))
                                     :else (if exon-idx
                                             {:type "exon" :idx exon-idx :count (count exon-ranges)}
                                             {:type "intron" :idx intron-idx :count (count intron-ranges)}
                                             ))
                                   )]

                 {:type (:type region-type) ; "exon", "intron", "UTR-5" or "UTR-3"
                  :index (if (:idx region-type)
                           (inc (:idx region-type))
                           nil)
                  :exon-idx exon-idx
                  :intron-idx intron-idx
                  :count (:count region-type)
                  :gene rg}))))))

;;; Calculation of CDS coordinate
;;;
;;; cf. http://varnomen.hgvs.org/bg-material/numbering/

(defn- nearest-edge-and-offset
  [pos {:keys [strand exon-ranges]}]
  (->> exon-ranges
       (map (fn [[s e]] [[s (- pos s)] [e (- pos e)]]))
       (apply concat)
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
    (+ base-pos (cond-> (:offset coord) (= strand :reverse) (-)))))
