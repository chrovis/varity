(ns varity.fusion
  "NOTE: This algorithm was initially developed at Cancer Precision
  Medicine Center, Japanese Foundation for Cancer Research and Xcoo
  re-implemented it with Clojure."
  (:require [clojure.string :as string]
            [cljam.io.vcf.util :as vcf-util]
            [cljam.util.sequence :as util-seq]
            [varity.ref-gene :as rg]
            [varity.codon :as codon]))

(defn- aligned-genes
  "Takes breakpoints, genes and inserted sequences. Returns `nil` if the genes
  are misaligned and makes no sense. Otherwise, returns the input arguments,
  rearranging them in the order of 5' to 3' of the transcript if necessary."
  [{r1 :retained :as bp1} {r2 :retained :as bp2}
   {s1 :strand :as gene1} {s2 :strand :as gene2} inserted-seq]
  (let [inserted-seq (cond-> inserted-seq
                       (and inserted-seq (= :reverse s1)) util-seq/revcomp)]
    (cond
      (and (= s1 (case r1 :L :forward :R :reverse))
           (= s2 (case r2 :L :reverse :R :forward)))
      [bp1 bp2 gene1 gene2 inserted-seq]

      (and (= s1 (case r1 :L :reverse :R :forward))
           (= s2 (case r2 :L :forward :R :reverse)))
      [bp2 bp1 gene2 gene1 inserted-seq])))

(defn- exon-intron-seq
  "Returns a lazy sequence of regions corresponding to each exon and intron in a
  gene. The exons and introns are ordered by their index, thus they're reversed
  in genomic coordinate if the refGene record is on the reverse strand."
  [gene]
  (->> (rg/exon-seq gene)
       (map #(assoc % :type :exon :index (:exon-index %) :gene gene))
       (partition-all 2 1)
       (mapcat (fn [[e1 e2]]
                 (cond-> [e1]
                   e2 (conj (let [[e1* e2*] (sort-by :start [e1 e2])]
                              (assoc e1
                                     :type :intron
                                     :start (inc (:end e1*))
                                     :end (dec (:start e2*))))))))))

(defn- splice-5'-part [{:keys [type start end]
                        {:keys [cds-start cds-end]} :gene :as m}]
  ;; exclude 5'-UTR and 3'-UTR of 5'-side part of the fusion gene
  (when (and (= :exon type) (<= cds-start end) (<= start cds-end))
    (-> m
        (update :start max cds-start)
        (update :end min cds-end))))

(defn- splice-3'-part [{:keys [type start end]
                        {:keys [cds-start cds-end strand]} :gene :as m}]
  ;; include 5'-UTR, exclude 3'-UTR of 3'-side part of the fusion gene
  (when (and (= :exon type) (case strand
                              :forward (<= start cds-end)
                              :reverse (<= cds-start end)))
    (case strand
      :forward (update m :end min cds-end)
      :reverse (update m :start max cds-start))))

(defn- in-5'-utr? [pos {:keys [strand tx-start tx-end cds-start cds-end]}]
  (case strand
    :forward (<= tx-start pos (dec cds-start))
    :reverse (<= (inc cds-end) pos tx-end)))

(defn- in-3'-utr? [pos {:keys [strand tx-start tx-end cds-start cds-end]}]
  (case strand
    :forward (<= (inc cds-end) pos tx-end)
    :reverse (<= tx-start pos (dec cds-start))))

(defn- utr
  [pos gene]
  (let [fp? (in-5'-utr? pos gene)
        tp? (in-3'-utr? pos gene)]
    (assert (not (and fp? tp?)))
    (cond fp? :5'-utr tp? :3'-utr)))

(defn- partial-5'-part-exon [pos {{:keys [strand]} :gene :as exon}]
  (when (= (:type exon) :exon)
    (cond-> exon
      (= strand :forward) (update :end min pos)
      (= strand :reverse) (update :start max pos))))

(defn- partial-3'-part-exon [pos {{:keys [strand]} :gene :as exon}]
  (when (= (:type exon) :exon)
    (cond-> exon
      (= strand :forward) (update :start max pos)
      (= strand :reverse) (update :end min pos))))

(defn- transcript-regions
  "Returns a map that consists of the following key-values:
  - `:fusion?` `true` iff actually forms a fused transcript
  - `:breakpoint-regions` [`overlapped-1` `overlapped-2`] with `:utr` in each
  - `:transcript-regions` A sequence of regions corresponding to each exon in a
  fusion gene. The sequence may contain {:inserted-seq \"...\"}."
  [exon-introns-1 pos-1 overlapped-1
   inserted-seq
   overlapped-2 pos-2 exon-introns-2]
  (let [exons-1 (vec (keep splice-5'-part exon-introns-1))
        exons-2 (vec (keep splice-3'-part exon-introns-2))
        broken-1 (partial-5'-part-exon pos-1 (splice-5'-part overlapped-1))
        broken-2 (partial-3'-part-exon pos-2 (splice-3'-part overlapped-2))
        utr-1 (utr pos-1 (:gene overlapped-1))
        utr-2 (utr pos-2 (:gene overlapped-2))
        inserted-seq (some->> inserted-seq (array-map :inserted-seq))]
    (merge
     {:breakpoint-regions [(assoc overlapped-1 :utr utr-1)
                           (assoc overlapped-2 :utr utr-2)]}
     (cond
       ;; 5'-UTR, no start codon
       (= utr-1 :5'-utr)
       {:fusion? false, :transcript-regions nil}

       ;; 3'-UTR, already seen a stop codon
       (= utr-1 :3'-utr)
       {:fusion? false, :transcript-regions
        (cond-> exons-1 broken-1 (conj broken-1))}

       ;; mate starts from 3'-UTR, can't determine if it stops
       (if broken-1 (= utr-2 :3'-utr) (empty? exons-2))
       {:fusion? false, :transcript-regions
        (cond-> exons-1 broken-1 (conj broken-1))}

       :else
       (case [(:type overlapped-1) (:type overlapped-2)]
         ;; partial exon + inserted seq + partial exon
         [:exon :exon]
         {:fusion? true,
          :transcript-regions
          (cond-> (conj exons-1 broken-1)
            inserted-seq (conj inserted-seq)
            true (into (cons broken-2 exons-2)))}

         ;; likely to have a stop codon in the intron
         [:exon :intron]
         {:fusion? false, :transcript-regions
          (cond-> (conj exons-1 broken-1)
            inserted-seq (conj inserted-seq))}

         ;; continues to the next exon of 3' part
         ;; ignoring inserted seq
         ([:intron :exon] [:intron :intron])
         {:fusion? true, :transcript-regions (into exons-1 exons-2)})))))

(defn- split-at-breakpoint [pos xs]
  (split-with (fn [{:keys [start end]}] (not (<= start pos end))) xs))

;;;; public

(defn fusion-genes
  "Returns a lazy sequence of annotations of all possible combinations of fusion
  genes. An annotation is a map containing the following key-values:
  - `:breakpoints` Given breakpoints in the order of 5'-part, 3'-part.
  - `:genes` Candidate genes for each of the `:breakpoints`.
  - `:fusion?` `true` iff the genes form a fused transcript.
  - `:breakpoint-regions` Genic regions for each of the `:breakpoints`.
  - `:transcript-regions` A sequence of regions corresponding to each exon of a
  resulting transcript. May not be fused or may be `nil` if no transcription is
  expected."
  [rg-index
   {[{c1 :chr p1 :pos :as bp1}
     {c2 :chr p2 :pos :as bp2}] :breakpoints
    is :inserted-seq}]
  (for [gene1* (sort-by :name (rg/ref-genes c1 p1 rg-index))
        gene2* (sort-by :name (rg/ref-genes c2 p2 rg-index))
        :let [[bp1 bp2 gene1 gene2 is] (aligned-genes bp1 bp2 gene1* gene2* is)]
        :when bp1
        :let [split #(split-at-breakpoint (:pos %1) (exon-intron-seq %2))
              [xs1 [x1 & _]] (split bp1 gene1)
              [_ [x2 & xs2]] (split bp2 gene2)]]
    (merge {:breakpoints [bp1 bp2], :genes [gene1 gene2]}
           (transcript-regions xs1 (:pos bp1) x1 is x2 (:pos bp2) xs2))))

(defn in-frame?
  "Checks if the reading frame of given regions are in-frame."
  [transcript-regions]
  (-> (map (fn [m]
             (if (contains? m :inserted-seq)
               (count (:inserted-seq m))
               (inc (- (:end m) (:start m))))))
      (transduce + 0 transcript-regions)
      (mod 3)
      zero?))

(defn transcript
  "Translates a sequence of regions into an amino acid sequence."
  [seq-reader transcript-regions]
  (->> transcript-regions
       (map #(or (:inserted-seq %) (rg/read-exon-sequence seq-reader %)))
       (apply str)
       codon/amino-acid-sequence))

(defn fusion-transcripts
  "Like `fusion-genes` but returns only in-frame and actually fused genes. An
  amino acid sequence is added as `:transcript` for each candidates."
  [seq-reader rg-index bp]
  (->> bp
       (fusion-genes rg-index)
       (filter (every-pred :fusion? (comp in-frame? :transcript-regions)))
       (map (juxt identity
                  (comp (partial transcript seq-reader) :transcript-regions)))
       (map (fn [[m t]] (assoc m :transcript t)))))

(defn parse-vcf-breakpoints
  "Parses a VCF-style breakend string into a map used in `varity.fusion`."
  [chr pos ref-base alt-bnd]
  (let [{:keys [bases join strand]
         mate-chr :chr
         mate-pos :pos} (vcf-util/parse-breakend alt-bnd)
        mate-retained (case [join strand]
                        [:after :forward] :R
                        [:after :reverse] :L
                        [:before :forward] :L
                        [:before :reverse] :R)]
    (when-not (case join
                :before (= (first ref-base) (last bases))
                :after (= (first ref-base) (first bases)))
      (throw (ex-info "Reference base mismatch"
                      {:chr chr, :pos pos, :ref ref-base, :alt alt-bnd})))
    {:breakpoints
     [{:chr chr, :pos pos, :retained (case join :after :L :before :R)}
      {:chr mate-chr, :pos mate-pos, :retained mate-retained}],
     :inserted-seq (not-empty
                    (case join
                      :before (subs bases 0 (dec (count bases)))
                      :after (subs bases 1)))}))

(defn variant->breakpoints
  "Converts each allele in a VCF-style variant into a map used in
  `varity.fusion`. Returns a lazy sequence of the same length with alleles."
  [{:keys [chr pos ref alt] {:keys [SVTYPE END]} :info}]
  (for [a alt]
    (let [{allele-type :type symbol-id :id} (vcf-util/inspect-allele ref a)
          [alt-type _subtype] (some-> symbol-id (string/split #":"))]
      (cond
        (and (integer? END)
             (or (and (= allele-type :id) (= "DUP" alt-type))
                 (= "DUP" SVTYPE)))
        {:breakpoints [{:chr chr, :pos pos, :retained :R}
                       {:chr chr, :pos END, :retained :L}],
         :inserted-seq nil}

        (and (integer? END)
             (or (and (= allele-type :id) (= "DEL" alt-type))
                 (= "DEL" SVTYPE)))
        {:breakpoints [{:chr chr, :pos pos, :retained :L}
                       {:chr chr, :pos END, :retained :R}],
         :inserted-seq nil}

        (= allele-type :breakend)
        (parse-vcf-breakpoints chr pos ref a)))))
