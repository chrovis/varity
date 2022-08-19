(ns varity.hgvs-to-vcf.coding-dna
  (:require clj-hgvs.coordinate
            clj-hgvs.mutation
            [cljam.io.sequence :as cseq]
            [cljam.util.sequence :as util-seq]
            [varity.ref-gene :as rg]))

(defn- cds-coord->genomic-pos
  [coord rg]
  (cond
    (instance? clj_hgvs.coordinate.CodingDNACoordinate coord)
    (rg/cds-coord->genomic-pos coord rg)

    (or (instance? clj_hgvs.coordinate.UnknownCoordinate coord)
        (instance? clj_hgvs.coordinate.UncertainCoordinate coord))
    (throw (ex-info "Ambiguous coordinate" {:type ::ambiguous-coordinate
                                            :coordinate coord}))

    :else
    (throw (IllegalArgumentException.
            "coord must be clj-hgvs.coordinate/CDNACoordinate record"))))

;; AGT => ["AGT"]
;; AGTAGT => ["AGT" "AGTAGT"]
(defn- repeat-units
  [s]
  (->> (range (count s))
       (map inc)
       (map #(partition-all % s))
       (filter #(apply = %))
       (map first)
       (mapv #(apply str %))))

(defn- repeat-count
  [seq-rdr chr pos unit strand]
  (if (seq unit)
    (->> (iterate #((case strand
                      :forward +
                      :reverse -)
                    % (count unit)) pos)
         (map #(cseq/read-sequence seq-rdr {:chr chr
                                            :start %
                                            :end (+ % (count unit) -1)}))
         (take-while #(= % unit))
         count)
    (throw (IllegalArgumentException. "empty unit"))))

(defmulti vcf-variant (fn [mut* seq-rdr rg] (class mut*)))

(defmethod vcf-variant clj_hgvs.mutation.DNASubstitution
  [mut* _ {:keys [chr strand] :as rg}]
  (if-let [pos (cds-coord->genomic-pos (:coord mut*) rg)]
    (let [alt (if (= (:type mut*) "=")
                (:ref mut*)
                (:alt mut*))]
      {:chr chr
       :pos pos
       :ref (cond-> (:ref mut*)
              (= strand :reverse) (util-seq/revcomp))
       :alt (cond-> alt
              (= strand :reverse) (util-seq/revcomp))})))

(defmethod vcf-variant clj_hgvs.mutation.DNADeletion
  [mut* seq-rdr {:keys [chr strand] :as rg}]
  (let [coord-end (or (:coord-end mut*) (:coord-start mut*))
        start (cds-coord->genomic-pos (case strand
                                        :forward (:coord-start mut*)
                                        :reverse coord-end)
                                      rg)
        end (cds-coord->genomic-pos (case strand
                                      :forward coord-end
                                      :reverse (:coord-start mut*))
                                    rg)]
    (if (and start end)
      (when-let [ref (cseq/read-sequence seq-rdr {:chr chr, :start (dec start), :end end})]
        {:chr chr
         :pos (dec start)
         :ref ref
         :alt (subs ref 0 1)}))))

(defmethod vcf-variant clj_hgvs.mutation.DNADuplication
  [mut* seq-rdr {:keys [chr strand] :as rg}]
  (let [coord-end (or (:coord-end mut*) (:coord-start mut*))
        start (cds-coord->genomic-pos (case strand
                                        :forward (:coord-start mut*)
                                        :reverse coord-end)
                                      rg)
        end (cds-coord->genomic-pos (case strand
                                      :forward coord-end
                                      :reverse (:coord-start mut*))
                                    rg)]
    (if (and start end)
      (when-let [dup (cseq/read-sequence seq-rdr {:chr chr, :start start, :end end})]
        (let [repeat (repeat-count seq-rdr chr start dup (case strand
                                                           :forward :reverse
                                                           :reverse :forward))
              base-pos (case strand
                         :forward (dec (- start (* (count dup) (dec repeat))))
                         :reverse (dec start))
              base (cseq/read-sequence seq-rdr {:chr chr, :start base-pos, :end base-pos})]
          {:chr chr
           :pos base-pos
           :ref base
           :alt (str base dup)})))))

(defmethod vcf-variant clj_hgvs.mutation.DNAInsertion
  [mut* seq-rdr {:keys [chr strand] :as rg}]
  (if-let [start (cds-coord->genomic-pos (case strand
                                           :forward (:coord-start mut*)
                                           :reverse (:coord-end mut*))
                                         rg)]
    (when-let [ref (cseq/read-sequence seq-rdr {:chr chr, :start start, :end start})]
      {:chr chr
       :pos start
       :ref ref
       :alt (str ref (cond-> (:alt mut*)
                       (= strand :reverse) (util-seq/revcomp)))})))

(defmethod vcf-variant clj_hgvs.mutation.DNAInversion
  [mut* seq-rdr {:keys [chr strand] :as rg}]
  (let [start (cds-coord->genomic-pos (case strand
                                        :forward (:coord-start mut*)
                                        :reverse (:coord-end mut*))
                                      rg)
        end (cds-coord->genomic-pos (case strand
                                      :forward (:coord-end mut*)
                                      :reverse (:coord-start mut*))
                                    rg)]
    (if (and start end)
      (when-let [ref (cseq/read-sequence seq-rdr {:chr chr, :start (dec start), :end end})]
        {:chr chr
         :pos (dec start)
         :ref ref
         :alt (str (first ref) (util-seq/revcomp (subs ref 1)))}))))

(defmethod vcf-variant clj_hgvs.mutation.DNAIndel
  [mut* seq-rdr {:keys [chr strand] :as rg}]
  (let [coord-end (or (:coord-end mut*) (:coord-start mut*))
        start (cds-coord->genomic-pos (case strand
                                        :forward (:coord-start mut*)
                                        :reverse coord-end)
                                      rg)
        end (cds-coord->genomic-pos (case strand
                                      :forward coord-end
                                      :reverse (:coord-start mut*))
                                    rg)]
    (if (and start end)
      (when-let [ref (cseq/read-sequence seq-rdr {:chr chr, :start (dec start), :end end})]
        (if (or (nil? (:ref mut*))
                (= (subs ref 1) (cond-> (:ref mut*)
                                  (= strand :reverse) (util-seq/revcomp))))
          {:chr chr
           :pos (dec start)
           :ref ref
           :alt (str (first ref) (cond-> (:alt mut*)
                                   (= strand :reverse) (util-seq/revcomp)))})))))

(defn- valid-repeat?
  [seq-rdr rg pos unit]
  (let [r-start ((case (:strand rg)
                   :forward -
                   :reverse +)
                 pos (count unit))
        r-end (dec (+ r-start (count unit)))]
    (and (= (count (repeat-units unit)) 1)
         (not= (cseq/read-sequence seq-rdr {:chr (:chr rg)
                                            :start r-start
                                            :end r-end})
               unit))))

(defmethod vcf-variant clj_hgvs.mutation.DNARepeatedSeqs
  [mut* seq-rdr {:keys [chr strand] :as rg}]
  (let [start* (cds-coord->genomic-pos (:coord-start mut*) rg)
        end* (cond
               (:coord-end mut*) (cds-coord->genomic-pos (:coord-end mut*) rg)
               (:ref mut*) (when start*
                             (+ start* (cond-> (dec (count (:ref mut*)))
                                         (= strand :reverse) -)))
               :else start*)]
    (if (and start* end*)
      (let [[start end] (cond-> [start* end*] (= strand :reverse) reverse)
            unit (cseq/read-sequence seq-rdr {:chr chr, :start start, :end end})]
        (if (valid-repeat? seq-rdr rg start unit)
          (let [ref-repeat (repeat-count seq-rdr chr start unit strand)
                alt-repeat (:ncopy mut*)
                base-pos (case strand
                           :forward (dec start)
                           :reverse (dec (- start (* (count unit) (dec ref-repeat)))))
                base (cseq/read-sequence seq-rdr {:chr chr, :start base-pos, :end base-pos})]
            {:chr chr
             :pos base-pos
             :ref (if (<= ref-repeat alt-repeat)
                    base
                    (apply str base (repeat (- ref-repeat alt-repeat) unit)))
             :alt (if (<= ref-repeat alt-repeat)
                    (apply str base (repeat (- alt-repeat ref-repeat) unit))
                    base)}))))))

(defn ->vcf-variant
  [hgvs seq-rdr rg]
  (vcf-variant (:mutation hgvs) seq-rdr rg))
