(ns varity.vcf-to-hgvs.coding-dna
  (:require [clojure.pprint :as pp]
            [clojure.string :as string]
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.core :as hgvs]
            [clj-hgvs.mutation :as mut]
            [cljam.util.sequence :as util-seq]
            [cljam.io.sequence :as cseq]
            [proton.string :as pstring]
            [varity.ref-gene :as rg]
            [varity.vcf-to-hgvs.common :refer [diff-bases] :as common]))

(defn- repeat-info-forward
  [seq-rdr rg pos alt type]
  (->> (common/read-sequence-stepwise-backward
        seq-rdr
        {:chr (:chr rg)
         :start (- (:tx-start rg) rg/max-tx-margin)
         :end (dec (cond-> pos
                     (= type :del) (+ (count alt))))}
        100)
       (map (fn [seq*]
              (let [nseq* (count seq*)
                    pos* (case type
                           :ins (inc nseq*)
                           :del (inc (- nseq* (count alt))))]
                (if-let [[unit ref-repeat :as ri] (common/repeat-info
                                                   seq* pos* alt type)]
                  (let [nunit (count unit)]
                    (if (> (* nunit ref-repeat) (- nseq* nunit))
                      false
                      ri))))))
       (remove false?)
       (first)))

(defn- repeat-info-backward
  [seq-rdr rg pos alt type]
  (->> (common/read-sequence-stepwise
        seq-rdr
        {:chr (:chr rg)
         :start pos
         :end (+ (:tx-end rg) rg/max-tx-margin)}
        100)
       (map (fn [seq*]
              (let [nseq* (count seq*)
                    pos* (case type
                           :ins (inc nseq*)
                           :del (inc (- nseq* (count alt))))]
                (if-let [[unit ref-repeat :as ri] (common/repeat-info
                                                   (util-seq/revcomp seq*)
                                                   pos*
                                                   (util-seq/revcomp alt)
                                                   type)]
                  (let [nunit (count unit)]
                    (if (> (* nunit ref-repeat) (- nseq* nunit))
                      false
                      ri))))))
       (remove false?)
       (first)))

(defn- repeat-info*
  [seq-rdr rg pos ref-only alt-only]
  (when-let [[alt type] (cond
                          (and (empty? ref-only) (seq alt-only)) [alt-only :ins]
                          (and (seq ref-only) (empty? alt-only)) [ref-only :del])]
    (case (:strand rg)
      :forward (repeat-info-forward seq-rdr rg pos alt type)
      :reverse (repeat-info-backward seq-rdr rg pos alt type))))

(defn- mutation-type
  [seq-rdr rg pos ref alt {:keys [prefer-deletion?]}]
  (if (re-matches #"[acgntACGNT]*" alt)
    (let [[ref-only alt-only offset _] (diff-bases ref alt)
          nrefo (count ref-only)
          nalto (count alt-only)
          [unit ref-repeat alt-repeat] (repeat-info* seq-rdr rg (+ pos offset)
                                                     ref-only alt-only)]
      (cond
        (or (= nrefo nalto 0) (= nrefo nalto 1)) :substitution
        (and prefer-deletion? (pos? nrefo) (zero? nalto)) :deletion
        (= ref-only (util-seq/revcomp alt-only)) :inversion
        (and (some? unit) (= ref-repeat 1) (= alt-repeat 2)) :duplication
        (and (some? unit) (pos? alt-repeat)
             ;; In the protein coding region, repeat descriptions are used only
             ;; for repeat units with a length which is a multiple of 3, i.e.
             ;; which can not affect the reading frame.
             ;; See https://varnomen.hgvs.org/recommendations/DNA/variant/repeated/#note
             (if (and (rg/in-exon? (+ pos offset) rg)
                      (rg/in-exon? (+ pos offset (count ref-only) -1) rg))
               (zero? (mod (count unit) 3))
               true)) :repeated-seqs
        (and (pos? nrefo) (zero? nalto)) :deletion
        (and (pos? nrefo) (pos? nalto)) :indel
        (some? unit) :duplication
        (and (zero? nrefo) (pos? nalto)) :insertion
        :else (throw (ex-info "Unsupported variant" {:type ::unsupported-variant}))))
    (throw (ex-info "Unsupported variant" {:type ::unsupported-variant}))))

(defn- dna-substitution
  [rg pos ref alt]
  (let [{:keys [strand]} rg
        type (if (= ref alt) "=" ">")]
    (mut/dna-substitution (rg/cds-coord pos rg)
                          (cond-> ref (= strand :reverse) util-seq/revcomp)
                          type
                          (when-not (= type "=")
                            (cond-> alt (= strand :reverse) util-seq/revcomp)))))

(defn- dna-deletion
  [rg pos ref alt]
  (let [{:keys [strand]} rg
        [del _ offset _] (diff-bases ref alt)
        ndel (count del)
        left (+ pos offset)
        right (+ pos offset ndel -1)]
    (mut/dna-deletion (rg/cds-coord (case strand
                                      :forward left
                                      :reverse right)
                                    rg)
                      (if (> ndel 1)
                        (rg/cds-coord (case strand
                                        :forward right
                                        :reverse left)
                                      rg))
                      (cond-> (subs ref offset) (= strand :reverse) util-seq/revcomp))))

(defn- dna-duplication
  [rg pos ref alt]
  (let [{:keys [strand]} rg
        [_ ins offset _] (diff-bases ref alt)
        start (case strand
                :forward (+ pos offset (- (count ins)))
                :reverse (+ pos offset (count ins) -1))
        end (case strand
              :forward (dec (+ start (count ins)))
              :reverse (inc (- start (count ins))))]
    (mut/dna-duplication (rg/cds-coord start rg)
                         (rg/cds-coord end rg)
                         (cond-> ins (= strand :reverse) util-seq/revcomp))))

(defn- dna-insertion
  [rg pos ref alt]
  (let [{:keys [strand]} rg
        [_ ins offset _] (diff-bases ref alt)
        start (cond-> (+ pos offset) (= strand :forward) dec)
        end (cond-> (+ pos offset) (= strand :reverse) dec)]
    (mut/dna-insertion (rg/cds-coord start rg)
                       (rg/cds-coord end rg)
                       (cond-> ins (= strand :reverse) util-seq/revcomp))))

(defn- dna-inversion
  [rg pos ref alt]
  (let [{:keys [strand]} rg
        [inv _ offset _] (diff-bases ref alt)
        left (+ pos offset)
        right (+ pos offset (count inv) -1)
        start (case strand
                :forward left
                :reverse right)
        end (case strand
              :forward right
              :reverse left)]
    (mut/dna-inversion (rg/cds-coord start rg)
                       (rg/cds-coord end rg))))

;; e.g. ref alt
;;      CAG CTC
;;      CA  CTC
;;      CAG CT
(defn- dna-indel
  [rg pos ref alt]
  (let [{:keys [strand]} rg
        [del ins offset _] (diff-bases ref alt)
        ndel (count del)
        left (+ pos offset)
        right (+ pos offset ndel -1)]
    (mut/dna-indel (rg/cds-coord (case strand
                                   :forward left
                                   :reverse right)
                                 rg)
                   (if (> ndel 1)
                     (rg/cds-coord (case strand
                                     :forward right
                                     :reverse left)
                                   rg))
                   (cond-> del (= strand :reverse) util-seq/revcomp)
                   (cond-> ins (= strand :reverse) util-seq/revcomp))))

(defn- dna-repeated-seqs
  [seq-rdr rg pos ref alt]
  (let [{:keys [strand]} rg
        [ref-only alt-only offset _] (diff-bases ref alt)
        [unit ref-repeat alt-repeat] (repeat-info* seq-rdr rg (+ pos offset)
                                                   ref-only alt-only)
        nunit (count unit)
        start (case strand
                :forward (+ pos offset (- (* nunit (min ref-repeat alt-repeat))))
                :reverse (+ pos offset (* nunit ref-repeat) -1))
        end (case strand
              :forward (dec (+ start nunit))
              :reverse (inc (- start nunit)))]
    (mut/dna-repeated-seqs (rg/cds-coord start rg)
                           (rg/cds-coord end rg)
                           unit
                           alt-repeat)))

(defn- mutation
  [seq-rdr rg pos ref alt options]
  (case (mutation-type seq-rdr rg pos ref alt options)
    :substitution (dna-substitution rg pos ref alt)
    :deletion (dna-deletion rg pos ref alt)
    :duplication (dna-duplication rg pos ref alt)
    :insertion (dna-insertion rg pos ref alt)
    :inversion (dna-inversion rg pos ref alt)
    :indel (dna-indel rg pos ref alt)
    :repeated-seqs (dna-repeated-seqs seq-rdr rg pos ref alt)))

(defn ->hgvs
  ([variant seq-rdr rg]
   (->hgvs variant seq-rdr rg {}))
  ([{:keys [pos ref alt]} seq-rdr rg options]
   (hgvs/hgvs (:name rg)
              :coding-dna
              (mutation seq-rdr rg pos ref alt options))))

(defn- sequence-pstring
  [ref-seq alt-seq start end {:keys [pos ref alt]} rg]
  (let [pos (case (:strand rg)
              :forward pos
              :reverse (if (= (count ref) (count alt))
                         (+ pos (count ref) -1)
                         (+ pos (count ref))))
        pos* (case (:strand rg)
               :forward (- pos start)
               :reverse (- end pos))
        ref-seq (cond-> ref-seq (= (:strand rg) :reverse) util-seq/revcomp)
        alt-seq (cond-> alt-seq (= (:strand rg) :reverse) util-seq/revcomp)
        [ref-up ref ref-down] (pstring/split-at ref-seq
                                                [pos*
                                                 (+ pos* (count ref))])
        [alt-up alt alt-down] (pstring/split-at alt-seq
                                                [pos*
                                                 (+ pos* (count alt))])
        nmut (max (count ref) (count alt))
        ticks (->> (iterate #(+ % 10) (inc (- start (mod start 10))))
                   (take-while #(<= % end))
                   (remove #(< % start)))
        ticks (cond->> ticks (= (:strand rg) :reverse) reverse)
        ticks-cds (map #(coord/format (rg/cds-coord % rg)) ticks)
        tick-intervals (conj
                        (case (:strand rg)
                          :forward (->> ticks
                                        (map #(+ % (if (< pos %)
                                                     (max 0 (- (count alt) (count ref)))
                                                     0)))
                                        (map #(- % start))
                                        (partition 2 1)
                                        (mapv #(- (second %) (first %))))
                          :reverse (->> ticks
                                        (map #(+ % (if (> pos %)
                                                     (max 0 (- (count alt) (count ref)))
                                                     0)))
                                        (map #(- end %))
                                        (partition 2 1)
                                        (mapv #(- (second %) (first %)))))
                        0)
        ticks-offset (case (:strand rg)
                       :forward (- (first ticks) start)
                       :reverse (- end (first ticks)))]
    (string/join
     \newline
     [(apply pp/cl-format nil
             (apply str "~V@T" (repeat (count ticks-cds) "~VA"))
             ticks-offset (interleave tick-intervals ticks-cds))
      (pp/cl-format nil "~A~VA~A" ref-up nmut ref ref-down)
      (pp/cl-format nil "~A~VA~A" alt-up nmut alt alt-down)
      (pp/cl-format nil "~V@T~V@{~A~:*~}" pos*  nmut "^")])))

(defn- lower-case-intron
  [seq* start rg]
  (->> (map vector seq* (range start (+ start (count seq*))))
       (map (fn [[c i]]
              (if (rg/in-exon? i rg) c (string/lower-case c))))
       (apply str)))

(defn debug-string
  [{:keys [pos ref alt]} seq-rdr rg]
  (let [start (max 1 (- pos 20))
        end (+ pos (count ref) 20)
        ref (lower-case-intron ref pos rg)
        alt (lower-case-intron alt pos rg)
        ref-seq (-> (cseq/read-sequence seq-rdr {:chr (:chr rg), :start start, :end end})
                    (lower-case-intron start rg))
        alt-seq (common/alt-sequence ref-seq start pos ref alt)
        [ticks refs alts muts] (string/split-lines
                                (sequence-pstring ref-seq alt-seq start end
                                                  {:pos pos :ref ref :alt alt} rg))]
    (string/join \newline [(str "         " ticks)
                           (str "  c.ref: " refs)
                           (str "  c.alt: " alts)
                           (str "         " muts)])))
