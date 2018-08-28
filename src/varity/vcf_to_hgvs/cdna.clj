(ns varity.vcf-to-hgvs.cdna
  (:require [clj-hgvs.core :as hgvs]
            [clj-hgvs.mutation :as mut]
            [cljam.util.sequence :as util-seq]
            [varity.ref-gene :as rg]
            [varity.vcf-to-hgvs.common :refer [diff-bases] :as common]))

(defn- repeat-info-forward
  [seq-rdr rg pos ins]
  (->> (common/read-sequence-stepwise-backward
        seq-rdr
        {:chr (:chr rg)
         :start (- (:tx-start rg) rg/max-tx-margin)
         :end (dec pos)}
        100)
       (map (fn [seq*]
              (let [nseq* (count seq*)]
                (if-let [[unit ref-repeat :as ri] (common/repeat-info seq* (inc nseq*) ins)]
                  (let [nunit (count unit)]
                    (if (> (* nunit ref-repeat) (- nseq* nunit))
                      false
                      ri))))))
       (remove false?)
       (first)))

(defn- repeat-info-backward
  [seq-rdr rg pos ins]
  (->> (common/read-sequence-stepwise
        seq-rdr
        {:chr (:chr rg)
         :start pos
         :end (+ (:tx-end rg) rg/max-tx-margin)}
        100)
       (map (fn [seq*]
              (let [nseq* (count seq*)]
                (if-let [[unit ref-repeat :as ri] (common/repeat-info (util-seq/revcomp seq*)
                                                                      (inc nseq*)
                                                                      (util-seq/revcomp ins))]
                  (let [nunit (count unit)]
                    (if (> (* nunit ref-repeat) (- nseq* nunit))
                      false
                      ri))))))
       (remove false?)
       (first)))

(defn- repeat-info*
  [seq-rdr rg pos ins]
  (case (:strand rg)
    :forward (repeat-info-forward seq-rdr rg pos ins)
    :reverse (repeat-info-backward seq-rdr rg pos ins)))

(defn- mutation-type
  [seq-rdr rg pos ref alt]
  (if (re-matches #"[acgntACGNT]*" alt)
    (let [[ref-only alt-only offset _] (diff-bases ref alt)
          nrefo (count ref-only)
          nalto (count alt-only)
          [unit ref-repeat ins-repeat] (repeat-info* seq-rdr rg (+ pos offset) alt-only)]
      (cond
        (and (= nrefo 1) (= nalto 1)) :substitution
        (= ref-only (util-seq/revcomp alt-only)) :inversion
        (and (pos? nrefo) (zero? nalto)) :deletion
        (and (pos? nrefo) (pos? nalto)) :indel
        (some? unit) (cond
                       (and (= ref-repeat 1) (= ins-repeat 1)) :duplication
                       (or (> ref-repeat 1) (> ins-repeat 1)) :repeated-seqs)
        (and (zero? nrefo) (pos? nalto)) :insertion
        :else (throw (IllegalArgumentException. "Unsupported variant"))))
    (throw (IllegalArgumentException. "Unsupported variant"))))

(defn- dna-substitution
  [rg pos ref alt]
  (let [{:keys [strand]} rg]
    (mut/dna-substitution (rg/cds-coord pos rg)
                          (cond-> ref (= strand :reverse) util-seq/revcomp)
                          (if (= ref alt) "=" ">")
                          (cond-> alt (= strand :reverse) util-seq/revcomp))))

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
        [_ ins offset _] (diff-bases ref alt)
        [unit ref-repeat ins-repeat] (repeat-info* seq-rdr rg (+ pos offset) ins)
        nunit (count unit)
        start (case strand
                :forward (+ pos offset (- (* nunit ref-repeat)))
                :reverse (+ pos offset (* nunit ref-repeat) -1))
        end (case strand
              :forward (dec (+ start nunit))
              :reverse (inc (- start nunit)))]
    (mut/dna-repeated-seqs (rg/cds-coord start rg)
                           (rg/cds-coord end rg)
                           unit
                           (+ ref-repeat ins-repeat))))

(defn- mutation
  [seq-rdr rg pos ref alt]
  (case (mutation-type seq-rdr rg pos ref alt)
    :substitution (dna-substitution rg pos ref alt)
    :deletion (dna-deletion rg pos ref alt)
    :duplication (dna-duplication rg pos ref alt)
    :insertion (dna-insertion rg pos ref alt)
    :inversion (dna-inversion rg pos ref alt)
    :indel (dna-indel rg pos ref alt)
    :repeated-seqs (dna-repeated-seqs seq-rdr rg pos ref alt)))

(defn ->hgvs
  [{:keys [pos ref alt]} seq-rdr rg]
  (hgvs/hgvs (:name rg)
             :cdna
             (mutation seq-rdr rg pos ref alt)))
