(ns varity.vcf-to-hgvs.cdna
  (:require [clj-hgvs.core :as hgvs]
            [clj-hgvs.mutation :as mut]
            [varity.ref-gene :as rg]
            [varity.util :refer [revcomp-bases]]
            [varity.vcf-to-hgvs.common :refer [diff-bases] :as common]))

(defn- repeat-info-forward
  [fa-rdr rg pos ins]
  (->> (common/read-sequence-stepwise-backward
        fa-rdr {:chr (:chr rg), :start (:tx-start rg), :end (dec pos)} 100)
       (map (fn [seq*]
              (let [nseq* (count seq*)]
                (if-let [[unit ref-repeat ins-repeat :as ri] (common/repeat-info seq* (inc nseq*) ins)]
                  (let [nunit (count unit)]
                    (if (> (* nunit ref-repeat) (- nseq* nunit))
                      false
                      ri))))))
       (remove false?)
       (first)))

(defn- repeat-info-backward
  [fa-rdr rg pos ins]
  (->> (common/read-sequence-stepwise
        fa-rdr {:chr (:chr rg), :start pos, :end (:tx-end rg)} 100)
       (map (fn [seq*]
              (let [nseq* (count seq*)]
                (if-let [[unit ref-repeat ins-repeat :as ri] (common/repeat-info (revcomp-bases seq*)
                                                                                 (inc nseq*)
                                                                                 (revcomp-bases ins))]
                  (let [nunit (count unit)]
                    (if (> (* nunit ref-repeat) (- nseq* nunit))
                      false
                      ri))))))
       (remove false?)
       (first)))

(defn- repeat-info*
  [fa-rdr rg pos ins]
  (case (:strand rg)
    "+" (repeat-info-forward fa-rdr rg pos ins)
    "-" (repeat-info-backward fa-rdr rg pos ins)))

(defn- mutation-type
  [fa-rdr rg pos ref alt]
  (if (re-matches #"[acgntACGNT]*" alt)
    (let [[ref-only alt-only offset _] (diff-bases ref alt)
          nrefo (count ref-only)
          nalto (count alt-only)
          [unit ref-repeat ins-repeat] (repeat-info* fa-rdr rg (+ pos offset) alt-only)]
      (cond
        (and (= nrefo 1) (= nalto 1)) :substitution
        (= ref-only (revcomp-bases alt-only)) :inversion
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
                          (cond-> ref (= strand "-") revcomp-bases)
                          (if (= ref alt) "=" ">")
                          (cond-> alt (= strand "-") revcomp-bases))))

(defn- dna-deletion
  [rg pos ref alt]
  (let [{:keys [strand]} rg
        [del _ offset _] (diff-bases ref alt)
        ndel (count del)
        left (+ pos offset)
        right (+ pos offset ndel -1)]
    (mut/dna-deletion (rg/cds-coord (case strand
                                      "+" left
                                      "-" right)
                                    rg)
                      (if (> ndel 1)
                        (rg/cds-coord (case strand
                                        "+" right
                                        "-" left)
                                      rg))
                      (cond-> (subs ref offset) (= strand "-") revcomp-bases))))

(defn- dna-duplication
  [rg pos ref alt]
  (let [{:keys [strand]} rg
        [_ ins offset _] (diff-bases ref alt)
        start (case strand
                "+" (+ pos offset (- (count ins)))
                "-" (+ pos offset (count ins) -1))
        end (case strand
              "+" (dec (+ start (count ins)))
              "-" (inc (- start (count ins))))]
    (mut/dna-duplication (rg/cds-coord start rg)
                         (rg/cds-coord end rg)
                         (cond-> ins (= strand "-") revcomp-bases))))

(defn- dna-insertion
  [rg pos ref alt]
  (let [{:keys [strand cds-start cds-end exon-ranges]} rg
        [_ ins offset _] (diff-bases ref alt)
        start (cond-> (+ pos offset) (= strand "+") dec)
        end (cond-> (+ pos offset) (= strand "-") dec)]
    (mut/dna-insertion (rg/cds-coord start rg)
                       (rg/cds-coord end rg)
                       (cond-> ins (= strand "-") revcomp-bases))))

(defn- dna-inversion
  [rg pos ref alt]
  (let [{:keys [strand]} rg
        [inv _ offset _] (diff-bases ref alt)
        left (+ pos offset)
        right (+ pos offset (count inv) -1)
        start (case strand
                "+" left
                "-" right)
        end (case strand
              "+" right
              "-" left)]
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
                                   "+" left
                                   "-" right)
                                 rg)
                   (if (> ndel 1)
                     (rg/cds-coord (case strand
                                     "+" right
                                     "-" left)
                                   rg))
                   (cond-> del (= strand "-") revcomp-bases)
                   (cond-> ins (= strand "-") revcomp-bases))))

(defn- dna-repeated-seqs
  [fa-rdr rg pos ref alt]
  (let [{:keys [strand]} rg
        [_ ins offset _] (diff-bases ref alt)
        [unit ref-repeat ins-repeat] (repeat-info* fa-rdr rg (+ pos offset) ins)
        nunit (count unit)
        start (case strand
                "+" (+ pos offset (- (* nunit ref-repeat)))
                "-" (+ pos offset (* nunit ref-repeat) -1))
        end (case strand
              "+" (dec (+ start nunit))
              "-" (inc (- start nunit)))]
    (mut/dna-repeated-seqs (rg/cds-coord start rg)
                           (rg/cds-coord end rg)
                           unit
                           (+ ref-repeat ins-repeat))))

(defn- mutation
  [fa-rdr rg pos ref alt]
  (case (mutation-type fa-rdr rg pos ref alt)
    :substitution (dna-substitution rg pos ref alt)
    :deletion (dna-deletion rg pos ref alt)
    :duplication (dna-duplication rg pos ref alt)
    :insertion (dna-insertion rg pos ref alt)
    :inversion (dna-inversion rg pos ref alt)
    :indel (dna-indel rg pos ref alt)
    :repeated-seqs (dna-repeated-seqs fa-rdr rg pos ref alt)))

(defn ->hgvs
  [{:keys [pos ref alt]} fa-rdr rg]
  (hgvs/hgvs (:name rg)
             :cdna
             (mutation fa-rdr rg pos ref alt)))
