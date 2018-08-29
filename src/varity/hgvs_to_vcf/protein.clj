(ns varity.hgvs-to-vcf.protein
  (:require [clj-hgvs.core :as hgvs]
            [clj-hgvs.mutation :as mut]
            [cljam.io.sequence :as cseq]
            [cljam.util.sequence :as util-seq]
            [varity.codon :as codon]
            [varity.ref-gene :as rg])
  (:import [clj_hgvs.mutation ProteinSubstitution]))

(defn- pos-candidates
  [coord rg]
  (let [cds-pos (inc (* (dec (:position coord)) 3))]
    (->> (range cds-pos (+ cds-pos 3))
         (keep #(rg/cds->genomic-pos % rg))
         (sort))))

(defn- one-character-substituted?
  [ref-seq alt-seq pos]
  (and (not= ref-seq alt-seq)
       (= (assoc (vec ref-seq) pos nil)
          (assoc (vec alt-seq) pos nil))))

(defn- vcf-variants
  [seq-rdr {:keys [chr strand] :as rg} mut*]
  (if (instance? ProteinSubstitution mut*)
    (let [pos-cands (pos-candidates (:coord mut*) rg)]
      (if (seq pos-cands)
        (let [ref-codon1 (->> pos-cands
                              (map #(cseq/read-sequence seq-rdr {:chr chr, :start %, :end %}))
                              (apply str))
              ref-codon (cond-> ref-codon1 (= strand "-") util-seq/revcomp)
              palt (mut/->short-amino-acid (:alt mut*))
              codon-cands (codon/amino-acid->codons palt)]
          (->> codon-cands
               (keep (fn [codon*]
                      (->> pos-cands
                           (map-indexed
                            (fn [idx pos]
                              (if (one-character-substituted? ref-codon codon* idx)
                                {:chr chr
                                 :pos pos
                                 :ref (cond-> (str (nth ref-codon idx))
                                        (= strand "-") util-seq/revcomp)
                                 :alt (cond-> (str (nth codon* idx))
                                        (= strand "-") util-seq/revcomp)})))
                           (remove nil?))))
               (flatten)))))
    (throw (IllegalArgumentException. "Unsupported mutation"))))

(defn ->vcf-variants
  [hgvs seq-rdr rg]
  (vcf-variants seq-rdr rg (:mutation hgvs)))

(defn- vcf-variants-with-cdna-hgvs
  [seq-rdr {:keys [chr strand] :as rg} mut*]
  (if (instance? ProteinSubstitution mut*)
    (let [pos-cands (pos-candidates (:coord mut*) rg)]
      (if (seq pos-cands)
        (let [ref-codon1 (->> pos-cands
                              (map #(cseq/read-sequence seq-rdr {:chr chr, :start %, :end %}))
                              (apply str))
              ref-codon (cond-> ref-codon1 (= strand "-") util-seq/revcomp)
              palt (mut/->short-amino-acid (:alt mut*))
              codon-cands (codon/amino-acid->codons palt)
              pos-cands (cond-> pos-cands (= strand "-") reverse)]
          (->> codon-cands
               (keep (fn [codon*]
                      (->> pos-cands
                           (map-indexed
                            (fn [idx pos]
                              (if (one-character-substituted? ref-codon codon* idx)
                                (let [ref (str (nth ref-codon idx))
                                      alt (str (nth codon* idx))]
                                 {:vcf {:chr chr
                                        :pos pos
                                        :ref (cond-> ref (= strand "-") util-seq/revcomp)
                                        :alt (cond-> alt (= strand "-") util-seq/revcomp)}
                                  :cdna (hgvs/hgvs (:name rg) :cdna (mut/dna-substitution (rg/cds-coord pos rg)
                                                                                  ref
                                                                                  (if (= ref alt) "=" ">")
                                                                                  alt))}))))
                           (remove nil?))))
               (flatten)))))
    (throw (IllegalArgumentException. "Unsupported mutation"))))

(defn ->vcf-variants-with-cdna-hgvs
  [hgvs seq-rdr rg]
  (vcf-variants-with-cdna-hgvs seq-rdr rg (:mutation hgvs)))
