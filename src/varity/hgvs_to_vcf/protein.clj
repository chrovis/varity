(ns varity.hgvs-to-vcf.protein
  (:require [clj-hgvs.core :as hgvs]
            [clj-hgvs.mutation :as mut]
            [cljam.fasta :as fa]
            [varity.codon :as codon]
            [varity.ref-gene :as rg]
            [varity.util :refer [revcomp-bases]])
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
  [fa-rdr {:keys [chr strand] :as rg} mut*]
  (if (instance? ProteinSubstitution mut*)
    (let [pos-cands (pos-candidates (:coord mut*) rg)]
      (if (seq pos-cands)
        (let [ref-codon1 (->> pos-cands
                              (map #(fa/read-sequence fa-rdr {:chr chr, :start %, :end %}))
                              (apply str))
              ref-codon (cond-> ref-codon1 (= strand "-") revcomp-bases)
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
                                        (= strand "-") revcomp-bases)
                                 :alt (cond-> (str (nth codon* idx))
                                        (= strand "-") revcomp-bases)})))
                           (remove nil?))))
               (flatten)))))
    (throw (IllegalArgumentException. "Unsupported mutation"))))

(defn ->vcf-variants
  [hgvs fa-rdr rg]
  (vcf-variants fa-rdr rg (:mutation hgvs)))

(defn- vcf-variants-with-cdna-hgvs
  [fa-rdr {:keys [chr strand] :as rg} mut*]
  (if (instance? ProteinSubstitution mut*)
    (let [pos-cands (pos-candidates (:coord mut*) rg)]
      (if (seq pos-cands)
        (let [ref-codon1 (->> pos-cands
                              (map #(fa/read-sequence fa-rdr {:chr chr, :start %, :end %}))
                              (apply str))
              ref-codon (cond-> ref-codon1 (= strand "-") revcomp-bases)
              palt (mut/->short-amino-acid (:alt mut*))
              codon-cands (codon/amino-acid->codons palt)]
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
                                        :ref (cond-> ref (= strand "-") revcomp-bases)
                                        :alt (cond-> alt (= strand "-") revcomp-bases)}
                                  :cdna (hgvs/hgvs nil :cdna (mut/dna-substitution (rg/cds-coord pos rg)
                                                                                   ref
                                                                                   (if (= ref alt) "=" ">")
                                                                                   alt))}))))
                           (remove nil?))))
               (flatten)))))
    (throw (IllegalArgumentException. "Unsupported mutation"))))

(defn ->vcf-variants-with-cdna-hgvs
  [hgvs fa-rdr rg]
  (vcf-variants-with-cdna-hgvs fa-rdr rg (:mutation hgvs)))
