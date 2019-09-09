(ns varity.hgvs-to-vcf.protein
  (:require [clojure.string :as string]
            [clojure.tools.logging :as log]
            [clj-hgvs.core :as hgvs]
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

(defn- trim-left-right [positions ref alt]
  (->> (map vector positions ref alt)
       (drop-while #(apply = (rest %)))
       reverse
       (drop-while #(apply = (rest %)))
       reverse
       (apply map vector)))

(defn- contiguous? [xs]
  (when-first [f xs]
    (let [l (last xs)]
      (or (= xs (range f (inc l) 1))
          (= xs (range f (dec l) -1))))))

(defn- codon->variant [positions ref alt]
  (when-not (= ref alt) ;; not a variant
    (let [[ps r a] (seq (trim-left-right positions ref alt))]
      (if (contiguous? ps)
        [ps (string/join r) (string/join a)]
        (log/warnf
         "A candidate variant crossing an intron is unsupported: %s %s %s"
         ps r a)))))

(defn- vcf-variants
  [seq-rdr {:keys [chr strand] :as rg} mut*]
  (if (instance? ProteinSubstitution mut*)
    (let [pos-cands (pos-candidates (:coord mut*) rg)]
      (if (seq pos-cands)
        (let [ref-codon1 (->> pos-cands
                              (map #(cseq/read-sequence
                                     seq-rdr
                                     {:chr chr, :start %, :end %}))
                              (apply str))
              reverse? (= strand :reverse)
              ref-codon (cond-> ref-codon1 reverse? util-seq/revcomp)
              ref-aa (codon/codon->amino-acid ref-codon)
              mut-ref-aa (mut/->short-amino-acid (:ref mut*))]
          (when (= mut-ref-aa ref-aa)
            (let [palt (mut/->short-amino-acid (:alt mut*))
                  codon-cands (codon/amino-acid->codons palt)
                  pos-cands* (cond-> pos-cands reverse? reverse)]
              (->> codon-cands
                   (keep
                    (fn [codon*]
                      (when-let [[ps ref alt] (codon->variant
                                               pos-cands* ref-codon codon*)]
                        {:chr chr,
                         :pos ((if reverse? last first) ps),
                         :ref (cond-> ref reverse? util-seq/revcomp),
                         :alt (cond-> alt reverse? util-seq/revcomp)})))))))))
    (throw (ex-info "Unsupported mutation" {:type ::unsupported-mutation}))))

(defn ->vcf-variants
  [hgvs seq-rdr rg]
  (vcf-variants seq-rdr rg (:mutation hgvs)))

(defn- vcf-variants-with-coding-dna-hgvs
  [seq-rdr {:keys [chr strand] :as rg} mut*]
  (if (instance? ProteinSubstitution mut*)
    (let [pos-cands (pos-candidates (:coord mut*) rg)]
      (if (seq pos-cands)
        (let [ref-codon1 (->> pos-cands
                              (map #(cseq/read-sequence
                                     seq-rdr
                                     {:chr chr, :start %, :end %}))
                              (apply str))
              reverse? (= strand :reverse)
              ref-codon (cond-> ref-codon1 reverse? util-seq/revcomp)
              ref-aa (codon/codon->amino-acid ref-codon)
              mut-ref-aa (mut/->short-amino-acid (:ref mut*))]
          (when (= mut-ref-aa ref-aa)
            (let [palt (mut/->short-amino-acid (:alt mut*))
                  codon-cands (codon/amino-acid->codons palt)
                  pos-cands* (cond-> pos-cands reverse? reverse)]
              (->> codon-cands
                   (keep
                    (fn [codon*]
                      (when-let [[ps ref alt] (codon->variant
                                               pos-cands* ref-codon codon*)]
                        (let [cds-pos (rg/cds-coord (first ps) rg)
                              cds-end (rg/cds-coord (last ps) rg)
                              mut (if-not (= cds-pos cds-end)
                                    (mut/dna-indel cds-pos cds-end ref alt)
                                    (mut/dna-substitution cds-pos ref ">" alt))]
                          {:vcf {:chr chr,
                                 :pos ((if reverse? last first) ps),
                                 :ref (cond-> ref reverse? util-seq/revcomp),
                                 :alt (cond-> alt reverse? util-seq/revcomp)},
                           :coding-dna (hgvs/hgvs (:name rg) :coding-dna mut)}))))))))))
    (throw (IllegalArgumentException. "Unsupported mutation"))))

(defn ->vcf-variants-with-coding-dna-hgvs
  [hgvs seq-rdr rg]
  (vcf-variants-with-coding-dna-hgvs seq-rdr rg (:mutation hgvs)))
