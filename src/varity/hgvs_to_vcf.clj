(ns varity.hgvs-to-vcf
  "Functions to convert HGVS into VCF-style variants."
  (:require [cljam.fasta :as fa]
            [varity.hgvs-to-vcf.cdna :as cdna]
            [varity.hgvs-to-vcf.protein :as prot]
            [varity.ref-gene :as rg]))

(defn- dispatch
  [ref-fa ref-gene]
  (cond
    (string? ref-fa) :fasta-path

    (instance? cljam.fasta.reader.FASTAReader ref-fa)
    (cond
      (string? ref-gene) :ref-gene-path
      (instance? varity.ref_gene.RefGeneIndex ref-gene) :ref-gene-index
      (map? ref-gene) :ref-gene-entity)))

(defn- cdna-hgvs->vcf-variants
  [hgvs ref-fa rgs]
  (distinct (keep #(cdna/->vcf-variant hgvs ref-fa %) rgs)))

(defn- protein-hgvs->vcf-variants
  [hgvs ref-fa rgs]
  (distinct (apply concat (keep #(prot/->vcf-variants hgvs ref-fa %) rgs))))

(defmulti hgvs->vcf-variants
  "Converts cDNA/protein HGVS into possible VCF-style variants. ref-fa must be a
  path to reference FASTA or cljam.fasta.reader.FASTAReader. ref-gene must be a
  path to refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned
  sequence consists of {:chr :pos :ref :alt}."
  {:arglists '([hgvs gene ref-fa ref-gene])}
  (fn [hgvs gene ref-fa ref-gene]
    (dispatch ref-fa ref-gene)))

(defmethod hgvs->vcf-variants :fasta-path
  [hgvs gene ref-fa ref-gene]
  (with-open [fa-rdr (fa/reader ref-fa)]
    (doall (hgvs->vcf-variants hgvs gene fa-rdr ref-gene))))

(defmethod hgvs->vcf-variants :ref-gene-path
  [hgvs gene fa-rdr ref-gene]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (hgvs->vcf-variants hgvs gene fa-rdr rgidx)))

(defmethod hgvs->vcf-variants :ref-gene-index
  [hgvs gene fa-rdr rgidx]
  (let [convert (case (:kind hgvs)
                  :cdna cdna-hgvs->vcf-variants
                  :protein protein-hgvs->vcf-variants)]
    (convert hgvs fa-rdr (rg/ref-genes gene rgidx))))

(defmethod hgvs->vcf-variants :ref-gene-entity
  [hgvs gene fa-rdr rg]
  (let [convert (case (:kind hgvs)
                  :cdna cdna/->vcf-variant
                  :protein prot/->vcf-variants)]
    (convert hgvs fa-rdr rg)))

(defmulti protein-hgvs->vcf-variants-with-cdna-hgvs
  "Converts protein HGVS into possible VCF-style variants and cDNA HGVS. ref-fa
  must be a path to reference FASTA or cljam.fasta.reader.FASTAReader. ref-gene
  must be a path to refGene.txt(.gz), ref-gene index, or a ref-gene entity. A
  returned sequence consists of {:vcf :cdna}."
  {:arglists '([hgvs gene ref-fa ref-gene])}
  (fn [hgvs gene ref-fa ref-gene]
    (dispatch ref-fa ref-gene)))

(defmethod protein-hgvs->vcf-variants-with-cdna-hgvs :fasta-path
  [hgvs gene ref-fa ref-gene]
  (with-open [fa-rdr (fa/reader ref-fa)]
    (doall (protein-hgvs->vcf-variants-with-cdna-hgvs hgvs gene fa-rdr ref-gene))))

(defmethod protein-hgvs->vcf-variants-with-cdna-hgvs :ref-gene-path
  [hgvs gene fa-rdr ref-gene]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (protein-hgvs->vcf-variants-with-cdna-hgvs hgvs gene fa-rdr rgidx)))

(defmethod protein-hgvs->vcf-variants-with-cdna-hgvs :ref-gene-index
  [hgvs gene fa-rdr rgidx]
  {:pre [(= (:kind hgvs) :protein)]}
  (->> (rg/ref-genes gene rgidx)
       (keep #(prot/->vcf-variants-with-cdna-hgvs hgvs fa-rdr %))
       (apply concat)))
