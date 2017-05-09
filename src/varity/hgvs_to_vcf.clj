(ns varity.hgvs-to-vcf
  "Functions to convert HGVS into VCF-style variants."
  (:require [clojure.string :as string]
            [cljam.fasta :as fa]
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
  "Converts cDNA/protein hgvs into possible VCF-style variants. ref-seq of hgvs
  (e.g. NM_005228) is used for ref-genes search. Alternatively, gene (e.g. EGFR)
  can be used if ref-seq does not exist. ref-fa must be a path to reference
  FASTA or cljam.fasta.reader.FASTAReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of {:chr :pos :ref :alt}."
  {:arglists '([hgvs ref-fa ref-gene] [hgvs gene ref-fa ref-gene])}
  (fn [& args]
    (apply dispatch (take-last 2 args))))

(defmethod hgvs->vcf-variants :fasta-path
  ([hgvs ref-fa ref-gene] (hgvs->vcf-variants hgvs nil ref-fa ref-gene))
  ([hgvs gene ref-fa ref-gene]
   (with-open [fa-rdr (fa/reader ref-fa)]
     (doall (hgvs->vcf-variants hgvs gene fa-rdr ref-gene)))))

(defmethod hgvs->vcf-variants :ref-gene-path
  ([hgvs fa-rdr ref-gene] (hgvs->vcf-variants hgvs nil fa-rdr ref-gene))
  ([hgvs gene fa-rdr ref-gene]
   (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
     (hgvs->vcf-variants hgvs gene fa-rdr rgidx))))

(defmethod hgvs->vcf-variants :ref-gene-index
  ([hgvs fa-rdr rgidx] (hgvs->vcf-variants hgvs nil fa-rdr rgidx))
  ([hgvs gene fa-rdr rgidx]
   (let [convert (case (:kind hgvs)
                   :cdna cdna-hgvs->vcf-variants
                   :protein protein-hgvs->vcf-variants)
         rgs (if-let [[ref-seq] (re-find #"^(NM|NR)_\d+" (str (:transcript hgvs)))]
               (rg/ref-genes ref-seq rgidx)
               (if-not (string/blank? gene)
                 (rg/ref-genes gene rgidx)
                 (throw (IllegalArgumentException. "ref-seq (NM_, NR_) or gene must be provided"))))]
     (convert hgvs fa-rdr rgs))))

(defmethod hgvs->vcf-variants :ref-gene-entity
  ([hgvs fa-rdr rg] (hgvs->vcf-variants hgvs nil fa-rdr rg))
  ([hgvs _ fa-rdr rg]
   (let [convert (case (:kind hgvs)
                   :cdna cdna/->vcf-variant
                   :protein prot/->vcf-variants)]
     (convert hgvs fa-rdr rg))))

(defmulti protein-hgvs->vcf-variants-with-cdna-hgvs
  "Converts protein HGVS into possible VCF-style variants and cDNA HGVS. ref-seq
  of hgvs (e.g. NM_005228) is used for ref-genes search. Alternatively, gene
  (e.g. EGFR) can be used if ref-seq does not exist. ref-fa must be a path to
  reference FASTA or cljam.fasta.reader.FASTAReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of {:vcf :cdna}."
  {:arglists '([hgvs ref-fa ref-gene] [hgvs gene ref-fa ref-gene])}
  (fn [& args]
    (apply dispatch (take-last 2 args))))

(defmethod protein-hgvs->vcf-variants-with-cdna-hgvs :fasta-path
  ([hgvs ref-fa ref-gene] (protein-hgvs->vcf-variants-with-cdna-hgvs nil ref-fa ref-gene))
  ([hgvs gene ref-fa ref-gene]
   (with-open [fa-rdr (fa/reader ref-fa)]
     (doall (protein-hgvs->vcf-variants-with-cdna-hgvs hgvs gene fa-rdr ref-gene)))))

(defmethod protein-hgvs->vcf-variants-with-cdna-hgvs :ref-gene-path
  ([hgvs fa-rdr ref-gene] (protein-hgvs->vcf-variants-with-cdna-hgvs nil fa-rdr ref-gene))
  ([hgvs gene fa-rdr ref-gene]
   (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
     (protein-hgvs->vcf-variants-with-cdna-hgvs hgvs gene fa-rdr rgidx))))

(defmethod protein-hgvs->vcf-variants-with-cdna-hgvs :ref-gene-index
  ([hgvs fa-rdr rgidx] (protein-hgvs->vcf-variants-with-cdna-hgvs nil fa-rdr rgidx))
  ([hgvs gene fa-rdr rgidx]
   {:pre [(= (:kind hgvs) :protein)]}
   (->> (rg/ref-genes gene rgidx)
        (keep #(prot/->vcf-variants-with-cdna-hgvs hgvs fa-rdr %))
        (apply concat))))
