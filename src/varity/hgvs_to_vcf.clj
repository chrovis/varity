(ns varity.hgvs-to-vcf
  "Functions to convert HGVS into VCF-style variants."
  (:require [cljam.fasta :as fa]
            [varity.hgvs-to-vcf.cdna :as cdna]
            [varity.hgvs-to-vcf.protein :as prot]
            [varity.ref-gene :as rg]))

(defn- cdna-hgvs->vcf-variants
  [ref-fa rgs hgvs]
  (keep #(cdna/rg->vcf-variant ref-fa % hgvs) rgs))

(defn- protein-hgvs->vcf-variants
  [ref-fa rgs hgvs]
  (apply concat (keep #(prot/rg->vcf-variants ref-fa % hgvs) rgs)))

(defmulti hgvs->vcf-variants
  "Converts cDNA/protein HGVS into possible VCF-style variants. The first
  argument, ref-fa, must be a path to reference FASTA or
  cljam.fasta.reader.FASTAReader. The second argument, ref-gene, must be a path
  to refGene.txt(.gz) or ref-gene index. A returned sequence consists of
  {:chr :pos :ref :alt}."
  {:arglists '([ref-fa ref-gene hgvs gene])}
  (fn [ref-fa ref-gene hgvs gene]
    (class ref-fa)))

(defmethod hgvs->vcf-variants String
  [ref-fa ref-gene hgvs gene]
  (with-open [fa-rdr (fa/reader ref-fa)]
    (doall (hgvs->vcf-variants fa-rdr ref-gene hgvs gene))))

(defmethod hgvs->vcf-variants cljam.fasta.reader.FASTAReader
  [ref-fa ref-gene hgvs gene]
  (let [rgidx (if (map? ref-gene)
                ref-gene
                (rg/index (rg/load-ref-genes ref-gene)))
        convert (case (:kind hgvs)
                  :cdna cdna-hgvs->vcf-variants
                  :protein protein-hgvs->vcf-variants)]
    (convert ref-fa (rg/ref-genes gene rgidx) hgvs)))

(defmulti protein-hgvs->vcf-variants-with-cdna-hgvs
  "Converts protein HGVS into possible VCF-style variants and cDNA HGVS. The
  first argument, ref-fa, must be a path to reference FASTA or
  cljam.fasta.reader.FASTAReader. The second argument, ref-gene, must be a path
  to refGene.txt(.gz) or ref-gene index. A returned sequence consists of
  {:vcf :cdna}."
  {:arglists '([ref-fa ref-gene hgvs gene])}
  (fn [ref-fa ref-gene hgvs gene]
    (class ref-fa)))

(defmethod protein-hgvs->vcf-variants-with-cdna-hgvs String
  [ref-fa ref-gene hgvs gene]
  (with-open [fa-rdr (fa/reader ref-fa)]
    (doall (protein-hgvs->vcf-variants-with-cdna-hgvs fa-rdr ref-gene hgvs gene))))

(defmethod protein-hgvs->vcf-variants-with-cdna-hgvs cljam.fasta.reader.FASTAReader
  [ref-fa ref-gene hgvs gene]
  {:pre [(= (:kind hgvs) :protein)]}
  (let [rgidx (if (map? ref-gene)
                ref-gene
                (rg/index (rg/load-ref-genes ref-gene)))]
    (->> (rg/ref-genes gene rgidx)
         (keep #(prot/rg->vcf-variants-with-cdna-hgvs ref-fa % hgvs))
         (apply concat))))
