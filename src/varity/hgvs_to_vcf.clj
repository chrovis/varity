(ns varity.hgvs-to-vcf
  "Functions to convert HGVS into VCF-style variants."
  (:require [clojure.string :as string]
            [cljam.io.sequence :as cseq]
            [cljam.io.util :as io-util]
            [varity.hgvs-to-vcf.cdna :as cdna]
            [varity.hgvs-to-vcf.protein :as prot]
            [varity.ref-gene :as rg]))

(defn- dispatch
  [ref-seq ref-gene]
  (cond
    (string? ref-seq) :ref-seq-path

    (io-util/sequence-reader? ref-seq)
    (cond
      (string? ref-gene) :ref-gene-path
      (instance? varity.ref_gene.RefGeneIndex ref-gene) :ref-gene-index
      (map? ref-gene) :ref-gene-entity)))

(defn- cdna-hgvs->vcf-variants
  [hgvs seq-rdr rgs]
  (distinct (keep #(cdna/->vcf-variant hgvs seq-rdr %) rgs)))

(defn- protein-hgvs->vcf-variants
  [hgvs seq-rdr rgs]
  (distinct (apply concat (keep #(prot/->vcf-variants hgvs seq-rdr %) rgs))))

(defmulti hgvs->vcf-variants
  "Converts cDNA/protein hgvs into possible VCF-style variants. Transcript of
  hgvs, such as NM_005228, is used for ref-genes search. Alternatively, gene,
  such as EGFR, can be used if transcript does not exist. ref-seq must be a path
  to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of {:chr :pos :ref :alt}."
  {:arglists '([hgvs ref-seq ref-gene] [hgvs gene ref-seq ref-gene])}
  (fn [& args]
    (apply dispatch (take-last 2 args))))

(defmethod hgvs->vcf-variants :ref-seq-path
  ([hgvs ref-seq ref-gene] (hgvs->vcf-variants hgvs nil ref-seq ref-gene))
  ([hgvs gene ref-seq ref-gene]
   (with-open [seq-rdr (cseq/reader ref-seq)]
     (doall (hgvs->vcf-variants hgvs gene seq-rdr ref-gene)))))

(defmethod hgvs->vcf-variants :ref-gene-path
  ([hgvs seq-rdr ref-gene] (hgvs->vcf-variants hgvs nil seq-rdr ref-gene))
  ([hgvs gene seq-rdr ref-gene]
   (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
     (hgvs->vcf-variants hgvs gene seq-rdr rgidx))))

(defmethod hgvs->vcf-variants :ref-gene-index
  ([hgvs seq-rdr rgidx] (hgvs->vcf-variants hgvs nil seq-rdr rgidx))
  ([hgvs gene seq-rdr rgidx]
   (let [convert (case (:kind hgvs)
                   :cdna cdna-hgvs->vcf-variants
                   :protein protein-hgvs->vcf-variants)
         rgs (if-let [[rs] (re-find #"^(NM|NR)_\d+" (str (:transcript hgvs)))]
               (rg/ref-genes rs rgidx)
               (if-not (string/blank? gene)
                 (rg/ref-genes gene rgidx)
                 (throw (ex-info "Transcript (NM_, NR_) or gene must be supplied."
                                 {:type ::ref-gene-clue-not-found}))))]
     (convert hgvs seq-rdr rgs))))

(defmethod hgvs->vcf-variants :ref-gene-entity
  ([hgvs seq-rdr rg] (hgvs->vcf-variants hgvs nil seq-rdr rg))
  ([hgvs _ seq-rdr rg]
   (let [convert (case (:kind hgvs)
                   :cdna cdna/->vcf-variant
                   :protein prot/->vcf-variants)]
     (convert hgvs seq-rdr rg))))

(defmulti protein-hgvs->vcf-variants-with-cdna-hgvs
  "Converts protein HGVS into possible VCF-style variants and cDNA HGVS.
  Transcript of hgvs, such as NM_005228, is used for ref-genes search.
  Alternatively, gene, such as EGFR, can be used if transcript does not exist.
  ref-seq must be a path to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of {:vcf :cdna}."
  {:arglists '([hgvs ref-seq ref-gene] [hgvs gene ref-seq ref-gene])}
  (fn [& args]
    (apply dispatch (take-last 2 args))))

(defmethod protein-hgvs->vcf-variants-with-cdna-hgvs :ref-seq-path
  ([hgvs ref-seq ref-gene] (protein-hgvs->vcf-variants-with-cdna-hgvs nil ref-seq ref-gene))
  ([hgvs gene ref-seq ref-gene]
   (with-open [seq-rdr (cseq/reader ref-seq)]
     (doall (protein-hgvs->vcf-variants-with-cdna-hgvs hgvs gene seq-rdr ref-gene)))))

(defmethod protein-hgvs->vcf-variants-with-cdna-hgvs :ref-gene-path
  ([hgvs seq-rdr ref-gene] (protein-hgvs->vcf-variants-with-cdna-hgvs nil seq-rdr ref-gene))
  ([hgvs gene seq-rdr ref-gene]
   (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
     (protein-hgvs->vcf-variants-with-cdna-hgvs hgvs gene seq-rdr rgidx))))

(defmethod protein-hgvs->vcf-variants-with-cdna-hgvs :ref-gene-index
  ([hgvs seq-rdr rgidx] (protein-hgvs->vcf-variants-with-cdna-hgvs nil seq-rdr rgidx))
  ([hgvs gene seq-rdr rgidx]
   {:pre [(= (:kind hgvs) :protein)]}
   (->> (rg/ref-genes gene rgidx)
        (keep #(prot/->vcf-variants-with-cdna-hgvs hgvs seq-rdr %))
        (apply concat))))
