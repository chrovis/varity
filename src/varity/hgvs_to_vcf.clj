(ns varity.hgvs-to-vcf
  "Functions to convert HGVS into VCF-style variants."
  (:require [clojure.string :as string]
            [cljam.io.sequence :as cseq]
            [cljam.io.util :as io-util]
            [varity.hgvs-to-vcf.coding-dna :as coding-dna]
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

(defn- coding-dna-hgvs->vcf-variants
  [hgvs seq-rdr rgs]
  (distinct
   (keep (fn [rg]
           (try
             (coding-dna/->vcf-variant hgvs seq-rdr rg)
             (catch Exception e
               (when-not (= (:type (ex-data e)) ::rg/invalid-coordinate)
                 (throw e)))))
         rgs)))

(defn- protein-hgvs->vcf-variants
  [hgvs seq-rdr rgs]
  (distinct (apply concat (keep #(prot/->vcf-variants hgvs seq-rdr %) rgs))))

(defn- supported-transcript?
  [s]
  (some? (re-matches #"^((NM|NR)_|ENS(T|P))\d+(\.\d+)?$" (str s))))

(defmulti hgvs->vcf-variants
  "Converts coding DNA/protein hgvs into possible VCF-style variants. Transcript of
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
  ([{:keys [kind transcript] :as hgvs} gene seq-rdr rgidx]
   (let [convert (case kind
                   :coding-dna coding-dna-hgvs->vcf-variants
                   :protein protein-hgvs->vcf-variants
                   (throw (ex-info "supported HGVS kinds are only `:coding-dna` and `:protein`"
                                   {:type ::unsupported-hgvs-kind
                                    :hgvs-kind kind})))
         rgs (if (supported-transcript? transcript)
               (rg/ref-genes (str transcript))
               (if-not (string/blank? gene)
                 (rg/ref-genes gene rgidx)
                 (throw (ex-info "Transcript (NM_, NR_, ENST, ENSP) or gene must be supplied."
                                 {:type ::ref-gene-clue-not-found}))))]
     (if (seq rgs)
       (convert hgvs seq-rdr rgs)
       (throw (ex-info "Could not find specified gene in `rgidx`"
                       {:type ::gene-not-found
                        :transcript transcript :gene gene}))))))

(defmethod hgvs->vcf-variants :ref-gene-entity
  ([hgvs seq-rdr rg] (hgvs->vcf-variants hgvs nil seq-rdr rg))
  ([hgvs _ seq-rdr rg]
   (let [convert (case (:kind hgvs)
                   :coding-dna coding-dna/->vcf-variant
                   :protein prot/->vcf-variants)]
     (convert hgvs seq-rdr rg))))

(defmulti protein-hgvs->vcf-variants-with-coding-dna-hgvs
  "Converts protein HGVS into possible VCF-style variants and coding DNA HGVS.
  Transcript of hgvs, such as NM_005228, is used for ref-genes search.
  Alternatively, gene, such as EGFR, can be used if transcript does not exist.
  ref-seq must be a path to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of {:vcf :coding-dna}."
  {:arglists '([hgvs ref-seq ref-gene] [hgvs gene ref-seq ref-gene])}
  (fn [& args]
    (apply dispatch (take-last 2 args))))

(defmethod protein-hgvs->vcf-variants-with-coding-dna-hgvs :ref-seq-path
  ([hgvs ref-seq ref-gene] (protein-hgvs->vcf-variants-with-coding-dna-hgvs nil ref-seq ref-gene))
  ([hgvs gene ref-seq ref-gene]
   (with-open [seq-rdr (cseq/reader ref-seq)]
     (doall (protein-hgvs->vcf-variants-with-coding-dna-hgvs hgvs gene seq-rdr ref-gene)))))

(defmethod protein-hgvs->vcf-variants-with-coding-dna-hgvs :ref-gene-path
  ([hgvs seq-rdr ref-gene] (protein-hgvs->vcf-variants-with-coding-dna-hgvs nil seq-rdr ref-gene))
  ([hgvs gene seq-rdr ref-gene]
   (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
     (protein-hgvs->vcf-variants-with-coding-dna-hgvs hgvs gene seq-rdr rgidx))))

(defmethod protein-hgvs->vcf-variants-with-coding-dna-hgvs :ref-gene-index
  ([hgvs seq-rdr rgidx] (protein-hgvs->vcf-variants-with-coding-dna-hgvs nil seq-rdr rgidx))
  ([hgvs gene seq-rdr rgidx]
   {:pre [(= (:kind hgvs) :protein)]}
   (->> (rg/ref-genes gene rgidx)
        (keep #(prot/->vcf-variants-with-coding-dna-hgvs hgvs seq-rdr %))
        (apply concat))))
