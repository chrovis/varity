(ns varity.hgvs
  (:require [cljam.io.sequence :as cseq]
            [cljam.io.util :as io-util]
            [varity.hgvs-to-vcf :as h2v]
            [varity.ref-gene :as rg]
            [varity.vcf-to-hgvs :as v2h]))

(def ^:private option-patterns
  [{:prefer-deletion? false}
   {:prefer-deletion? true}
   {:prefer-insertion? false}
   {:prefer-insertion? true}
   {:prefer-deletion-insertion? false}
   {:prefer-deletion-insertion? true}
   {:prefer-extension-for-initial-codon-alt? false}
   {:prefer-extension-for-initial-codon-alt? true}])

(defn- dispatch
  [ref-seq ref-gene]
  (cond
    (string? ref-seq) :ref-seq-path

    (io-util/sequence-reader? ref-seq)
    (cond
      (string? ref-gene) :ref-gene-path
      (instance? varity.ref_gene.GeneAnnotationIndex ref-gene) :ref-gene-index
      (map? ref-gene) :ref-gene-entity)))

(defmulti find-aliases
  "Returns alternative expressions of `hgvs`, including the given expression.

  There can be multiple HGVS expressions for the same variant. For example,
  \"NM_123.4:c.4_6[1]\" and \"NM_123.4:c.7_9del\" can be equivalent in some
  cases. When you give either of the two HGVS, `find-aliases` returns both HGVS.

  `hgvs` must be a coding DNA HGVS with a transcript. `ref-seq` must be a path
  to reference or an instance which implements
  `cljam.io.protocols/ISequenceReader`. `ref-gene` must be a path to
  refGene.txt(.gz), a ref-gene index, or a ref-gene entity.

  You can select type of returned variant by variant-type option. Default is coding-dna."
  {:arglists '([hgvs ref-seq ref-gene opts])}
  (fn [_hgvs & args]
    (apply dispatch (take 2 args))))

(defmethod find-aliases :ref-seq-path
  [hgvs ref-seq ref-gene & {:as opts}]
  (with-open [seq-rdr (cseq/reader ref-seq)]
    (doall (find-aliases hgvs seq-rdr ref-gene opts))))

(defmethod find-aliases :ref-gene-path
  [hgvs seq-rdr ref-gene & {:as opts}]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (find-aliases hgvs seq-rdr rgidx opts)))

(defmethod find-aliases :ref-gene-index
  [hgvs seq-rdr rgidx & {:keys [variant-type] :or {variant-type :coding-dna}}]
  {:pre [(= (:kind hgvs) :coding-dna)]}
  (let [rgidx* (rg/index (rg/ref-genes (:transcript hgvs) rgidx))
        variants (h2v/hgvs->vcf-variants hgvs seq-rdr rgidx*)
        vcf-variant->hgvs* (case variant-type
                             :coding-dna v2h/vcf-variant->coding-dna-hgvs
                             :protein v2h/vcf-variant->protein-hgvs)
        vcf-variant->hgvs (fn [variant]
                            (mapcat #(vcf-variant->hgvs* variant seq-rdr rgidx* %)
                                    option-patterns))]
    (if (seq variants)
      (->> variants
           (mapcat vcf-variant->hgvs)
           distinct)
      (throw (ex-info "The VCF variant is not found."
                      {:type ::invalid-variant
                       :hgvs hgvs})))))

(defmethod find-aliases :ref-gene-entity
  [hgvs seq-rdr rg & {:keys [variant-type] :or {variant-type :coding-dna}}]
  {:pre [(= (:kind hgvs) :coding-dna)]}
  (if-let [variant (h2v/hgvs->vcf-variants hgvs seq-rdr rg)]
    (let [vcf-variant->hgvs (case variant-type
                              :coding-dna v2h/vcf-variant->coding-dna-hgvs
                              :protein v2h/vcf-variant->protein-hgvs)]
      (->> option-patterns
           (map #(vcf-variant->hgvs variant seq-rdr rg %))
           distinct))
    (throw (ex-info "The VCF variant is not found."
                    {:type ::invalid-variant
                     :hgvs hgvs}))))
