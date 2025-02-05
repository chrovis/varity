(ns varity.vcf-to-hgvs
  "Functions to convert a VCF-style variant into HGVS."
  (:require [clojure.tools.logging :as log]
            [cljam.io.sequence :as cseq]
            [cljam.io.util :as io-util]
            [cljam.util.chromosome :refer [normalize-chromosome-key]]
            [proton.string :as pstring]
            [varity.ref-gene :as rg]
            [varity.vcf-to-hgvs.genome :as genome]
            [varity.vcf-to-hgvs.coding-dna :as coding-dna]
            [varity.vcf-to-hgvs.common :refer [normalize-variant]]
            [varity.vcf-to-hgvs.protein :as prot]))

(defn- valid-ref?
  [seq-rdr chr pos ref]
  (= (cseq/read-sequence seq-rdr {:chr chr, :start pos, :end (+ pos (count ref) -1)}) ref))

(defn- dispatch
  [ref-seq ref-gene]
  (cond
    (string? ref-seq) :ref-seq-path

    (io-util/sequence-reader? ref-seq)
    (cond
      (string? ref-gene) :ref-gene-path
      (instance? varity.ref_gene.GeneAnnotationIndex ref-gene) :ref-gene-index
      (map? ref-gene) :ref-gene-entity)))

(defn- coding-dna-ref-gene? [rg]
  (some? (re-matches #"(NM_|ENST)\d+(\.\d+)?" (:name rg))))

(defn select-variant
  [var seq-rdr rg & {:keys [three-prime-rule]}]
  (if-let [nvar (normalize-variant var seq-rdr rg)]
    (let [var-start-cds-coord (rg/cds-coord (:pos var) rg)
          var-end-cds-coord (rg/cds-coord (+ (:pos var) (max (count (:ref var)) (count (:alt var)))) rg)
          nvar-start-cds-coord (rg/cds-coord (:pos nvar) rg)
          nvar-end-cds-coord (rg/cds-coord (+ (:pos nvar) (max (count (:ref nvar)) (count (:alt nvar)))) rg)
          restrict-cds (:restrict-cds three-prime-rule)]
      (if (or (= (:region var-start-cds-coord) (:region nvar-start-cds-coord)
                  (:region var-end-cds-coord) (:region nvar-end-cds-coord))
              (not restrict-cds))
        nvar
        var))
    var))

(defn- cds-affected?
  [var rg]
  (and (<= (:cds-start rg) (+ (:pos var) (count (:ref var))))
       (<= (:pos var) (:cds-end rg))))

(defn- print-debug-info
  [var seq-rdr rg]
  (try
    (println " variant:" (-> (select-keys var [:chr :pos :ref :alt])
                             (update :ref pstring/prune 20)
                             (update :alt pstring/prune 20)
                             str))
    (println "ref-gene:" (str (select-keys rg [:name :name2 :strand])))
    (newline)
    (println (genome/debug-string var seq-rdr rg))
    (newline)
    (println (coding-dna/debug-string var seq-rdr rg))
    (when (cds-affected? var rg)
      (newline)
      (println (prot/debug-string var seq-rdr rg)))
    (newline)
    (catch Exception e
      (log/warn "Debug printing throws error" e))))

(def ^:private default-options
  {:prefer-deletion? false
   :prefer-insertion? false
   :prefer-extension-for-initial-codon-alt? false
   :tx-margin 5000
   :verbose? false
   :three-prime-rule {:restrict-cds true}})

;;; -> Coding DNA HGVS

(defmulti vcf-variant->coding-dna-hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into coding DNA HGVS.
  alt must be a single alternation such as \"TG\". \"TG,T\", for example, is not
  allowed. ref-seq must be a path to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of coding DNA HGVS defined in clj-hgvs.

  Options:

    :prefer-deletion?  Prefer deletion (e.g. \"c.7_9del\") to repeated
                       sequences (e.g. \"c.4_6[1]\"), default false.

    :prefer-insertion? Prefer insertion (e.g. \"c.9_10insAGG\") to repeated
                       sequences (e.g. \"c.4_6[3]\"), default false.

    :tx-margin         The length of transcription margin, up to a maximum of
                       10000, default 5000.

    :three-prime-rule  Three prime rule applied condition
                       default {:restrict-cds true}.

    :verbose?          Print debug information, default false."
  {:arglists '([variant ref-seq ref-gene]
               [variant ref-seq ref-gene options])}
  (fn [_ ref-seq ref-gene & _]
    (dispatch ref-seq ref-gene)))

(defmethod vcf-variant->coding-dna-hgvs :ref-seq-path
  [variant ref-seq ref-gene & [options]]
  (with-open [seq-rdr (cseq/reader ref-seq)]
    (doall (vcf-variant->coding-dna-hgvs variant seq-rdr ref-gene options))))

(defmethod vcf-variant->coding-dna-hgvs :ref-gene-path
  [variant seq-rdr ref-gene & [options]]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (vcf-variant->coding-dna-hgvs variant seq-rdr rgidx options)))

(defmethod vcf-variant->coding-dna-hgvs :ref-gene-index
  [{:keys [chr pos ref alt]} seq-rdr rgidx & [options]]
  (let [options (merge default-options options)
        chr (normalize-chromosome-key chr)]
    (if (valid-ref? seq-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx (:tx-margin options))
           (filter coding-dna-ref-gene?)
           (map (fn [rg]
                  (assoc (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                         seq-rdr rg
                                         :three-prime-rule (:three-prime-rule options))
                         :rg rg)))
           (map (fn [{:keys [rg] :as m}]
                  (when (:verbose? options)
                    (print-debug-info m seq-rdr rg))
                  (coding-dna/->hgvs m seq-rdr rg options)))
           distinct)
      (throw (ex-info "ref is not found on the position."
                      {:type ::invalid-ref
                       :variant {:chr chr, :pos pos, :ref ref, :alt alt}})))))

(defmethod vcf-variant->coding-dna-hgvs :ref-gene-entity
  [{:keys [pos ref alt]} seq-rdr {:keys [chr] :as rg} & [options]]
  (let [options (merge default-options options)]
    (if (valid-ref? seq-rdr chr pos ref)
      (let [nv (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                               seq-rdr rg
                               :three-prime-rule (:three-prime-rule options))]
        (when (:verbose? options)
          (print-debug-info nv seq-rdr rg))
        (coding-dna/->hgvs (assoc nv :rg rg) seq-rdr rg options))
      (throw (ex-info "ref is not found on the position."
                      {:type ::invalid-ref
                       :variant {:chr chr, :pos pos, :ref ref, :alt alt}})))))

;;; -> protein HGVS

(defmulti vcf-variant->protein-hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into protein HGVS.
  alt must be a single alternation such as \"TG\". \"TG,T\", for example, is not
  allowed. ref-seq must be a path to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of protein HGVS defined in clj-hgvs.

  Options:

    :prefer-deletion?  Prefer deletion (e.g. \"p.P7_H8del\") to repeated
                       sequences (e.g. \"p.P5_H6[1]\"), default false.

    :prefer-insertion? Prefer insertion (e.g. \"c.H9_L10insRPH\") to repeated
                       sequences (e.g. \"c.R4_H6[3]\"), default false.

    :prefer-extension-for-initial-codon-alt? Prefer extension to protein unknown variant
                                             that affects initial codon, default false.

    :three-prime-rule  Three prime rule applied condition
                       default {:restrict-cds true}.

    :verbose?          Print debug information, default false."
  {:arglists '([variant ref-seq ref-gene]
               [variant ref-seq ref-gene options])}
  (fn [_ ref-seq ref-gene & _]
    (dispatch ref-seq ref-gene)))

(defmethod vcf-variant->protein-hgvs :ref-seq-path
  [variant ref-seq ref-gene & [options]]
  (with-open [seq-rdr (cseq/reader ref-seq)]
    (doall (vcf-variant->protein-hgvs variant seq-rdr ref-gene options))))

(defmethod vcf-variant->protein-hgvs :ref-gene-path
  [variant seq-rdr ref-gene & [options]]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (vcf-variant->protein-hgvs variant seq-rdr rgidx options)))

(defmethod vcf-variant->protein-hgvs :ref-gene-index
  [{:keys [chr pos ref alt]} seq-rdr rgidx & [options]]
  (let [options (merge default-options options)
        chr (normalize-chromosome-key chr)]
    (if (valid-ref? seq-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx)
           (filter coding-dna-ref-gene?)
           (map (fn [rg]
                  (assoc (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                         seq-rdr rg
                                         :three-prime-rule (:three-prime-rule options))
                         :rg rg)))
           (filter #(cds-affected? % (:rg %)))
           (keep (fn [{:keys [rg] :as m}]
                   (when (:verbose? options)
                     (print-debug-info m seq-rdr rg))
                   (prot/->hgvs m seq-rdr rg options)))
           distinct)
      (throw (ex-info "ref is not found on the position."
                      {:type ::invalid-ref
                       :variant {:chr chr, :pos pos, :ref ref, :alt alt}})))))

(defmethod vcf-variant->protein-hgvs :ref-gene-entity
  [{:keys [pos ref alt]} seq-rdr {:keys [chr] :as rg} & [options]]
  (let [options (merge default-options options)]
    (if (valid-ref? seq-rdr chr pos ref)
      (let [nv (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                               seq-rdr rg
                               :three-prime-rule (:three-prime-rule options))]
        (when (cds-affected? nv rg)
          (when (:verbose? options)
            (print-debug-info nv seq-rdr rg))
          (prot/->hgvs (assoc nv :rg rg) seq-rdr rg options)))
      (throw (ex-info "ref is not found on the position."
                      {:type ::invalid-ref
                       :variant {:chr chr, :pos pos, :ref ref, :alt alt}})))))

;;; -> Multiple HGVS

(defmulti vcf-variant->hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into HGVS. alt must
  be a single alternation such as \"TG\". \"TG,T\", for example, is not allowed.
  ref-seq must be a path to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of maps, each having :coding-dna and :protein HGVS defined in
  clj-hgvs.

  Options:

    :prefer-deletion?  Prefer deletion (e.g. \"c.7_9del\") to repeated
                       sequences (e.g. \"c.4_6[1]\"), default false.

    :prefer-insertion? Prefer insertion (e.g. \"c.9_10insAGG\") to repeated
                       sequences (e.g. \"c.4_6[3]\"), default false.

    :prefer-extension-for-initial-codon-alt? Prefer extension to protein unknown variant
                                             that affects initial codon, default false.

    :tx-margin         The length of transcription margin, up to a maximum of
                       10000, default 5000.

    :three-prime-rule  Three prime rule applied condition
                       default {:coding-dna {:restrict-cds true} :protein {:restrict-cds true}}.

    :verbose?          Print debug information, default false."
  {:arglists '([variant ref-seq ref-gene]
               [variant ref-seq ref-gene options])}
  (fn [_ ref-seq ref-gene & _]
    (dispatch ref-seq ref-gene)))

(def ^:private vcf-variant->hgvs-three-prime-rule-default-option
  {:three-prime-rule {:coding-dna {:restrict-cds true}
                      :protein {:restrict-cds true}}})

(defmethod vcf-variant->hgvs :ref-seq-path
  [variant ref-seq ref-gene & [options]]
  (with-open [seq-rdr (cseq/reader ref-seq)]
    (doall (vcf-variant->hgvs variant seq-rdr ref-gene options))))

(defmethod vcf-variant->hgvs :ref-gene-path
  [variant seq-rdr ref-gene & [options]]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (vcf-variant->hgvs variant seq-rdr rgidx options)))

(defmethod vcf-variant->hgvs :ref-gene-index
  [{:keys [chr pos ref alt]} seq-rdr rgidx & [options]]
  (let [options (merge default-options vcf-variant->hgvs-three-prime-rule-default-option options)
        chr (normalize-chromosome-key chr)]
    (if (valid-ref? seq-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx (:tx-margin options))
           (filter coding-dna-ref-gene?)
           (map (fn [rg]
                  (let [coding-dna-three-prime-rule (get-in options [:three-prime-rule :coding-dna])
                        protein-three-prime-rule (get-in options [:three-prime-rule :protein])
                        nv (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                           seq-rdr rg
                                           :three-prime-rule coding-dna-three-prime-rule)]
                    {:coding-dna (assoc nv
                                        :rg rg
                                        :type :coding-dna)
                     :protein (assoc (if (= coding-dna-three-prime-rule protein-three-prime-rule)
                                       nv
                                       (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                                       seq-rdr rg
                                                       :three-prime-rule protein-three-prime-rule))
                                     :rg rg
                                     :type :protein)})))
           (map (fn [{:keys [coding-dna protein]}]
                  (when (:verbose? options)
                    (print-debug-info coding-dna seq-rdr (:rg coding-dna))
                    (print-debug-info protein seq-rdr (:rg protein)))
                  {:coding-dna (coding-dna/->hgvs coding-dna seq-rdr (:rg coding-dna) options)
                   :protein (when (cds-affected? coding-dna (:rg coding-dna))
                                (prot/->hgvs protein seq-rdr (:rg protein) options))}))
           distinct)
      (throw (ex-info "ref is not found on the position."
                      {:type ::invalid-ref
                       :variant {:chr chr, :pos pos, :ref ref, :alt alt}})))))

(defmethod vcf-variant->hgvs :ref-gene-entity
  [{:keys [pos ref alt]} seq-rdr {:keys [chr] :as rg} & [options]]
  (let [options (merge default-options vcf-variant->hgvs-three-prime-rule-default-option options)]
    (if (valid-ref? seq-rdr chr pos ref)
      (let [coding-dna-three-prime-rule (get-in options [:three-prime-rule :coding-dna])
            protein-three-prime-rule (get-in options [:three-prime-rule :protein])
            nv (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                               seq-rdr rg
                               :three-prime-rule coding-dna-three-prime-rule)
            dnv (assoc nv
                       :type :coding-dna)
            pnv (assoc (if (= coding-dna-three-prime-rule protein-three-prime-rule)
                         nv
                         (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                         seq-rdr rg
                                         :three-prime-rule protein-three-prime-rule))
                       :type :protein)]
        (when (:verbose? options)
          (print-debug-info dnv seq-rdr rg)
          (print-debug-info pnv seq-rdr rg))
        (println (:three-prime-rule options))
        {:coding-dna (coding-dna/->hgvs dnv seq-rdr rg options)
         :protein (when (cds-affected? dnv rg)
                    (prot/->hgvs pnv seq-rdr rg options))})
      (throw (ex-info "ref is not found on the position."
                      {:type ::invalid-ref
                       :variant {:chr chr, :pos pos, :ref ref, :alt alt}})))))
