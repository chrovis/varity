(ns varity.vcf-to-hgvs
  "Functions to convert a VCF-style variant into HGVS."
  (:require [clojure.tools.logging :as log]
            [cljam.io.sequence :as cseq]
            [cljam.io.util :as io-util]
            [cljam.util.chromosome :refer [normalize-chromosome-key]]
            [proton.string :as pstring]
            [varity.ref-gene :as rg]
            [varity.vcf-to-hgvs.genome :as genome]
            [varity.vcf-to-hgvs.cdna :as cdna]
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
      (instance? varity.ref_gene.RefGeneIndex ref-gene) :ref-gene-index
      (map? ref-gene) :ref-gene-entity)))

(defn- cdna-ref-gene? [rg]
  (some? (re-matches #"NM_\d+(\.\d+)?" (:name rg))))

(defn select-variant
  [var seq-rdr rg]
  (if-let [nvar (normalize-variant var seq-rdr rg)]
    (let [var-start-cds-coord (rg/cds-coord (:pos var) rg)
          var-end-cds-coord (rg/cds-coord (+ (:pos var) (max (count (:ref var)) (count (:alt var)))) rg)
          nvar-start-cds-coord (rg/cds-coord (:pos nvar) rg)
          nvar-end-cds-coord (rg/cds-coord (+ (:pos nvar) (max (count (:ref nvar)) (count (:alt nvar)))) rg)]
      (if (= (:region var-start-cds-coord) (:region nvar-start-cds-coord)
             (:region var-end-cds-coord) (:region nvar-end-cds-coord))
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
    (println (cdna/debug-string var seq-rdr rg))
    (when (cds-affected? var rg)
      (newline)
      (println (prot/debug-string var seq-rdr rg)))
    (newline)
    (catch Exception e
      (log/warn "Debug printing throws error" e))))

;;; -> cDNA HGVS

(defmulti vcf-variant->cdna-hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into cDNA HGVS. alt
  must be a single alternation such as \"TG\". \"TG,T\", for example, is not
  allowed. ref-seq must be a path to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of cDNA HGVS defined in clj-hgvs.

  Options:

    :tx-margin  The length of transcription margin, up to a maximum of 10000,
                default 5000.
    :verbose?   Print debug information, default false."
  {:arglists '([variant ref-seq ref-gene]
               [variant ref-seq ref-gene options])}
  (fn [_ ref-seq ref-gene & _]
    (dispatch ref-seq ref-gene)))

(defmethod vcf-variant->cdna-hgvs :ref-seq-path
  [variant ref-seq ref-gene & [options]]
  (with-open [seq-rdr (cseq/reader ref-seq)]
    (doall (vcf-variant->cdna-hgvs variant seq-rdr ref-gene options))))

(defmethod vcf-variant->cdna-hgvs :ref-gene-path
  [variant seq-rdr ref-gene & [options]]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (vcf-variant->cdna-hgvs variant seq-rdr rgidx options)))

(defmethod vcf-variant->cdna-hgvs :ref-gene-index
  [{:keys [chr pos ref alt]} seq-rdr rgidx & [options]]
  (let [{:keys [tx-margin] :or {tx-margin 5000}} options
        chr (normalize-chromosome-key chr)]
    (if (valid-ref? seq-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx tx-margin)
           (filter cdna-ref-gene?)
           (map (fn [rg]
                  (assoc (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                         seq-rdr rg)
                         :rg rg)))
           (map (fn [{:keys [rg] :as m}]
                  (when (:verbose? options)
                    (print-debug-info m seq-rdr rg))
                  (cdna/->hgvs m seq-rdr rg)))
           distinct)
      (throw (ex-info "ref is not found on the position."
                      {:type ::invalid-ref
                       :variant {:chr chr, :pos pos, :ref ref, :alt alt}})))))

(defmethod vcf-variant->cdna-hgvs :ref-gene-entity
  [{:keys [pos ref alt]} seq-rdr {:keys [chr] :as rg} & [options]]
  (if (valid-ref? seq-rdr chr pos ref)
    (let [nv (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                             seq-rdr rg)]
      (when (:verbose? options)
        (print-debug-info nv seq-rdr rg))
      (cdna/->hgvs (assoc nv :rg rg) seq-rdr rg))
    (throw (ex-info "ref is not found on the position."
                    {:type ::invalid-ref
                     :variant {:chr chr, :pos pos, :ref ref, :alt alt}}))))

;;; -> protein HGVS

(defmulti vcf-variant->protein-hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into protein HGVS.
  alt must be a single alternation such as \"TG\". \"TG,T\", for example, is not
  allowed. ref-seq must be a path to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of protein HGVS defined in clj-hgvs.

  Options:

    :verbose?  Print debug information, default false."
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
  (let [chr (normalize-chromosome-key chr)]
    (if (valid-ref? seq-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx)
           (filter cdna-ref-gene?)
           (map (fn [rg]
                  (assoc (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                         seq-rdr rg)
                         :rg rg)))
           (filter #(cds-affected? % (:rg %)))
           (keep (fn [{:keys [rg] :as m}]
                   (when (:verbose? options)
                     (print-debug-info m seq-rdr rg))
                   (prot/->hgvs m seq-rdr rg)))
           distinct)
      (throw (ex-info "ref is not found on the position."
                      {:type ::invalid-ref
                       :variant {:chr chr, :pos pos, :ref ref, :alt alt}})))))

(defmethod vcf-variant->protein-hgvs :ref-gene-entity
  [{:keys [pos ref alt]} seq-rdr {:keys [chr] :as rg} & [options]]
  (if (valid-ref? seq-rdr chr pos ref)
    (let [nv (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                             seq-rdr rg)]
      (when (cds-affected? nv rg)
        (when (:verbose? options)
          (print-debug-info nv seq-rdr rg))
        (prot/->hgvs (assoc nv :rg rg) seq-rdr rg)))
    (throw (ex-info "ref is not found on the position."
                    {:type ::invalid-ref
                     :variant {:chr chr, :pos pos, :ref ref, :alt alt}}))))

;;; -> Multiple HGVS

(defmulti vcf-variant->hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into HGVS. alt must
  be a single alternation such as \"TG\". \"TG,T\", for example, is not allowed.
  ref-seq must be a path to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of maps, each having :cdna and :protein HGVS defined in clj-hgvs.

  Options:

    :tx-margin  The length of transcription margin, up to a maximum of 10000,
                default 5000.
    :verbose?   Print debug information, default false."
  {:arglists '([variant ref-seq ref-gene]
               [variant ref-seq ref-gene options])}
  (fn [_ ref-seq ref-gene & _]
    (dispatch ref-seq ref-gene)))

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
  (let [{:keys [tx-margin] :or {tx-margin 5000}} options
        chr (normalize-chromosome-key chr)]
    (if (valid-ref? seq-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx tx-margin)
           (filter cdna-ref-gene?)
           (map (fn [rg]
                  (assoc (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                         seq-rdr rg)
                         :rg rg)))
           (map (fn [{:keys [rg] :as m}]
                  (when (:verbose? options)
                    (print-debug-info m seq-rdr rg))
                  {:cdna (cdna/->hgvs m seq-rdr rg)
                   :protein (if (cds-affected? m rg)
                              (prot/->hgvs m seq-rdr rg))}))
           distinct)
      (throw (ex-info "ref is not found on the position."
                      {:type ::invalid-ref
                       :variant {:chr chr, :pos pos, :ref ref, :alt alt}})))))

(defmethod vcf-variant->hgvs :ref-gene-entity
  [{:keys [pos ref alt]} seq-rdr {:keys [chr] :as rg} & options]
  (if (valid-ref? seq-rdr chr pos ref)
    (let [nv (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                             seq-rdr rg)]
      (when (:verbose? options)
        (print-debug-info nv seq-rdr rg))

      {:cdna (cdna/->hgvs nv seq-rdr rg)
       :protein (if (cds-affected? nv rg)
                  (prot/->hgvs nv seq-rdr rg))})
    (throw (ex-info "ref is not found on the position."
                    {:type ::invalid-ref
                     :variant {:chr chr, :pos pos, :ref ref, :alt alt}}))))
