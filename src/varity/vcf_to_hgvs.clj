(ns varity.vcf-to-hgvs
  "Functions to convert a VCF-style variant into HGVS."
  (:require [cljam.io.sequence :as cseq]
            [cljam.io.util :as io-util]
            [cljam.util.chromosome :refer [normalize-chromosome-key]]
            [varity.ref-gene :as rg]
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
  (let [nvar (normalize-variant var seq-rdr rg)
        nvar-cds-coord (rg/cds-coord (+ (:pos nvar) (max (count (:ref nvar)) (count (:alt nvar)))) rg)]
    (if (nil? (:region nvar-cds-coord))
      nvar
      var)))

;;; -> cDNA HGVS

(defmulti vcf-variant->cdna-hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into cDNA HGVS. alt
  must be a single alternation such as \"TG\". \"TG,T\", for example, is not
  allowed. ref-seq must be a path to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of cDNA HGVS defined in clj-hgvs."
  {:arglists '([variant ref-seq ref-gene])}
  (fn [variant ref-seq ref-gene]
    (dispatch ref-seq ref-gene)))

(defmethod vcf-variant->cdna-hgvs :ref-seq-path
  [variant ref-seq ref-gene]
  (with-open [seq-rdr (cseq/reader ref-seq)]
    (doall (vcf-variant->cdna-hgvs variant seq-rdr ref-gene))))

(defmethod vcf-variant->cdna-hgvs :ref-gene-path
  [variant seq-rdr ref-gene]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (vcf-variant->cdna-hgvs variant seq-rdr rgidx)))

(defmethod vcf-variant->cdna-hgvs :ref-gene-index
  [{:keys [chr pos ref alt]} seq-rdr rgidx]
  (let [chr (normalize-chromosome-key chr)]
    (if (valid-ref? seq-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx)
           (filter cdna-ref-gene?)
           (map (fn [rg]
                  (assoc (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                         seq-rdr rg)
                         :rg rg)))
           (map #(cdna/->hgvs % seq-rdr (:rg %)))
           distinct)
      (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos))))))

(defmethod vcf-variant->cdna-hgvs :ref-gene-entity
  [{:keys [pos ref alt]} seq-rdr {:keys [chr] :as rg}]
  (if (valid-ref? seq-rdr chr pos ref)
    (let [nv (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                             seq-rdr rg)]
      (cdna/->hgvs (assoc nv :rg rg) seq-rdr rg))
    (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos)))))

;;; -> protein HGVS

(defmulti vcf-variant->protein-hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into protein HGVS.
  alt must be a single alternation such as \"TG\". \"TG,T\", for example, is not
  allowed. ref-seq must be a path to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of protein HGVS defined in clj-hgvs."
  {:arglists '([variant ref-seq ref-gene])}
  (fn [variant ref-seq ref-gene]
    (dispatch ref-seq ref-gene)))

(defmethod vcf-variant->protein-hgvs :ref-seq-path
  [variant ref-seq ref-gene]
  (with-open [seq-rdr (cseq/reader ref-seq)]
    (doall (vcf-variant->protein-hgvs variant seq-rdr ref-gene))))

(defmethod vcf-variant->protein-hgvs :ref-gene-path
  [variant seq-rdr ref-gene]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (vcf-variant->protein-hgvs variant seq-rdr rgidx)))

(defmethod vcf-variant->protein-hgvs :ref-gene-index
  [{:keys [chr pos ref alt]} seq-rdr rgidx]
  (let [chr (normalize-chromosome-key chr)]
    (if (valid-ref? seq-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx)
           (filter cdna-ref-gene?)
           (map (fn [rg]
                  (assoc (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                         seq-rdr rg)
                         :rg rg)))
           (filter #(rg/in-cds? (:pos %) (:rg %)))
           (keep #(prot/->hgvs % seq-rdr (:rg %)))
           distinct)
      (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos))))))

(defmethod vcf-variant->protein-hgvs :ref-gene-entity
  [{:keys [pos ref alt]} seq-rdr {:keys [chr] :as rg}]
  (if (valid-ref? seq-rdr chr pos ref)
    (let [{:keys [pos] :as nv} (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                               seq-rdr rg)]
      (if (rg/in-cds? pos rg)
        (prot/->hgvs (assoc nv :rg rg) seq-rdr rg)))
    (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos)))))

;;; -> Multiple HGVS

(defmulti vcf-variant->hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into HGVS. alt must
  be a single alternation such as \"TG\". \"TG,T\", for example, is not allowed.
  ref-seq must be a path to reference or an instance which implements
  cljam.io.protocols/ISequenceReader. ref-gene must be a path to
  refGene.txt(.gz), ref-gene index, or a ref-gene entity. A returned sequence
  consists of maps, each having :cdna and :protein HGVS defined in clj-hgvs."
  {:arglists '([variant ref-seq ref-gene])}
  (fn [variant ref-seq ref-gene]
    (dispatch ref-seq ref-gene)))

(defmethod vcf-variant->hgvs :ref-seq-path
  [variant ref-seq ref-gene]
  (with-open [seq-rdr (cseq/reader ref-seq)]
    (doall (vcf-variant->hgvs variant seq-rdr ref-gene))))

(defmethod vcf-variant->hgvs :ref-gene-path
  [variant seq-rdr ref-gene]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (vcf-variant->hgvs variant seq-rdr rgidx)))

(defmethod vcf-variant->hgvs :ref-gene-index
  [{:keys [chr pos ref alt]} seq-rdr rgidx]
  (let [chr (normalize-chromosome-key chr)]
    (if (valid-ref? seq-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx)
           (filter cdna-ref-gene?)
           (map (fn [rg]
                  (assoc (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                         seq-rdr rg)
                         :rg rg)))
           (map (fn [{:keys [pos rg] :as m}]
                  {:cdna (cdna/->hgvs m seq-rdr rg)
                   :protein (if (rg/in-cds? pos rg)
                              (prot/->hgvs m seq-rdr rg))}))
           distinct)
      (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos))))))

(defmethod vcf-variant->hgvs :ref-gene-entity
  [{:keys [pos ref alt]} seq-rdr {:keys [chr] :as rg}]
  (if (valid-ref? seq-rdr chr pos ref)
    (let [{:keys [pos] :as nv} (select-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                               seq-rdr rg)]
      {:cdna (cdna/->hgvs nv seq-rdr rg)
       :protein (if (rg/in-cds? pos rg)
                  (prot/->hgvs nv seq-rdr rg))})
    (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos)))))
