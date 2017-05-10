(ns varity.vcf-to-hgvs
  "Functions to convert a VCF-style variant into HGVS."
  (:require [cljam.fasta :as fa]
            [cljam.util.chromosome :refer [normalize-chromosome-key]]
            [varity.ref-gene :as rg]
            [varity.vcf-to-hgvs.cdna :as cdna]
            [varity.vcf-to-hgvs.common :refer [normalize-variant]]
            [varity.vcf-to-hgvs.protein :as prot]))

(defn- valid-ref?
  [fa-rdr chr pos ref]
  (= (fa/read-sequence fa-rdr {:chr chr, :start pos, :end (+ pos (count ref) -1)}) ref))

(defn- dispatch
  [ref-fa ref-gene]
  (cond
    (string? ref-fa) :fasta-path

    (instance? cljam.fasta.reader.FASTAReader ref-fa)
    (cond
      (string? ref-gene) :ref-gene-path
      (instance? varity.ref_gene.RefGeneIndex ref-gene) :ref-gene-index
      (map? ref-gene) :ref-gene-entity)))

;;; -> cDNA HGVS

(defmulti vcf-variant->cdna-hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into cDNA HGVS. alt
  must be a single alternation such as \"TG\". \"TG,T\", for example, is not
  allowed. ref-fa must be a path to reference FASTA or
  cljam.fasta.reader.FASTAReader. ref-gene must be a path to refGene.txt(.gz),
  ref-gene index, or a ref-gene entity. A returned sequence consists of cDNA
  HGVS defined in clj-hgvs."
  {:arglists '([variant ref-fa ref-gene])}
  (fn [variant ref-fa ref-gene]
    (dispatch ref-fa ref-gene)))

(defmethod vcf-variant->cdna-hgvs :fasta-path
  [variant ref-fa ref-gene]
  (with-open [fa-rdr (fa/reader ref-fa)]
    (doall (vcf-variant->cdna-hgvs variant fa-rdr ref-gene))))

(defmethod vcf-variant->cdna-hgvs :ref-gene-path
  [variant fa-rdr ref-gene]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (vcf-variant->cdna-hgvs variant fa-rdr rgidx)))

(defmethod vcf-variant->cdna-hgvs :ref-gene-index
  [{:keys [chr pos ref alt]} fa-rdr rgidx]
  (let [chr (normalize-chromosome-key chr)]
    (if (valid-ref? fa-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx)
           (map (fn [rg]
                  (assoc (normalize-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                            fa-rdr rg)
                         :rg rg)))
           (map #(cdna/->hgvs % fa-rdr (:rg %)))
           distinct)
      (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos))))))

(defmethod vcf-variant->cdna-hgvs :ref-gene-entity
  [{:keys [pos ref alt]} fa-rdr {:keys [chr] :as rg}]
  (if (valid-ref? fa-rdr chr pos ref)
    (let [nv (normalize-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                fa-rdr rg)]
      (cdna/->hgvs (assoc nv :rg rg) fa-rdr rg))
    (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos)))))

;;; -> protein HGVS

(defmulti vcf-variant->protein-hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into protein HGVS.
  alt must be a single alternation such as \"TG\". \"TG,T\", for example, is not
  allowed. ref-fa must be a path to reference FASTA or
  cljam.fasta.reader.FASTAReader. ref-gene must be a path to refGene.txt(.gz),
  ref-gene index, or a ref-gene entity. A returned sequence consists of protein
  HGVS defined in clj-hgvs."
  {:arglists '([variant ref-fa ref-gene])}
  (fn [variant ref-fa ref-gene]
    (dispatch ref-fa ref-gene)))

(defmethod vcf-variant->protein-hgvs :fasta-path
  [variant ref-fa ref-gene]
  (with-open [fa-rdr (fa/reader ref-fa)]
    (doall (vcf-variant->protein-hgvs variant fa-rdr ref-gene))))

(defmethod vcf-variant->protein-hgvs :ref-gene-path
  [variant fa-rdr ref-gene]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (vcf-variant->cdna-hgvs variant fa-rdr rgidx)))

(defmethod vcf-variant->protein-hgvs :ref-gene-index
  [{:keys [chr pos ref alt]} fa-rdr rgidx]
  (let [chr (normalize-chromosome-key chr)]
    (if (valid-ref? fa-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx)
           (map (fn [rg]
                  (assoc (normalize-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                            fa-rdr rg)
                         :rg rg)))
           (filter #(rg/in-cds? (:pos %) (:rg %)))
           (keep #(prot/->hgvs % fa-rdr (:rg %)))
           distinct)
      (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos))))))

(defmethod vcf-variant->protein-hgvs :ref-gene-entity
  [{:keys [pos ref alt]} fa-rdr {:keys [chr] :as rg}]
  (if (valid-ref? fa-rdr chr pos ref)
    (let [{:keys [pos] :as nv} (normalize-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                                   fa-rdr rg)]
      (if (rg/in-cds? pos rg)
        (prot/->hgvs (assoc nv :rg rg) fa-rdr rg)))
    (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos)))))

;;; -> Multiple HGVS

(defmulti vcf-variant->hgvs
  "Converts a VCF-style variant (:chr, :pos, :ref, and :alt) into HGVS. alt must
  be a single alternation such as \"TG\". \"TG,T\", for example, is not allowed.
  ref-fa must be a path to reference FASTA or cljam.fasta.reader.FASTAReader.
  ref-gene must be a path to refGene.txt(.gz), ref-gene index, or a ref-gene
  entity. A returned sequence consists of maps, each having :cdna and :protein
  HGVS defined in clj-hgvs."
  {:arglists '([variant ref-fa ref-gene])}
  (fn [variant ref-fa ref-gene]
    (dispatch ref-fa ref-gene)))

(defmethod vcf-variant->hgvs :fasta-path
  [variant ref-fa ref-gene]
  (with-open [fa-rdr (fa/reader ref-fa)]
    (doall (vcf-variant->hgvs variant fa-rdr ref-gene))))

(defmethod vcf-variant->hgvs :ref-gene-path
  [variant fa-rdr ref-gene]
  (let [rgidx (rg/index (rg/load-ref-genes ref-gene))]
    (vcf-variant->hgvs variant fa-rdr rgidx)))

(defmethod vcf-variant->hgvs :ref-gene-index
  [{:keys [chr pos ref alt]} fa-rdr rgidx]
  (let [chr (normalize-chromosome-key chr)]
    (if (valid-ref? fa-rdr chr pos ref)
      (->> (rg/ref-genes chr pos rgidx)
           (map (fn [rg]
                  (assoc (normalize-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                            fa-rdr rg)
                         :rg rg)))
           (map (fn [{:keys [pos rg] :as m}]
                  {:cdna (cdna/->hgvs m fa-rdr rg)
                   :protein (if (rg/in-cds? pos rg)
                              (prot/->hgvs m fa-rdr rg))}))
           distinct)
      (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos))))))

(defmethod vcf-variant->hgvs :ref-gene-entity
  [{:keys [pos ref alt]} fa-rdr {:keys [chr] :as rg}]
  (if (valid-ref? fa-rdr chr pos ref)
    (let [{:keys [pos] :as nv} (normalize-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                                  fa-rdr rg)]
      {:cdna (cdna/->hgvs nv fa-rdr rg)
       :protein (if (rg/in-cds? pos rg)
                  (prot/->hgvs nv fa-rdr rg))})
    (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos)))))
