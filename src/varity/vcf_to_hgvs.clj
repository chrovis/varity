(ns varity.vcf-to-hgvs
  "Functions to convert a VCF-style variant into HGVS."
  (:require [cljam.fasta :as fa]
            [varity.ref-gene :as rg]
            [varity.vcf-to-hgvs.cdna :as cdna]
            [varity.vcf-to-hgvs.common :refer [normalize-variant]]
            [varity.vcf-to-hgvs.protein :as prot]))

(defn- valid-ref?
  [fa-rdr chr pos ref]
  (= (fa/read-sequence fa-rdr {:chr chr, :start pos, :end (+ pos (count ref) -1)}) ref))

;;; -> cDNA HGVS

(defmulti vcf-variant->cdna-hgvs
  "Converts a VCF-style variant (chr, pos, ref, and alt) into cDNA HGVS. The
  first argument, ref-fa, must be a path to reference FASTA or
  cljam.fasta.reader.FASTAReader. The second argument, ref-gene, must be a path
  to refGene.txt(.gz) or ref-gene index. alt must be a single alternation such
  as \"TG\". \"TG,T\", for example, is not allowed. A returned sequence consists
  of cDNA HGVS defined in clj-hgvs."
  {:arglists '([ref-fa ref-gene chr pos ref alt])}
  (fn [ref-fa ref-gene chr pos ref alt]
    (class ref-fa)))

(defmethod vcf-variant->cdna-hgvs String
  [ref-fa ref-gene chr pos ref alt]
  (with-open [fa-rdr (fa/reader ref-fa)]
    (doall (vcf-variant->cdna-hgvs fa-rdr ref-gene chr pos ref alt))))

(defmethod vcf-variant->cdna-hgvs cljam.fasta.reader.FASTAReader
  [ref-fa ref-gene chr pos ref alt]
  (if (valid-ref? ref-fa chr pos ref)
    (let [rgidx (if (map? ref-gene)
                  ref-gene
                  (rg/index (rg/load-ref-genes ref-gene)))]
      (->> (rg/ref-genes chr pos rgidx)
           (map (fn [rg]
                  (assoc (normalize-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                            ref-fa rg)
                         :rg rg)))
           (map (fn [{:keys [rg chr pos ref alt]}]
                  (cdna/rg->hgvs ref-fa rg pos ref alt)))
           distinct))
    (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos)))))

;;; -> protein HGVS

(defmulti vcf-variant->protein-hgvs
  "Converts a VCF-style variant (chr, pos, ref, and alt) into protein HGVS. The
  first argument, ref-fa, must be a path to reference FASTA or
  cljam.fasta.reader.FASTAReader. The second argument, ref-gene, must be a path
  to refGene.txt(.gz) or ref-gene index. alt must be a single alternation such
  as \"TG\". \"TG,T\", for example, is not allowed. A returned sequence consists
  of protein HGVS defined in clj-hgvs."
  {:arglists '([ref-fa ref-gene chr pos ref alt])}
  (fn [ref-fa ref-gene chr pos ref alt]
    (class ref-fa)))

(defmethod vcf-variant->protein-hgvs String
  [ref-fa ref-gene chr pos ref alt]
  (with-open [fa-rdr (fa/reader ref-fa)]
    (doall (vcf-variant->protein-hgvs fa-rdr ref-gene chr pos ref alt))))

(defmethod vcf-variant->protein-hgvs cljam.fasta.reader.FASTAReader
  [ref-fa ref-gene chr pos ref alt]
  (if (valid-ref? ref-fa chr pos ref)
    (let [rgidx (if (map? ref-gene)
                  ref-gene
                  (rg/index (rg/load-ref-genes ref-gene)))]
      (->> (rg/ref-genes chr pos rgidx)
           (map (fn [rg]
                  (assoc (normalize-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                            ref-fa rg)
                         :rg rg)))
           (filter #(rg/in-cds? pos (:rg %)))
           (keep (fn [{:keys [rg chr pos ref alt]}]
                   (prot/rg->hgvs ref-fa rg pos ref alt)))
           distinct))
    (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos)))))

;;; -> Multiple HGVS

(defmulti vcf-variant->hgvs
  "Converts a VCF-style variant (chr, pos, ref, and alt) into HGVS. The first
  argument, ref-fa, must be a path to reference FASTA or
  cljam.fasta.reader.FASTAReader. The second argument, ref-gene, must be a path
  to refGene.txt(.gz) or ref-gene index. alt must be a single alternation such
  as \"TG\". \"TG,T\", for example, is not allowed. A returned sequence consists
  of maps, each having :cdna and :protein HGVS defined in clj-hgvs."
  {:arglists '([ref-fa ref-gene chr pos ref alt])}
  (fn [ref-fa ref-gene chr pos ref alt]
    (class ref-fa)))

(defmethod vcf-variant->hgvs String
  [ref-fa ref-gene chr pos ref alt]
  (with-open [fa-rdr (fa/reader ref-fa)]
    (doall (vcf-variant->hgvs fa-rdr ref-gene chr pos ref alt))))

(defmethod vcf-variant->hgvs cljam.fasta.reader.FASTAReader
  [ref-fa ref-gene chr pos ref alt]
  (if (valid-ref? ref-fa chr pos ref)
    (let [rgidx (if (map? ref-gene)
                  ref-gene
                  (rg/index (rg/load-ref-genes ref-gene)))]
      (->> (rg/ref-genes chr pos rgidx)
           (map (fn [rg]
                  (assoc (normalize-variant {:chr chr, :pos pos, :ref ref, :alt alt}
                                            ref-fa rg)
                         :rg rg)))
           (map (fn [{:keys [rg chr pos ref alt]}]
                  {:cdna (cdna/rg->hgvs ref-fa rg pos ref alt)
                   :protein (if (rg/in-cds? pos rg)
                              (prot/rg->hgvs ref-fa rg pos ref alt))}))
           distinct))
    (throw (Exception. (format "\"%s\" is not found on %s:%d" ref chr pos)))))
