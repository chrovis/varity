(ns varity.vcf-to-hgvs-bench
  (:require [cljam.io.sequence :as cseq]
            [libra.bench :refer :all]
            [libra.criterium :as c]
            [varity.ref-gene :as rg]
            [varity.vcf-to-hgvs :as v2h]
            [varity.t-common :refer :all]))

(def variants
  [{:chr "chr7", :pos 55191822, :ref "T", :alt "G"}
   {:chr "chr1", :pos 247815239, :ref "AAGG", :alt "A"}
   {:chr "chr2", :pos 26254257, :ref "G", :alt "GACT"}
   {:chr "chr3", :pos 122740443, :ref "G", :alt "GAGA"}
   {:chr "chr2", :pos 47445589, :ref "CTTACTGAT", :alt "CCC"}
   {:chr "chr1", :pos 11796319, :ref "C", :alt "CGGCGGC"}
   {:chr "chr1", :pos 69567, :ref "A", :alt "AT"}
   {:chr "chr2", :pos 189011772, :ref "T", :alt "C"}])

(defbench ^:slow vcf-variant->protein-hgvs-bench
  (prepare-cavia!)
  (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
    (with-open [seq-rdr (cseq/reader test-ref-seq-file)]
      (is (c/quick-bench (doseq [v variants]
                           (doall (v2h/vcf-variant->protein-hgvs v seq-rdr rgidx))))))))
