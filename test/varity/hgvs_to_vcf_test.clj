(ns varity.hgvs-to-vcf-test
  (:require [clojure.test :refer :all]
            [clj-hgvs.core :as hgvs]
            [varity.ref-gene :as rg]
            [varity.hgvs-to-vcf :refer :all]
            [varity.t-common :refer :all]))

(deftest ^:slow hgvs->vcf-variants-test
  (cavia-testing "cDNA HGVS to vcf variants"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [hgvs* gene e]
          (= (hgvs->vcf-variants (hgvs/parse hgvs*) gene test-fa-file rgidx) e)
        "c.2573T>G" "EGFR" '({:chr "chr7", :pos 55191822, :ref "T", :alt "G"}) ; cf. rs121434568
        "c.665C>T" "MTHFR" '({:chr "chr1", :pos 11796321, :ref "G", :alt "A"}) ; cf. rs1801133
        )))
  (cavia-testing "protein HGVS to possible vcf variants"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [hgvs* gene e]
          (= (hgvs->vcf-variants (hgvs/parse hgvs*) gene test-fa-file rgidx) e)
        "p.L858R" "EGFR" '({:chr "chr7", :pos 55191822, :ref "T", :alt "G"}) ; cf. rs121434568
        "p.A222V" "MTHFR" '({:chr "chr1", :pos 11796321, :ref "G", :alt "A"}) ; cf. rs1801133
        ))))
