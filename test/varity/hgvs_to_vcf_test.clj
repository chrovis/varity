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
        ;; substitution
        "c.2573T>G" "EGFR" '({:chr "chr7", :pos 55191822, :ref "T", :alt "G"}) ; cf. rs121434568
        "c.665C>T" "MTHFR" '({:chr "chr1", :pos 11796321, :ref "G", :alt "A"}) ; cf. rs1801133

        ;; deletion
        "c.1157_1158delCG" "KLHL17" '({:chr "chr1", :pos 963222, :ref "GCG", :alt "G"})
        "c.156-107326_156-107324delGTG" "LSAMP"'({:chr "chr3", :pos 116193879, :ref "ACAC", :alt "A"}) ; cf. rs17358

        ;; duplication
        "c.4062_4065dupGACT" "MSH6" '({:chr "chr2", :pos 47806842, :ref "T", :alt "TGACT"}) ; cf. rs55740729 (+)
        "c.4_6dupACT" "HADHB" '({:chr "chr2", :pos 26254260, :ref "T", :alt "TACT"}
                                {:chr "chr2", :pos 26261026, :ref "A", :alt "AACA"}) ; cf. rs3839049 (+)
        "c.1383_1389dupGAACTCC" "P3H1" '({:chr "chr1", :pos 42752620, :ref "T", :alt "TGGAGTTC"}) ; cf. rs137853953 (-)

        ;; insertion
        "c.477_478insT" "OR4F5" '({:chr "chr1", :pos 69567, :ref "A", :alt "AT"})
        "c.1368_1369insTCT" "HSPBAP1" '({:chr "chr3", :pos 122740443, :ref "G", :alt "GAGA"}
                                        {:chr "chr3", :pos 122740359, :ref "C", :alt "CAGA"}) ; cf. rs16338 (-)

        ;; inversion
        "c.4002-31_4002-8inv" "MSH6" '({:chr "chr2", :pos 47806747, :ref "AAAACTTTTTTTTTTTTTTTTTTAA", :alt "ATTAAAAAAAAAAAAAAAAAAGTTT"}) ; cf. rs267608133 (+)

        ;; indel
        "c.385_386delAGinsGTT" "MLH1" '({:chr "chr3", :pos 37006994, :ref "AAG", :alt "AGTT"}) ; cf. rs63751710 (+)
        "c.862_863delCAinsG" "HSPG2" '({:chr "chr1", :pos 21887514, :ref "CTG", :alt "CC"}) ; cf. rs2010297 (-)

        ;; repeated sequences
        "c.2571_2573[3]" "EGFR" '({:chr "chr7", :pos 55191822, :ref "T", :alt "TGCTGCT"})
        "c.1510-122_1510-121[3]" "UROC1" '({:chr "chr3", :pos 126492636, :ref "C", :alt "CCTCT"}) ; cf. rs2307882 (-)
        "c.6063+6[9]" "COL6A3" '({:chr "chr2", :pos 237363246, :ref "A", :alt "AA"}
                                 {:chr "chr2", :pos 237341025, :ref "T", :alt "TGGGGG"}
                                 {:chr "chr2", :pos 237353343, :ref "G", :alt "GTTTTTT"}) ; cf. rs11385011 (-)
        )))
  (cavia-testing "protein HGVS to possible vcf variants"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [hgvs* gene e]
          (= (hgvs->vcf-variants (hgvs/parse hgvs*) gene test-fa-file rgidx) e)
        "p.L858R" "EGFR" '({:chr "chr7", :pos 55191822, :ref "T", :alt "G"}) ; cf. rs121434568
        "p.A222V" "MTHFR" '({:chr "chr1", :pos 11796321, :ref "G", :alt "A"}) ; cf. rs1801133
        ))))
