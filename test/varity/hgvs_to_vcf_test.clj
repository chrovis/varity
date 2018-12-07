(ns varity.hgvs-to-vcf-test
  (:require [clojure.test :refer :all]
            [clj-hgvs.core :as hgvs]
            [varity.ref-gene :as rg]
            [varity.hgvs-to-vcf :refer :all]
            [varity.t-common :refer :all]))

(defslowtest hgvs->vcf-variants-test
  (cavia-testing "cDNA HGVS to vcf variants"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [hgvs* e]
          (= (hgvs->vcf-variants (hgvs/parse hgvs*) test-ref-seq-file rgidx) e)
        ;; substitution
        "NM_005228:c.2573T>G" '({:chr "chr7", :pos 55191822, :ref "T", :alt "G"}) ; cf. rs121434568
        "NM_005957:c.665C>T" '({:chr "chr1", :pos 11796321, :ref "G", :alt "A"}) ; cf. rs1801133
        "NM_000546:c.215C=" '({:chr "chr17", :pos 7676154, :ref "G", :alt "G"}) ; cf. rs1042522
        "NM_198253:c.-124C>T" '({:chr "chr5", :pos 1295113, :ref "G", :alt "A"}) ; promoter

        ;; deletion
        "NM_198317:c.1157_1158delCG" '({:chr "chr1", :pos 963222, :ref "GCG", :alt "G"})
        "NM_002338:c.156-107326_156-107324delGTG" '({:chr "chr3", :pos 116193879, :ref "ACAC", :alt "A"}) ; cf. rs17358

        ;; duplication
        "NM_000179:c.4062_4065dupGACT" '({:chr "chr2", :pos 47806842, :ref "T", :alt "TGACT"}) ; cf. rs55740729 (+)
        "NM_000183:c.4_6dupACT" '({:chr "chr2", :pos 26254260, :ref "T", :alt "TACT"}) ; cf. rs3839049 (+)
        "NM_022356:c.1383_1389dupGAACTCC" '({:chr "chr1", :pos 42752620, :ref "T", :alt "TGGAGTTC"}) ; cf. rs137853953 (-)

        ;; insertion
        "NM_001005484:c.477_478insT" '({:chr "chr1", :pos 69567, :ref "A", :alt "AT"})
        "NM_024610:c.1368_1369insTCT" '({:chr "chr3", :pos 122740443, :ref "G", :alt "GAGA"}) ; cf. rs16338 (-)

        ;; inversion
        "NM_000179:c.4002-31_4002-8inv" '({:chr "chr2", :pos 47806747, :ref "AAAACTTTTTTTTTTTTTTTTTTAA", :alt "ATTAAAAAAAAAAAAAAAAAAGTTT"}) ; cf. rs267608133 (+)

        ;; indel
        "NM_000249:c.385_386delAGinsGTT" '({:chr "chr3", :pos 37006994, :ref "AAG", :alt "AGTT"}) ; cf. rs63751710 (+)
        "NM_005529:c.862_863delCAinsG" '({:chr "chr1", :pos 21887514, :ref "CTG", :alt "CC"}) ; cf. rs2010297 (-)

        ;; repeated sequences
        "NM_005228:c.2571_2573[3]" '({:chr "chr7", :pos 55191822, :ref "T", :alt "TGCTGCT"})
        "NM_144639:c.1510-122_1510-121[3]" '({:chr "chr3", :pos 126492636, :ref "C", :alt "CCTCT"}) ; cf. rs2307882 (-)
        "NM_004369:c.6063+6[9]" '({:chr "chr2", :pos 237363246, :ref "A", :alt "AA"}) ; cf. rs11385011 (-)
        )))
  (cavia-testing "cDNA HGVS with gene to vcf variants"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [hgvs* gene e]
          (= (hgvs->vcf-variants (hgvs/parse hgvs*) gene test-ref-seq-file rgidx) e)
        ;; substitution
        "c.2573T>G" "EGFR" '({:chr "chr7", :pos 55191822, :ref "T", :alt "G"}) ; cf. rs121434568
        "c.665C>T" "MTHFR" '({:chr "chr1", :pos 11796321, :ref "G", :alt "A"}) ; cf. rs1801133
        "c.215C=" "TP53" '({:chr "chr17", :pos 7674920, :ref "G", :alt "G"}
                           {:chr "chr17", :pos 7674271, :ref "G", :alt "G"}
                           {:chr "chr17", :pos 7676154, :ref "G", :alt "G"}
                           {:chr "chr17", :pos 7676037, :ref "G", :alt "G"}) ; cf. rs1042522
        "c.-124C>T" "TERT" '({:chr "chr5", :pos 1295113, :ref "G", :alt "A"}) ; promoter

        ;; deletion
        "c.1157_1158delCG" "KLHL17" '({:chr "chr1", :pos 963222, :ref "GCG", :alt "G"})
        "c.156-107326_156-107324delGTG" "LSAMP" '({:chr "chr3", :pos 116193879, :ref "ACAC", :alt "A"}) ; cf. rs17358

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
  (cavia-testing "protein HGVS with gene to possible vcf variants"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [hgvs* gene e]
          (= (hgvs->vcf-variants (hgvs/parse hgvs*) gene test-ref-seq-file rgidx) e)
        "p.L858R" "EGFR" '({:chr "chr7", :pos 55191822, :ref "T", :alt "G"}) ; cf. rs121434568
        "p.A222V" "MTHFR" '({:chr "chr1", :pos 11796321, :ref "G", :alt "A"}) ; cf. rs1801133
        "p.Q61K" "NRAS" '({:chr "chr1", :pos 114713909, :ref "G", :alt "T"}) ; cf. rs121913254
        "p.Q61K" "KRAS" '({:chr "chr12", :pos 25227343, :ref "G", :alt "T"}) ; cf. rs121913238
        )))
  (cavia-testing "protein HGVS with gene to possible vcf variants with cDNA HGVS"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [hgvs* gene e]
          (= (protein-hgvs->vcf-variants-with-cdna-hgvs (hgvs/parse hgvs*) gene test-ref-seq-file rgidx) e)
        "p.T790M" "EGFR" `({:vcf {:chr "chr7", :pos 55181378, :ref "C", :alt "T"} ; cf. rs121434569
                            :cdna ~(hgvs/parse "NM_005228:c.2369C>T")})
        "p.L1196M" "ALK" `({:vcf {:chr "chr2", :pos 29220765, :ref "G", :alt "T"} ; cf. rs1057519784
                            :cdna ~(hgvs/parse "NM_004304:c.3586C>A")})))))
