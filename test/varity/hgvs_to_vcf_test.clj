(ns varity.hgvs-to-vcf-test
  (:require [clojure.test :refer :all]
            [clj-hgvs.core :as hgvs]
            [cljam.io.sequence :as cseq]
            [varity.ref-gene :as rg]
            [varity.hgvs-to-vcf
             :refer [hgvs->vcf-variants
                     protein-hgvs->vcf-variants-with-coding-dna-hgvs]
             :as h2v]
            [varity.hgvs-to-vcf.coding-dna :as h2v-coding-dna]
            [varity.t-common :refer [cavia-testing defslowtest
                                     test-ref-gene-file
                                     test-ncbi-ref-seq-file
                                     test-ref-seq-file]]
            [varity.vcf-to-hgvs :as v2h]))

(deftest supported-transcript?-test
  (testing "supported transcript"
    (are [transcript] (true? (#'h2v/supported-transcript? transcript))
      "NM_001.3"
      "NM_112"
      "NR_001.4"
      "NR_002"

      "ENST00000497784.2"
      "ENST00000497784"
      "ENSP0001"
      "ENSP0001.11"))
  (testing "unsupported transcript"
    (are [transcript] (false? (#'h2v/supported-transcript? transcript))
      "xNM_001.3"
      "NM001"
      "NM3.3"
      "NM_.3"
      "NM_"
      "NR"
      "NMR_3.3"
      "XM_024451963.1"

      "ENST"
      "ENSP"
      "ENSTP001"

      "NM_111.1ENST001.1")))

(defn load-ref-seqs [f] (rg/load-ref-seqs f {:filter-fns [#(rg/rna-accession? (:name %))]}))

(defslowtest hgvs->vcf-variants-test
  (cavia-testing "coding DNA HGVS to vcf variants"
    (let [ref-gene-idx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [hgvs* e]
          (= (hgvs->vcf-variants (hgvs/parse hgvs*) test-ref-seq-file ref-gene-idx) e)
        ;; substitution
        "NM_005228:c.2573T>G" '({:chr "chr7", :pos 55191822, :ref "T", :alt "G"}) ; cf. rs121434568
        "NM_005957:c.665C>T" '({:chr "chr1", :pos 11796321, :ref "G", :alt "A"}) ; cf. rs1801133
        "NM_000546:c.215C=" '({:chr "chr17", :pos 7676154, :ref "G", :alt "G"}) ; cf. rs1042522
        "NM_198253:c.-124C>T" '({:chr "chr5", :pos 1295113, :ref "G", :alt "A"}) ; promoter

        ;; deletion
        "NM_198317:c.1157_1158delCG" '({:chr "chr1", :pos 963222, :ref "GCG", :alt "G"})
        "NM_002338:c.156-107326_156-107324delGTG" '({:chr "chr3", :pos 116193879, :ref "ACAC", :alt "A"}) ; cf. rs17358

        ;; duplication
        "NM_000179:c.4062_4065dupGACT" '({:chr "chr2", :pos 47806838, :ref "T", :alt "TGACT"}) ; cf. rs55740729 (+)
        "NM_000183:c.4_6dupACT" '({:chr "chr2", :pos 26254257, :ref "G", :alt "GACT"}) ; cf. rs3839049 (+)
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
        "NM_005228:c.2571_2573[3]" '({:chr "chr7", :pos 55191819, :ref "G", :alt "GGCTGCT"})
        "NM_005228:c.2571GCT[3]" '({:chr "chr7", :pos 55191819, :ref "G", :alt "GGCTGCT"})
        "NM_144639:c.1510-122_1510-121[3]" '({:chr "chr3", :pos 126492636, :ref "C", :alt "CCTCT"}) ; cf. rs2307882 (-)
        "NM_144639:c.1510-122AG[3]" '({:chr "chr3", :pos 126492636, :ref "C", :alt "CCTCT"})
        "NM_000059:c.18AG[2]" '({:chr "chr13", :pos 32316477, :ref "AAG", :alt "A"}) ; cf. rs397507623 (+)
        "NM_004333:c.-95_-90[3]" '({:chr "chr7", :pos 140924774, :ref "GGGAGGC", :alt "G"}) ; cf. rs727502907 (-)
        ))
    (let [ncbi-ref-seq-idx (rg/index (load-ref-seqs test-ncbi-ref-seq-file))]
      (are [hgvs* e]
          (= (hgvs->vcf-variants (hgvs/parse hgvs*) test-ref-seq-file ncbi-ref-seq-idx) e)
        ;; substitution
        "NM_005228.5:c.2573T>G" '({:chr "chr7", :pos 55191822, :ref "T", :alt "G"}) ; cf. rs121434568
        "NM_005957.5:c.665C>T" '({:chr "chr1", :pos 11796321, :ref "G", :alt "A"}) ; cf. rs1801133
        "NM_000546.6:c.215C=" '({:chr "chr17", :pos 7676154, :ref "G", :alt "G"}) ; cf. rs1042522
        "NM_198253.3:c.-124C>T" '({:chr "chr5", :pos 1295113, :ref "G", :alt "A"}) ; promoter

        ;; deletion
        "NM_198317.3:c.1157_1158delCG" '({:chr "chr1", :pos 963222, :ref "GCG", :alt "G"})
        "NM_002338.5:c.156-107326_156-107324delGTG" '({:chr "chr3", :pos 116193879, :ref "ACAC", :alt "A"}) ; cf. rs17358

        ;; duplication
        "NM_000179.3:c.4062_4065dupGACT" '({:chr "chr2", :pos 47806838, :ref "T", :alt "TGACT"}) ; cf. rs55740729 (+)
        "NM_000183.3:c.4_6dupACT" '({:chr "chr2", :pos 26254257, :ref "G", :alt "GACT"}) ; cf. rs3839049 (+)
        "NM_022356.4:c.1383_1389dupGAACTCC" '({:chr "chr1", :pos 42752620, :ref "T", :alt "TGGAGTTC"}) ; cf. rs137853953 (-)

        ;; insertion
        "NM_001005484.2:c.477_478insT" '({:chr "chr1", :pos 69504, :ref "G", :alt "GT"})
        "NM_024610.6:c.1368_1369insTCT" '({:chr "chr3", :pos 122740443, :ref "G", :alt "GAGA"}) ; cf. rs16338 (-)

        ;; inversion
        "NM_000179.3:c.4002-31_4002-8inv" '({:chr "chr2", :pos 47806747, :ref "AAAACTTTTTTTTTTTTTTTTTTAA", :alt "ATTAAAAAAAAAAAAAAAAAAGTTT"}) ; cf. rs267608133 (+)

        ;; indel
        "NM_000249.4:c.385_386delAGinsGTT" '({:chr "chr3", :pos 37006994, :ref "AAG", :alt "AGTT"}) ; cf. rs63751710 (+)
        "NM_005529.7:c.862_863delCAinsG" '({:chr "chr1", :pos 21887514, :ref "CTG", :alt "CC"}) ; cf. rs2010297 (-)

        ;; repeated sequences
        "NM_005228.5:c.2571_2573[3]" '({:chr "chr7", :pos 55191819, :ref "G", :alt "GGCTGCT"})
        "NM_005228.5:c.2571GCT[3]" '({:chr "chr7", :pos 55191819, :ref "G", :alt "GGCTGCT"})
        "NM_144639.3:c.1510-122_1510-121[3]" '({:chr "chr3", :pos 126492636, :ref "C", :alt "CCTCT"}) ; cf. rs2307882 (-)
        "NM_144639.3:c.1510-122AG[3]" '({:chr "chr3", :pos 126492636, :ref "C", :alt "CCTCT"})
        "NM_000059.4:c.18AG[2]" '({:chr "chr13", :pos 32316477, :ref "AAG", :alt "A"}) ; cf. rs397507623 (+)
        "NM_004333.6:c.-95_-90[3]" '({:chr "chr7", :pos 140924774, :ref "GGGAGGC", :alt "G"}) ; cf. rs727502907 (-)
        )))
  (cavia-testing "coding DNA HGVS with gene to vcf variants"
    (let [ref-gene-idx (rg/index (rg/load-ref-genes test-ref-gene-file))
          ncbi-ref-seq-idx (rg/index (load-ref-seqs test-ncbi-ref-seq-file))]
      (are [hgvs* gene e]
          (= (hgvs->vcf-variants (hgvs/parse hgvs*) gene test-ref-seq-file ref-gene-idx)
             e)
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
        "c.4062_4065dupGACT" "MSH6" '({:chr "chr2", :pos 47806838, :ref "T", :alt "TGACT"}) ; cf. rs55740729 (+)
        "c.4_6dupACT" "HADHB" '({:chr "chr2", :pos 26254257, :ref "G", :alt "GACT"}
                                {:chr "chr2", :pos 26261023, :ref "G", :alt "GACA"}) ; cf. rs3839049 (+)
        "c.1383_1389dupGAACTCC" "P3H1" '({:chr "chr1", :pos 42752620, :ref "T", :alt "TGGAGTTC"}) ; cf. rs137853953 (-)

        ;; insertion
        "c.477_478insT" "OR4F5" '({:chr "chr1", :pos 69567, :ref "A", :alt "AT"})
        "c.1368_1369insTCT" "HSPBAP1" '({:chr "chr3", :pos 122740443, :ref "G", :alt "GAGA"}
                                        {:chr "chr3", :pos 122740359, :ref "C", :alt "CAGA"}) ; cf. rs16338 (-)
        "c.1368_1369insNNN" "HSPBAP1" '({:chr "chr3", :pos 122740443, :ref "G", :alt "GNNN"}
                                        {:chr "chr3", :pos 122740359, :ref "C", :alt "CNNN"})
        "c.1368_1369insN[3]" "HSPBAP1" '({:chr "chr3", :pos 122740443, :ref "G", :alt "GNNN"}
                                         {:chr "chr3", :pos 122740359, :ref "C", :alt "CNNN"})

        ;; inversion
        "c.4002-31_4002-8inv" "MSH6" '({:chr "chr2", :pos 47806747, :ref "AAAACTTTTTTTTTTTTTTTTTTAA", :alt "ATTAAAAAAAAAAAAAAAAAAGTTT"}) ; cf. rs267608133 (+)

        ;; indel
        "c.385_386delAGinsGTT" "MLH1" '({:chr "chr3", :pos 37006994, :ref "AAG", :alt "AGTT"}) ; cf. rs63751710 (+)
        "c.862_863delCAinsG" "HSPG2" '({:chr "chr1", :pos 21887514, :ref "CTG", :alt "CC"}) ; cf. rs2010297 (-)
        "c.862_863delCAinsNNN" "HSPG2" '({:chr "chr1", :pos 21887514, :ref "CTG", :alt "CNNN"})
        "c.862_863delCAinsN[3]" "HSPG2" '({:chr "chr1", :pos 21887514, :ref "CTG", :alt "CNNN"})

        ;; repeated sequences
        "c.2571_2573[3]" "EGFR" '({:chr "chr7", :pos 55191819, :ref "G", :alt "GGCTGCT"})
        "c.2571GCT[3]" "EGFR" '({:chr "chr7", :pos 55191819, :ref "G", :alt "GGCTGCT"})
        "c.1510-122_1510-121[3]" "UROC1" '({:chr "chr3", :pos 126492636, :ref "C", :alt "CCTCT"}) ; cf. rs2307882 (-)
        "c.1510-122AG[3]" "UROC1" '({:chr "chr3", :pos 126492636, :ref "C", :alt "CCTCT"})
        "c.18AG[2]" "BRCA2" '({:chr "chr13", :pos 32316477, :ref "AAG", :alt "A"}) ; cf. rs397507623 (+)
        "c.-95_-90[3]" "BRAF" '({:chr "chr7", :pos 140924774, :ref "GGGAGGC", :alt "G"}) ; cf. rs727502907 (-)
        )
      (are [hgvs* gene e]
          (= (hgvs->vcf-variants (hgvs/parse hgvs*) gene test-ref-seq-file ncbi-ref-seq-idx)
             e)
        ;; substitution
        "c.2573T>G" "EGFR" '({:chr "chr7", :pos 55198723, :ref "T", :alt "G"}
                             {:chr "chr7", :pos 55191822, :ref "T", :alt "G"}
                             {:chr "chr7", :pos 55205358, :ref "T", :alt "G"}
                             {:chr "chr7", :pos 55198747, :ref "T", :alt "G"})
        "c.665C>T" "MTHFR" '({:alt "A", :chr "chr1", :pos 11800256, :ref "G"}
                             {:alt "A", :chr "chr1", :pos 11796321, :ref "G"})
        "c.215C=" "TP53" '({:chr "chr17", :pos 7674920, :ref "G", :alt "G"}
                           {:chr "chr17", :pos 7674271, :ref "G", :alt "G"}
                           {:chr "chr17", :pos 7676037, :ref "G", :alt "G"}
                           {:chr "chr17", :pos 7676154, :ref "G", :alt "G"})
        "c.-124C>T" "TERT" '({:chr "chr5", :pos 1295113, :ref "G", :alt "A"})

        ;; deletion
        "c.1157_1158delCG" "KLHL17" '({:chr "chr1", :pos 963222, :ref "GCG", :alt "G"})
        "c.156-107326_156-107324delGTG" "LSAMP" '({:chr "chr3", :pos 116193879, :ref "ACAC", :alt "A"})

        ;; duplication
        "c.4062_4065dupGACT" "MSH6" '({:chr "chr2", :pos 47806838, :ref "T", :alt "TGACT"})
        "c.4_6dupACT" "HADHB" '({:chr "chr2", :pos 26254257, :ref "G", :alt "GACT"}
                                {:chr "chr2", :pos 26261023, :ref "G", :alt "GACA"})
        "c.1383_1389dupGAACTCC" "P3H1" '({:chr "chr1", :pos 42752620, :ref "T", :alt "TGGAGTTC"})

        ;; insertion
        "c.477_478insT" "OR4F5" '({:alt "GT", :chr "chr1", :pos 69504, :ref "G"})
        "c.1368_1369insTCT" "HSPBAP1" '({:chr "chr3", :pos 122740443, :ref "G", :alt "GAGA"}
                                        {:chr "chr3", :pos 122740359, :ref "C", :alt "CAGA"})
        "c.1368_1369insNNN" "HSPBAP1" '({:chr "chr3", :pos 122740443, :ref "G", :alt "GNNN"}
                                        {:chr "chr3", :pos 122740359, :ref "C", :alt "CNNN"})
        "c.1368_1369insN[3]" "HSPBAP1" '({:chr "chr3", :pos 122740443, :ref "G", :alt "GNNN"}
                                         {:chr "chr3", :pos 122740359, :ref "C", :alt "CNNN"})

        ;; inversion
        "c.4002-31_4002-8inv" "MSH6" '({:chr "chr2", :pos 47806747, :ref "AAAACTTTTTTTTTTTTTTTTTTAA", :alt "ATTAAAAAAAAAAAAAAAAAAGTTT"})

        ;; indel
        "c.385_386delAGinsGTT" "MLH1" '({:alt "GGTT", :chr "chr3", :pos 37008843, :ref "GAG"}
                                        {:alt "AGTT", :chr "chr3", :pos 37006994, :ref "AAG"}
                                        {:alt "CGTT", :chr "chr3", :pos 37026005, :ref "CAG"})
        "c.862_863delCAinsG" "HSPG2" '({:chr "chr1", :pos 21887514, :ref "CTG", :alt "CC"})
        "c.862_863delCAinsNNN" "HSPG2" '({:chr "chr1", :pos 21887514, :ref "CTG", :alt "CNNN"})
        "c.862_863delCAinsN[3]" "HSPG2" '({:chr "chr1", :pos 21887514, :ref "CTG", :alt "CNNN"})

        ;; repeated sequences
        "c.2571_2573[3]" "EGFR" '({:alt "TGACGAC", :chr "chr7", :pos 55198720, :ref "T"}
                                  {:alt "GGCTGCT", :chr "chr7", :pos 55191819, :ref "G"}
                                  {:alt "ACTACTA", :chr "chr7", :pos 55205355, :ref "A"}
                                  {:alt "TTGGTGG", :chr "chr7", :pos 55198744, :ref "T"})
        "c.2571GCT[3]" "EGFR" '({:alt "TGACGAC", :chr "chr7", :pos 55198720, :ref "T"}
                                {:alt "GGCTGCT", :chr "chr7", :pos 55191819, :ref "G"}
                                {:alt "ACTACTA", :chr "chr7", :pos 55205355, :ref "A"}
                                {:alt "TTGGTGG", :chr "chr7", :pos 55198744, :ref "T"})
        "c.1510-122_1510-121[3]" "UROC1" '({:chr "chr3", :pos 126492636, :ref "C", :alt "CCTCT"})
        "c.1510-122AG[3]" "UROC1" '({:chr "chr3", :pos 126492636, :ref "C", :alt "CCTCT"})
        "c.18AG[2]" "BRCA2" '({:chr "chr13", :pos 32316477, :ref "AAG", :alt "A"})
        "c.-95_-90[3]" "BRAF" '({:alt "G", :chr "chr7", :pos 140924774, :ref "GGGAGGC"}
                                {:alt "TCGTGACCGTGAC", :chr "chr7", :pos 140924256, :ref "T"}))))
  (cavia-testing "conversion failure"
    (let [ref-gene-idx (rg/index (rg/load-ref-genes test-ref-gene-file))
          ncbi-ref-seq-idx (rg/index (load-ref-seqs test-ncbi-ref-seq-file))]
      (are [hgvs* error-type] (thrown-with-error-type?
                               error-type
                               (hgvs->vcf-variants (hgvs/parse hgvs*)
                                                   test-ref-seq-file ref-gene-idx))
        "NM_007294:c.1-?_80+?del"   ::h2v-coding-dna/ambiguous-coordinate
        "NM_000546:c.-202_-29+?dup" ::h2v-coding-dna/ambiguous-coordinate)
      (are [hgvs* gene error-type] (thrown-with-error-type?
                                    error-type
                                    (hgvs->vcf-variants (hgvs/parse hgvs*) gene
                                                        test-ref-seq-file ref-gene-idx))
        "c.1-?_80+?del"   "BRCA1" ::h2v-coding-dna/ambiguous-coordinate
        "c.-202_-29+?dup" "TP53"  ::h2v-coding-dna/ambiguous-coordinate)
      (testing "error"
        (are [error-type idx hgvs] (thrown-with-error-type?
                                    error-type
                                    (hgvs->vcf-variants (hgvs/parse hgvs)
                                                        test-ref-seq-file
                                                        idx))
          ::h2v/gene-not-found ref-gene-idx "NM_007294.4:c.1-?_80+?del"
          ::h2v/gene-not-found ncbi-ref-seq-idx  "NM_007294:c.1-?_80+?del"
          ::h2v/unsupported-hgvs-kind ref-gene-idx "NM_004006:r.6_8del"
          ::h2v/unsupported-hgvs-kind ncbi-ref-seq-idx "NM_004006.2:o.6_8del"))))
  (cavia-testing "protein HGVS with gene to possible vcf variants"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))
          ncbi-ref-seq-idx (rg/index (load-ref-seqs test-ncbi-ref-seq-file))]
      (are [hgvs* gene e]
          (= (hgvs->vcf-variants (hgvs/parse hgvs*) gene test-ref-seq-file rgidx)
             (hgvs->vcf-variants (hgvs/parse hgvs*) gene test-ref-seq-file ncbi-ref-seq-idx)
             e)
        "p.L858R" "EGFR" '({:chr "chr7", :pos 55191822, :ref "TG", :alt "GC"}
                           {:chr "chr7", :pos 55191822, :ref "TG", :alt "GT"}  ; cf. rs1057519848
                           {:chr "chr7", :pos 55191822, :ref "T", :alt "G"}    ; cf. rs121434568
                           {:chr "chr7", :pos 55191822, :ref "TG", :alt "GA"}
                           {:chr "chr7", :pos 55191821, :ref "CTG", :alt "AGA"}
                           {:chr "chr7", :pos 55191821, :ref "CT", :alt "AG"}) ; cf. rs1057519847
        "p.A222V" "MTHFR" '({:chr "chr1", :pos 11796320, :ref "GG", :alt "CA"}
                            {:chr "chr1", :pos 11796320, :ref "GG", :alt "AA"}
                            {:chr "chr1", :pos 11796320, :ref "GG", :alt "TA"}
                            {:chr "chr1", :pos 11796321, :ref "G", :alt "A"})  ; cf. rs1801133
        "p.Q61K" "NRAS" '({:chr "chr1", :pos 114713907, :ref "TTG", :alt "CTT"}
                          {:chr "chr1", :pos 114713909, :ref "G", :alt "T"})   ; cf. rs121913254
        "p.Q61K" "KRAS" '({:chr "chr12", :pos 25227341, :ref "TTG", :alt "CTT"}
                          {:chr "chr12", :pos 25227343, :ref "G", :alt "T"})   ; cf. rs121913238
        "p.K652T" "FGFR3" '({:chr "chr4", :pos 1806163, :ref "AG", :alt "CA"}
                            {:chr "chr4", :pos 1806163, :ref "AG", :alt "CC"}
                            {:chr "chr4", :pos 1806163, :ref "A", :alt "C"}    ; cf. rs121913105
                            {:chr "chr4", :pos 1806163, :ref "AG", :alt "CT"})))))

(defslowtest protein-hgvs->vcf-variants-with-coding-dna-hgvs-test
  (cavia-testing "protein HGVS with gene to possible vcf variants with coding DNA HGVS"
    (let [ref-gene-idx (rg/index (rg/load-ref-genes test-ref-gene-file))
          ncbi-ref-seq-idx (rg/index (load-ref-seqs test-ncbi-ref-seq-file))]
      (are [idx hgvs* gene e]
          (= (protein-hgvs->vcf-variants-with-coding-dna-hgvs (hgvs/parse hgvs*) gene test-ref-seq-file idx)
             e)
        ref-gene-idx "p.T790M" "EGFR" `({:vcf {:chr "chr7", :pos 55181378, :ref "C", :alt "T"} ; cf. rs121434569
                                         :coding-dna ~(hgvs/parse "NM_005228:c.2369C>T")})
        ref-gene-idx "p.L1196M" "ALK" `({:vcf {:chr "chr2", :pos 29220765, :ref "G", :alt "T"} ; cf. rs1057519784
                                         :coding-dna ~(hgvs/parse "NM_004304:c.3586C>A")})
        ref-gene-idx "p.Q61L" "NRAS" `({:vcf {:chr "chr1", :pos 114713907, :ref "TT", :alt "CA"}, ; cf. rs1057519695
                                        :coding-dna ~(hgvs/parse "NM_002524:c.182_183delAAinsTG")}
                                       {:vcf {:chr "chr1", :pos 114713907, :ref "TT", :alt "GA"},
                                        :coding-dna ~(hgvs/parse "NM_002524:c.182_183delAAinsTC")}
                                       {:vcf {:chr "chr1", :pos 114713907, :ref "TT", :alt "AA"},
                                        :coding-dna ~(hgvs/parse "NM_002524:c.182_183delAAinsTT")}
                                       {:vcf {:chr "chr1", :pos 114713908, :ref "T", :alt "A"}, ; cf. rs11554290
                                        :coding-dna ~(hgvs/parse "NM_002524:c.182A>T")}
                                       {:vcf {:chr "chr1", :pos 114713907, :ref "TTG", :alt "CAA"},
                                        :coding-dna ~(hgvs/parse "NM_002524:c.181_183delCAAinsTTG")}
                                       {:vcf {:chr "chr1", :pos 114713908, :ref "TG", :alt "AA"},
                                        :coding-dna ~(hgvs/parse "NM_002524:c.181_182delCAinsTT")})
        ref-gene-idx "p.K652T" "FGFR3" `({:vcf {:chr "chr4", :pos 1806163, :ref "AG", :alt "CA"},
                                          :coding-dna ~(hgvs/parse "NM_001163213:c.1955_1956delAGinsCA")}
                                         {:vcf {:chr "chr4", :pos 1806163, :ref "AG", :alt "CC"},
                                          :coding-dna ~(hgvs/parse "NM_001163213:c.1955_1956delAGinsCC")}
                                         {:vcf {:chr "chr4", :pos 1806163, :ref "A", :alt "C"},
                                          :coding-dna ~(hgvs/parse "NM_001163213:c.1955A>C")} ; cf. rs121913105
                                         {:vcf {:chr "chr4", :pos 1806163, :ref "AG", :alt "CT"},
                                          :coding-dna ~(hgvs/parse "NM_001163213:c.1955_1956delAGinsCT")})

        ncbi-ref-seq-idx "p.T790M" "EGFR" `({:vcf {:chr "chr7", :pos 55181378, :ref "C", :alt "T"}
                                             :coding-dna ~(hgvs/parse "NM_001346898.2:c.2369C>T")}
                                            {:vcf {:chr "chr7", :pos 55181378, :ref "C", :alt "T"} ; cf. rs121434569
                                             :coding-dna ~(hgvs/parse "NM_005228.5:c.2369C>T")})
        ncbi-ref-seq-idx "p.L1196M" "ALK" `({:vcf {:chr "chr2", :pos 29220765, :ref "G", :alt "T"} ; cf. rs1057519784
                                             :coding-dna ~(hgvs/parse "NM_004304.5:c.3586C>A")})
        ncbi-ref-seq-idx "p.Q61L" "NRAS" `({:vcf {:chr "chr1", :pos 114713907, :ref "TT", :alt "CA"}, ; cf. rs1057519695
                                            :coding-dna ~(hgvs/parse "NM_002524.5:c.182_183delAAinsTG")}
                                           {:vcf {:chr "chr1", :pos 114713907, :ref "TT", :alt "GA"},
                                            :coding-dna ~(hgvs/parse "NM_002524.5:c.182_183delAAinsTC")}
                                           {:vcf {:chr "chr1", :pos 114713907, :ref "TT", :alt "AA"},
                                            :coding-dna ~(hgvs/parse "NM_002524.5:c.182_183delAAinsTT")}
                                           {:vcf {:chr "chr1", :pos 114713908, :ref "T", :alt "A"}, ; cf. rs11554290
                                            :coding-dna ~(hgvs/parse "NM_002524.5:c.182A>T")}
                                           {:vcf {:chr "chr1", :pos 114713907, :ref "TTG", :alt "CAA"},
                                            :coding-dna ~(hgvs/parse "NM_002524.5:c.181_183delCAAinsTTG")}
                                           {:vcf {:chr "chr1", :pos 114713908, :ref "TG", :alt "AA"},
                                            :coding-dna ~(hgvs/parse "NM_002524.5:c.181_182delCAinsTT")})
        ncbi-ref-seq-idx "p.K652T" "FGFR3" `({:vcf {:chr "chr4", :pos 1806163, :ref "AG", :alt "CA"},
                                              :coding-dna ~(hgvs/parse "NM_001163213.2:c.1955_1956delAGinsCA")}
                                             {:vcf {:chr "chr4", :pos 1806163, :ref "AG", :alt "CC"},
                                              :coding-dna ~(hgvs/parse "NM_001163213.2:c.1955_1956delAGinsCC")}
                                             {:vcf {:chr "chr4", :pos 1806163, :ref "A", :alt "C"},
                                              :coding-dna ~(hgvs/parse "NM_001163213.2:c.1955A>C")} ; cf. rs121913105
                                             {:vcf {:chr "chr4", :pos 1806163, :ref "AG", :alt "CT"},
                                              :coding-dna ~(hgvs/parse "NM_001163213.2:c.1955_1956delAGinsCT")})))))

(defslowtest hgvs->vcf->hgvs-test
  (cavia-testing "protein HGVS with gene to possible vcf variants which gives the same protein HGVS"
    (let [ref-gene-idx (rg/index (rg/load-ref-genes test-ref-gene-file))
          ncbi-ref-seq-idx (rg/index (load-ref-seqs test-ncbi-ref-seq-file))]
      (with-open [r (cseq/reader test-ref-seq-file)]
        (are [?idx ?protein-hgvs ?gene-symbol]
            (->> (map :name (rg/ref-genes ?gene-symbol ?idx))
                 (mapcat
                  (fn [nm]
                    ;; hgvs -> variants with a specific accession number
                    (protein-hgvs->vcf-variants-with-coding-dna-hgvs
                     (hgvs/parse ?protein-hgvs) nm r ?idx)))
                 (map
                  (fn [{:keys [vcf] {:keys [transcript]} :coding-dna :as v}]
                    (let [hgvs (->> ?idx
                                    (rg/ref-genes transcript)
                                    first
                                    (v2h/vcf-variant->protein-hgvs vcf r))
                          fmt (hgvs/format hgvs {:amino-acid-format :short})]
                      ;; 1 variant & 1 transcript => 1 canonical hgvs
                      (assoc v :hgvs (assoc hgvs :format fmt)))))
                 (map (comp :format :hgvs))
                 (apply = ?protein-hgvs))
          ref-gene-idx "p.K652T" "FGFR3"
          ncbi-ref-seq-idx "p.K652T" "FGFR3"
          ncbi-ref-seq-idx "p.T790M" "EGFR")))))
