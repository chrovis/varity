(ns varity.vcf-to-hgvs-test
  (:require [clojure.test :refer :all]
            [clj-hgvs.core :as hgvs]
            [varity.ref-gene :as rg]
            [varity.vcf-to-hgvs :refer :all]
            [varity.t-common :refer :all]))

(use-fixtures :once disable-log-fixture)

(defn- vcf-variant->cdna-hgvs-texts
  [variant seq-rdr rgidx & [options]]
  (map #(hgvs/format % {:show-bases? true
                        :range-format :coord})
       (vcf-variant->cdna-hgvs variant seq-rdr rgidx options)))

(defslowtest vcf-variant->cdna-hgvs-test
  (cavia-testing "returns cDNA HGVS strings"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [chr pos ref alt e]
          (= (vcf-variant->cdna-hgvs-texts {:chr chr, :pos pos, :ref ref, :alt alt}
                                           test-ref-seq-file rgidx) e)
        ;; Substitution
        "chr7" 55191822 "T" "G" '("NM_005228:c.2573T>G") ; cf. rs121434568 (+)
        "chr1" 11796321 "G" "A" '("NM_005957:c.665C>T") ; cf. rs1801133 (-)
        "chr12" 25245350 "C" "A" '("NM_004985:c.35G>T"
                                   "NM_033360:c.35G>T") ; cf. rs121913529 (-)
        "chr17" 7676147 "G" "A" '("NM_001126115:c.-932C>T"
                                  "NM_001126116:c.-932C>T"
                                  "NM_001126117:c.-932C>T"
                                  "NM_001276697:c.-1013C>T"
                                  "NM_001276698:c.-1013C>T"
                                  "NM_001276699:c.-1013C>T"
                                  "NM_000546:c.222C>T"
                                  "NM_001126112:c.222C>T"
                                  "NM_001126113:c.222C>T"
                                  "NM_001126114:c.222C>T"
                                  "NM_001126118:c.105C>T"
                                  "NM_001276695:c.105C>T"
                                  "NM_001276696:c.105C>T"
                                  "NM_001276760:c.105C>T"
                                  "NM_001276761:c.105C>T") ; rs786201577 (synonymous)
        "chr17" 7676154 "G" "G" '("NM_001126115:c.-939C="
                                  "NM_001126116:c.-939C="
                                  "NM_001126117:c.-939C="
                                  "NM_001276697:c.-1020C="
                                  "NM_001276698:c.-1020C="
                                  "NM_001276699:c.-1020C="
                                  "NM_000546:c.215C="
                                  "NM_001126112:c.215C="
                                  "NM_001126113:c.215C="
                                  "NM_001126114:c.215C="
                                  "NM_001126118:c.98C="
                                  "NM_001276695:c.98C="
                                  "NM_001276696:c.98C="
                                  "NM_001276760:c.98C="
                                  "NM_001276761:c.98C=") ; cf. rs1042522

        ;; Deletion
        "chr1" 963222 "GCG" "G" '("NM_015658:c.-3984_-3983delCG"
                                  "NM_198317:c.1157_1158delCG"
                                  "NM_001160184:c.-3309_-3308delCG"
                                  "NM_032129:c.-3309_-3308delCG")
        "chr3" 116193879 "ACAC" "A" '("NM_001318915:c.156-107326_156-107324delGTG"
                                      "NM_002338:c.156-107326_156-107324delGTG") ; cf. rs17358
        "chr16" 282999 "GTCTC" "G" '("NM_001176:c.*190_*193delTCTC"
                                     "NM_006849:c.-166_-163delTCTC"
                                     "NM_003502:c.*5119_*5122delGAGA"
                                     "NM_181050:c.*5119_*5122delGAGA") ; cf. rs28365940 (deletion in UTR)
        "chr6" 33086236 "TA" "T" '("NM_002121:c.776delA") ; cf. rs67523850 (deletion in border of UTR)

        ;; Duplication
        "chr2" 47806842 "T" "TGACT" '("NM_000179:c.4062_4065dupGACT"
                                      "NM_001281492:c.3672_3675dupGACT"
                                      "NM_001281493:c.3156_3159dupGACT"
                                      "NM_001281494:c.3156_3159dupGACT"
                                      "NM_025133:c.*1276_*1279dupAGTC"
                                      "NM_001190274:c.*1276_*1279dupAGTC") ; cf. rs55740729 (+)
        "chr2" 26254257 "G" "GACT" '("NM_000183:c.4_6dupACT"
                                     "NM_001281512:c.4_6dupACT"
                                     "NM_001281513:c.-146_-144dupACT") ; cf. rs3839049 (+)
        "chr1" 42752620 "T" "TGGAGTTC" '("NM_001146289:c.1383_1389dupGAACTCC"
                                         "NM_001243246:c.1383_1389dupGAACTCC"
                                         "NM_022356:c.1383_1389dupGAACTCC") ; cf. rs137853953 (-)
        "chr7" 152247986 "G" "GT" '("NM_170606:c.2447dupA") ; cf. rs150073007 (-)

        ;; Insertion
        "chr1" 69567 "A" "AT" '("NM_001005484:c.477_478insT")
        "chr3" 122740443 "G" "GAGA" '("NM_024610:c.1368_1369insTCT"
                                      "NM_001320728:c.1284_1285insTCT") ; cf. rs16338 (-)

        ;; inversion
        "chr2" 47806747 "AAAACTTTTTTTTTTTTTTTTTTAA" "ATTAAAAAAAAAAAAAAAAAAGTTT"
        '("NM_000179:c.4002-31_4002-8inv"
          "NM_001281492:c.3612-31_3612-8inv"
          "NM_001281493:c.3096-31_3096-8inv"
          "NM_001281494:c.3096-31_3096-8inv"
          "NM_025133:c.*1347_*1370inv"
          "NM_001190274:c.*1347_*1370inv") ; cf. rs267608133 (+)
        ;; NOTE: strand - example is not found on dbSNP

        ;; indel
        "chr3" 37006994 "AAG" "AGTT" '("NM_000249:c.385_386delAGinsGTT"
                                       "NM_001258271:c.385_386delAGinsGTT"
                                       "NM_001258273:c.-339_-338delAGinsGTT"
                                       "NM_001167617:c.91_92delAGinsGTT"
                                       "NM_001167618:c.-339_-338delAGinsGTT"
                                       "NM_001167619:c.-247_-246delAGinsGTT"
                                       "NM_001258274:c.-339_-338delAGinsGTT") ; cf. rs63751710 (+)
        "chr1" 21887514 "CTG" "CC" '("NM_001291860:c.862_863delCAinsG"
                                     "NM_005529:c.862_863delCAinsG") ; cf. rs2010297 (-)

        ;; repeated sequences
        "chr7" 55191822 "T" "TGCTGCT" '("NM_005228:c.2571_2573[3]") ; not actual example (+)
        "chr3" 126492636 "C" "CCTCT" '("NM_001165974:c.1690-122_1690-121[3]"
                                       "NM_144639:c.1510-122_1510-121[3]") ; cf. rs2307882 (-)
        "chr2" 237363239 "T" "TA" '("NM_004369:c.6063+6[9]"
                                    "NM_057166:c.4242+6[9]"
                                    "NM_057167:c.5445+6[9]") ; cf. rs11385011 (-)
        )))
  (cavia-testing "tx-margin"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [chr pos ref alt tx-margin e]
          (= (vcf-variant->cdna-hgvs-texts {:chr chr, :pos pos, :ref ref, :alt alt}
                                           test-ref-seq-file rgidx
                                           {:tx-margin tx-margin}) e)
        "chr5" 1295113 "G" "A" 5000 '("NM_001193376:c.-124C>T"
                                      "NM_198253:c.-124C>T")
        "chr5" 1295113 "G" "A" 0 '())))
  (cavia-testing "throws Exception if inputs are illegal"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (is (thrown? Exception
                   (vcf-variant->cdna-hgvs {:chr "chr7", :pos 55191823, :ref "T", :alt "G"}
                                           test-ref-seq-file rgidx))))))

(defn- vcf-variant->protein-hgvs-texts
  [variant seq-rdr rgidx]
  (map #(hgvs/format % {:amino-acid-format :short
                        :ter-format :short})
       (vcf-variant->protein-hgvs variant seq-rdr rgidx)))

(defslowtest vcf-variant->protein-hgvs-test
  (cavia-testing "returns protein HGVS strings"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [chr pos ref alt e]
          (= (vcf-variant->protein-hgvs-texts {:chr chr, :pos pos, :ref ref, :alt alt}
                                              test-ref-seq-file rgidx) e)
        ;; Substitution
        "chr7" 55191822 "T" "G" '("p.L858R") ; cf. rs121434568
        "chr1" 11796321 "G" "A" '("p.A222V") ; cf. rs1801133
        "chr12" 25245350 "C" "A" '("p.G12V") ; cf. rs121913529 (-)
        "chr17" 7676147 "G" "A" '("p.A74="
                                  "p.A35=") ; cf. rs786201577 (synonymous)
        "chr6" 33086236 "TA" "T" '("p.*259=") ; cf. rs67523850 (deletion in border of UTR)
        "chr7" 152247986 "G" "GT" '("p.Y816*") ; cf. rs150073007 (-, nonsense mutation)
        "chr17" 31159027 "TGC" "T" '("p.A75*") ; not actual example (+, nonsense in del case)
        "chr2" 47478341 "TG" "T" '("p.L762*" "p.L696*") ;; rs786204050 (+) frameshift with termination

        ;; deletion
        "chr1" 1286041 "CCTT" "C" '("p.F227del")
        "chr3" 198153259 "GGCAGCAGCA" "G" '("p.Q82_Q84del"); cf. rs56683636 (+)
        "chr1" 247815239 "AAGG" "A" '("p.S163del") ; cf. rs35979231 (-)
        "chr1" 84574315 "CGCAGCGCCA" "C" '("p.L31_L33del"); cf. rs3217269 (-)

        ;; Duplication
        "chr2" 26254257 "G" "GACT" '("p.T2dup") ; cf. rs3839049 (+)
        "chr2" 233521093 "T" "TTTC" '("p.K752dup") ; cf. rs59586144 (-)

        ;; Insertion
        "chr3" 73062352 "T" "TTGG" '("p.=" "p.L91_E92insV") ; cf. rs143235716 (+)
        "chr3" 122740443 "G" "GAGA" '("p.P456_Q457insS"
                                      "p.P428_Q429insS") ; cf. rs71270423 (-)

        ;; indel
        "chr2" 47445589 "CTTACTGAT" "CCC" '("p.L440_D442delinsP" "p.L374_D376delinsP") ; cf. rs63749931 (+)
        "chr1" 152111364 "TGC" "TCG" '("p.E617_Q618delinsDE") ; cf. rs35444647 (-)

        ;; repeated sequences
        "chr7" 55191823 "G" "GCTGCTG" '("p.L858[3]") ; not actual example (+)
        "chr1" 11796319 "C" "CGGCGGC" '("p.A222[3]") ; not actual example (-)

        ;; Frame shift
        "chr1" 69567 "A" "AT" '("p.L160Sfs*7")
        "chr1" 963222 "GCG" "G" '("p.A386Gfs*12")
        "chr2" 47478341 "T" "TGG" '("p.L762Gfs*2" "p.L696Gfs*2")

        ;; Extension
        "chr2" 189011772 "T" "C" '("p.*1467Qext*45") ; cf. ClinVar 101338
        ;; NOTE: There are very few correct examples...

        ;; no effect
        "chr7" 55181876 "A" "T" '("p.=") ; not actual example (+)
        "chr7" 55181874 "TGAT" "T" '("p.=") ; not actual example (+)
        "chr7" 55181876 "A" "AGGT" '("p.=") ; not actual example (+)

        ;; unknown
        "chr12" 40393453 "G" "A" '("p.?") ; not actual example (+)
        )))
  (cavia-testing "throws Exception if inputs are illegal"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (is (thrown? Exception
                   (vcf-variant->protein-hgvs {:chr "chr7", :pos 55191823, :ref "T", :alt "G"}
                                              test-ref-seq-file rgidx))))))
