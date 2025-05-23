(ns varity.vcf-to-hgvs-test
  (:require [clojure.test :refer :all]
            [clj-hgvs.core :as hgvs]
            [varity.ref-gene :as rg]
            [varity.vcf-to-hgvs :refer :all :as v2h]
            [varity.t-common :refer [test-ref-seq-file
                                     test-ref-gene-file
                                     cavia-testing
                                     disable-log-fixture
                                     defslowtest]]))

(deftest var->start-end-cds-coord-test
  (let [forward-rg {:strand :forward
                    :tx-start 10
                    :tx-end 40
                    :cds-start 15
                    :cds-end 35
                    :exon-ranges [[10 20] [24 40]]}
        reverse-rg {:strand :reverse
                    :tx-start 10
                    :tx-end 40
                    :cds-start 15
                    :cds-end 35
                    :exon-ranges [[10 20] [24 40]]}
        start-end-coord->pos-region (fn [{:keys [start-cds-coord end-cds-coord]}]
                                      (let [->pos-region #(select-keys % [:position :region])]
                                        {:start-cds-coord (->pos-region start-cds-coord)
                                         :end-cds-coord (->pos-region end-cds-coord)}))]
    (testing "forward"
      (are [v p] (= p (start-end-coord->pos-region (#'v2h/var->start-end-cds-coord v forward-rg)))
        {:pos 15 :ref "A" :alt "G"} {:start-cds-coord {:position 1 :region nil}
                                     :end-cds-coord {:position 1 :region nil}}
        {:pos 15 :ref "A" :alt "AG"} {:start-cds-coord {:position 1 :region nil}
                                      :end-cds-coord {:position 1 :region nil}}
        {:pos 14 :ref "GATG" :alt "G"} {:start-cds-coord {:position 1 :region nil}
                                        :end-cds-coord {:position 3 :region nil}}
        {:pos 15 :ref "A" :alt "TG"} {:start-cds-coord {:position 1 :region nil}
                                      :end-cds-coord {:position 1 :region nil}}
        {:pos 14 :ref "G" :alt "GT"} {:start-cds-coord {:position 1 :region :upstream}
                                      :end-cds-coord {:position 1 :region :upstream}}
        {:pos 34 :ref "A" :alt "AT"} {:start-cds-coord {:position 17 :region nil}
                                      :end-cds-coord {:position 17 :region nil}}
        {:pos 35 :ref "A" :alt "AT"} {:start-cds-coord {:position 1 :region :downstream}
                                      :end-cds-coord {:position 1 :region :downstream}}))
    (testing "reverse"
      (are [v p] (= p (start-end-coord->pos-region (#'v2h/var->start-end-cds-coord v reverse-rg)))
        {:pos 35 :ref "T" :alt "G"} {:start-cds-coord {:position 1 :region nil}
                                     :end-cds-coord {:position 1 :region nil}}
        {:pos 34 :ref "A" :alt "ACC"} {:start-cds-coord {:position 1 :region nil}
                                       :end-cds-coord {:position 1 :region nil}}
        {:pos 33 :ref "CAT" :alt "C"} {:start-cds-coord {:position 2 :region nil}
                                       :end-cds-coord {:position 1 :region nil}}
        {:pos 35 :ref "T" :alt "AG"} {:start-cds-coord {:position 1 :region nil}
                                      :end-cds-coord {:position 1 :region nil}}
        {:pos 35 :ref "T" :alt "TC"} {:start-cds-coord {:position 1 :region :upstream}
                                      :end-cds-coord {:position 1 :region :upstream}}
        {:pos 15 :ref "T" :alt "TT"} {:start-cds-coord {:position 17 :region nil}
                                      :end-cds-coord {:position 17 :region nil}}
        {:pos 14 :ref "T" :alt "TT"} {:start-cds-coord {:position 1 :region :downstream}
                                      :end-cds-coord {:position 1 :region :downstream}}))))

(use-fixtures :once disable-log-fixture)

(defn- vcf-variant->coding-dna-hgvs-texts
  [variant seq-rdr rgidx & [options]]
  (map #(hgvs/format % {:show-bases? true
                        :range-format :coord})
       (vcf-variant->coding-dna-hgvs variant seq-rdr rgidx options)))

(defslowtest vcf-variant->coding-dna-hgvs-test
  (cavia-testing "returns coding DNA HGVS strings"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [chr pos ref alt e]
           (= (vcf-variant->coding-dna-hgvs-texts {:chr chr, :pos pos, :ref ref, :alt alt}
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
        "chr1" 963222 "GCG" "G" '("NM_015658:c.-3983_-3982delGC"
                                  "NM_198317:c.1157_1158delCG"
                                  "NM_001160184:c.-3309_-3308delCG"
                                  "NM_032129:c.-3309_-3308delCG")
        "chr7" 140800463 "CT" "C" '("NM_004333:c.878delA")
        "chr16" 280533 "CTCTCTGCCGG" "C" '("NM_001286485:c.-4632_-4623delCCGGCAGAGA"
                                           "NM_001286486:c.-5019_-5010delCCGGCAGAGA"
                                           "NM_003834:c.-4706_-4697delCCGGCAGAGA"
                                           "NM_183337:c.-4632_-4623delCCGGCAGAGA"
                                           "NM_001176:c.-147_-138delTCTCTGCCGG"
                                           "NM_006849:c.-2636_-2627delTCTCTGCCGG") ; cf. rs1460727826 (deletion in UTR)
        "chr6" 33086236 "TA" "T" '("NM_002121:c.776delA") ; cf. rs67523850 (deletion in border of UTR)
        "chr17" 7675170 "CAA" "C" '("NM_001126115:c.44_45delTT"
                                    "NM_001126116:c.44_45delTT"
                                    "NM_001126117:c.44_45delTT"
                                    "NM_001276697:c.-38_-37delTT"
                                    "NM_001276698:c.-38_-37delTT"
                                    "NM_001276699:c.-38_-37delTT"
                                    "NM_000546:c.440_441delTT"
                                    "NM_001126112:c.440_441delTT"
                                    "NM_001126113:c.440_441delTT"
                                    "NM_001126114:c.440_441delTT"
                                    "NM_001126118:c.323_324delTT"
                                    "NM_001276695:c.323_324delTT"
                                    "NM_001276696:c.323_324delTT"
                                    "NM_001276760:c.323_324delTT"
                                    "NM_001276761:c.323_324delTT")

        ;; Duplication
        "chr2" 47806844 "T" "TGATT" '("NM_000179:c.4068_4071dupGATT"
                                      "NM_001281492:c.3678_3681dupGATT"
                                      "NM_001281493:c.3162_3165dupGATT"
                                      "NM_001281494:c.3162_3165dupGATT"
                                      "NM_025133:c.*1272_*1275dupTCAA"
                                      "NM_001190274:c.*1272_*1275dupTCAA") ; cf. rs55740729 (+)
        "chr2" 26254257 "G" "GACT" '("NM_000183:c.5_7dupCTA"
                                     "NM_001281512:c.5_7dupCTA"
                                     "NM_001281513:c.-145_-143dupCTA") ; cf. rs3839049 (+)
        "chr1" 42752620 "T" "TGGAGTTC" '("NM_001146289:c.1383_1389dupGAACTCC"
                                         "NM_001243246:c.1383_1389dupGAACTCC"
                                         "NM_022356:c.1383_1389dupGAACTCC") ; cf. rs137853953 (-)
        "chr7" 152247986 "G" "GT" '("NM_170606:c.2447dupA") ; cf. rs150073007 (-)
        "chr5" 112839958 "A" "AA" '("NM_001127511:c.4310dupA"
                                    "NM_000038:c.4364dupA"
                                    "NM_001127510:c.4364dupA") ; cf. COSV57323270 (neither repeated sequences nor insertion)

        ;; Insertion
        "chr1" 69567 "A" "AT" '("NM_001005484:c.477_478insT")
        "chr3" 122740443 "G" "GAGA" '("NM_024610:c.1368_1369insTCT"
                                      "NM_001320728:c.1284_1285insTCT") ; cf. rs16338 (-)
        "chr3" 122740443 "G" "GNNN" '("NM_024610:c.1368_1369insNNN"
                                      "NM_001320728:c.1284_1285insNNN")

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
        "chr1" 21887514 "CTG" "CN" '("NM_001291860:c.862_863delCAinsN"
                                     "NM_005529:c.862_863delCAinsN")

        ;; repeated sequences
        "chr7" 55191822 "T" "TGCTGCT" '("NM_005228:c.2572_2574[3]") ; not actual example (+)
        "chr3" 126492636 "C" "CCTCT" '("NM_001165974:c.1690-122_1690-121[3]"
                                       "NM_144639:c.1510-122_1510-121[3]") ; cf. rs2307882 (-)
        "chr2" 237363239 "T" "TA" '("NM_004369:c.6063+6[9]"
                                    "NM_057166:c.4242+6[9]"
                                    "NM_057167:c.5445+6[9]") ; cf. rs11385011 (-)
        )))
  (cavia-testing "options"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [chr pos ref alt opts e] (= (vcf-variant->coding-dna-hgvs-texts
                                        {:chr chr, :pos pos, :ref ref, :alt alt}
                                        test-ref-seq-file rgidx opts)
                                       e)
        ;; prefer-deletion?, cf. rs727502907 (-)
        "chr7" 140924774 "GGGAGGC" "G" {:prefer-deletion? false} '("NM_004333:c.-95_-90[3]")
        "chr7" 140924774 "GGGAGGC" "G" {:prefer-deletion? true} '("NM_004333:c.-77_-72delGCCTCC")
        "chr20" 58854572 "CCGCCCCAGCCGATCCCGACTCCGGGGCGGCCCCTGA" "C" {:prefer-deletion? true}
        '("NM_016592:c.*42+13738_*42+13773delCGACTCCGGGGCGGCCCCTGACGCCCCAGCCGATCC"
          "NM_001309861:c.-39+12749_-39+12784delCGACTCCGGGGCGGCCCCTGACGCCCCAGCCGATCC"
          "NM_001309883:c.1172_1207delCGACTCCGGGGCGGCCCCTGACGCCCCAGCCGATCC"
          "NM_001077490:c.1172_1207delCGACTCCGGGGCGGCCCCTGACGCCCCAGCCGATCC"
          "NM_080425:c.1359_1394delCGACTCCGGGGCGGCCCCTGACGCCCCAGCCGATCC")
        "chr7" 141764158 "CATG" "C" {:prefer-deletion? false} '("NM_016943:c.1_3[1]") ; not actual example (+)
        "chr7" 141764158 "CATG" "C" {:prefer-deletion? true} '("NM_016943:c.4_6delATG") ; not actual example (+)

        ;; prefer-insertion?, cf. rs2307882 (-)
        "chr3" 126492636 "C" "CCTCT" {:prefer-insertion? false} '("NM_001165974:c.1690-122_1690-121[3]"
                                                                  "NM_144639:c.1510-122_1510-121[3]")
        "chr3" 126492636 "C" "CCTCT" {:prefer-insertion? true} '("NM_001165974:c.1690-121_1690-120insAGAG"
                                                                 "NM_144639:c.1510-121_1510-120insAGAG")

        ;; prefer-deletion-inserion?
        ;; inversion cf. rs267608133 (+)
        "chr2" 47806747 "AAAACTTTTTTTTTTTTTTTTTTAA" "ATTAAAAAAAAAAAAAAAAAAGTTT" {:prefer-deletion-insertion? true} '("NM_000179:c.4002-31_4002-8delAAACTTTTTTTTTTTTTTTTTTAAinsTTAAAAAAAAAAAAAAAAAAGTTT"
                                                                                                                     "NM_001281492:c.3612-31_3612-8delAAACTTTTTTTTTTTTTTTTTTAAinsTTAAAAAAAAAAAAAAAAAAGTTT"
                                                                                                                     "NM_001281493:c.3096-31_3096-8delAAACTTTTTTTTTTTTTTTTTTAAinsTTAAAAAAAAAAAAAAAAAAGTTT"
                                                                                                                     "NM_001281494:c.3096-31_3096-8delAAACTTTTTTTTTTTTTTTTTTAAinsTTAAAAAAAAAAAAAAAAAAGTTT"
                                                                                                                     "NM_025133:c.*1347_*1370delTTAAAAAAAAAAAAAAAAAAGTTTinsAAACTTTTTTTTTTTTTTTTTTAA"
                                                                                                                     "NM_001190274:c.*1347_*1370delTTAAAAAAAAAAAAAAAAAAGTTTinsAAACTTTTTTTTTTTTTTTTTTAA")

        ;; tx-margin
        "chr5" 1295113 "G" "A" {:tx-margin 5000} '("NM_001193376:c.-124C>T"
                                                   "NM_198253:c.-124C>T")
        "chr5" 1295113 "G" "A" {:tx-margin 0} '()

        ;; three-prime-rule
        "chr13" 24421115 "TGACTTAGCCT" "T" {:three-prime-rule {:restrict-cds true}} '("NM_006437:c.5169_*3delAGGCTAAGTC") ;; not actual example (-)
        "chr13" 24421115 "TGACTTAGCCT" "T" {:three-prime-rule {:restrict-cds false}} '("NM_006437:c.5170_*4delGGCTAAGTCA") ;; not actual example (-)
        "chr3" 53495165 "GGAT" "G" {:three-prime-rule {:restrict-cds true} :prefer-insertion? true :prefer-deletion? true} '("NM_000720:c.-1_2delGAT"
                                                                                                                             "NM_001128839:c.-1_2delGAT"
                                                                                                                             "NM_001128840:c.-1_2delGAT") ;; not actual example (+)
        "chr3" 53495165 "GGAT" "G" {:three-prime-rule {:restrict-cds false} :prefer-insertion? true :prefer-deletion? true} '("NM_000720:c.20_22delTGA"
                                                                                                                              "NM_001128839:c.20_22delTGA"
                                                                                                                              "NM_001128840:c.20_22delTGA")))) ;; not actual example (+)
  (cavia-testing "throws Exception if inputs are illegal"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (is (thrown? Exception
                   (vcf-variant->coding-dna-hgvs {:chr "chr7", :pos 55191823, :ref "T", :alt "G"}
                                                 test-ref-seq-file rgidx))))))

(defn- vcf-variant->protein-hgvs-texts
  [variant seq-rdr rgidx & [options]]
  (map #(hgvs/format % {:amino-acid-format :short
                        :show-ter-site? true
                        :ter-format :short})
       (vcf-variant->protein-hgvs variant seq-rdr rgidx options)))

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
        "chr17" 7676147 "G" "A" '("p.A74=" "p.A35=") ; cf. rs786201577 (synonymous)
        "chr11" 14279340 "G" "A" '("p.F127=" "p.F204=" "p.F169=") ; not actual example (-)
        "chr3" 149520808 "C" "CTTAA" '("p.*401=") ; not actual example (-)
        "chr7" 152247986 "G" "GT" '("p.Y816*") ; cf. rs150073007 (-, nonsense mutation)
        "chr18" 51048782 "C" "CAGT" '("p.Y117*") ; cf. not actual example (+, inframe nonsense mutation)
        "chr17" 31159027 "TGC" "T" '("p.A75*") ; not actual example (+, nonsense in del case)
        "chr2" 47478341 "TG" "T" '("p.L762*" "p.L696*") ;; rs786204050 (+) frameshift with termination
        "chr17" 7676202 "T" "TGTCCCTTAGTCTT" '("p.P58*" "p.P19*") ; cf. not actual example (-, frameshift with termination)
        "chr13" 24421115 "TGACTTAGCC" "T" '("p.G1724K") ; not actual example (-, del includes termination codon)
        "chrX" 15823233 "ATCCAAATAAAC" "A" '("p.S481*") ; not actual example (+, del includes termination codon)

        ;; Deletion
        "chr1" 240092288 "AGTC" "A" '("p.S61del") ; cf. rs772088733 (+)
        "chr7" 55174771 "AGGAATTAAGAGAAGC" "A" '("p.E746_A750del") ; cf. rs121913421 (+)
        "chr1" 247815239 "AAGG" "A" '("p.S163del") ; cf. rs35979231 (-)
        "chr2" 29223408 "AAGCAGT" "A" '("p.Y1096_C1097del") ; cf. rs776101205 (-)
        "chr11" 108259071 "GT" "G" '("p.L822*") ; https://github.com/chrovis/varity/issues/49
        "chr17" 7676206 "TGAACCATTGTTCAATATCGTCCGGGGACAGCATCAAATCATCCATTGCTTGGGACGGCAAGGGG" "T" '("p.P34Lfs*68" "p.?") ; https://mutalyzer.nl/normalizer/NC_000017.11:g.7676209_7676272del

        ;; deletion includes initiation codon
        "chr13" 24503770 "TCACCATC" "T" '("p.V2_M3del") ; not actual example(-)
        "chr13" 24503770 "TCACCA" "T" '("p.?") ; not actual example(-)

        ;; Duplication
        "chr2" 26254257 "G" "GACT" '("p.T2dup") ; cf. rs3839049 (+)
        "chr2" 233521093 "T" "TTTC" '("p.K752dup") ; cf. rs59586144 (-)

        ;; Insertion
        "chr3" 73062352 "T" "TTGG" '("p.=" "p.L91_E92insV") ; cf. rs143235716 (+)
        "chr3" 122740443 "G" "GAGA" '("p.P456_Q457insS"
                                      "p.P428_Q429insS") ; cf. rs71270423 (-)
        "chr3" 122740443 "G" "GNNN" '("p.P456_Q457insX"
                                      "p.P428_Q429insX")

        ;; deletion includes stop codon deletion
        "chr17" 9771484 "GCAGTTACC" "G" '("p.E310_G311insAGGQMGHPLEIKVFLA" "p.E311_G312insAGGQMGHPLEIKVFLA") ; not actual example(-)

        ;; indel
        "chr2" 47445589 "CTTACTGAT" "CCC" '("p.L440_D442delinsP" "p.L374_D376delinsP") ; cf. rs63749931 (+)
        "chr1" 152111364 "TGC" "TCG" '("p.E617_Q618delinsDE") ; cf. rs35444647 (-)
        "chr2" 47445589 "CTTACTGAT" "CNN" '("p.L440_D442delinsX" "p.L374_D376delinsX")

        ;; indel includes stop codon deletion
        "chr8" 116847497 "TCCTTATATAATATGGAACCTTGGTCCAGGTGTTGCGATGATGTCACTGTA" "T" '("p.Y617_I631delinsS")
        "chr17" 43045678 "TCAGTAG" "T" '("p.H758delinsQLQPATGTEPQDPKNELTKWPFQALGAPLTLQSFYCPG"
                                         "p.H1815delinsQLQPATGTEPQDPKNELTKWPFQALGAPLTLQSFYCPG"
                                         "p.H1862delinsQLQPATGTEPQDPKNELTKWPFQALGAPLTLQSFYCPG"
                                         "p.H1883delinsQLQPATGTEPQDPKNELTKWPFQALGAPLTLQSFYCPG") ; not actual example (-)
        "chr17" 43045679 "CAGTAGT" "C" '("p.H758delinsRLQPATGTEPQDPKNELTKWPFQALGAPLTLQSFYCPG"
                                         "p.H1815delinsRLQPATGTEPQDPKNELTKWPFQALGAPLTLQSFYCPG"
                                         "p.H1862delinsRLQPATGTEPQDPKNELTKWPFQALGAPLTLQSFYCPG"
                                         "p.H1883delinsRLQPATGTEPQDPKNELTKWPFQALGAPLTLQSFYCPG") ; not actual example (-)
        "chr10" 87965467 "GTCT" "G" '("p.V403delinsGIFFYQEG" "p.V576delinsGIFFYQEG" "p.V206delinsGIFFYQEG") ; not actual example (+)
        "chr10" 87965464 "AAAGTCT" "A" '("p.K402_V403delinsRIFFYQEG" "p.K575_V576delinsRIFFYQEG" "p.K205_V206delinsRIFFYQEG") ; not actual example (+)
        "chr13" 24421120 "TAGC" "T" '("p.G1724delinsEVK") ; not actual example (-)
        "chr13" 24421120 "TAGCCTT" "T" '("p.G1724delinsVK") ; not actual example (-)

        ;; repeated sequences
        "chr1" 47438996 "T" "TCCGCAC" '("p.P286_H287[5]") ; cf. rs3046924 (+)
        "chr1" 11796319 "C" "CGGCGGC" '("p.A222[3]") ; not actual example (-)

        ;; Frame shift
        "chr1" 69567 "A" "AT" '("p.L160Sfs*7")
        "chr1" 963222 "GCG" "G" '("p.A386Gfs*12")
        "chr2" 47478341 "T" "TGG" '("p.L762Gfs*2" "p.L696Gfs*2")
        "chr5" 112839837 "AGTGGCAT" "A" '("p.S1397Ifs*2"
                                          "p.S1415Ifs*2") ; https://github.com/chrovis/varity/issues/58
        "chr17" 31261816 "CC" "C" '("p.N1542Tfs*11" "p.N1563Tfs*11") ; cf. rs1555619041 (+)
        "chr1" 16138274 "CG" "C" '("p.R327Dfs*66")
        "chr3" 149520807 "ACAGC" "A" '("p.W399Cfs*62") ; not actual example (-)

        ;; initiation site is affected in DNA level but initiation site is not changed in protein level
        "chr7" 55019279 "TGC" "T" '("p.R2Tfs*34") ; not actual example (+)
        "chr11" 32396363 "C" "CGACCGTACAA" '("p.A170Cfs*6" "p.A153Cfs*6" "p.A365Cfs*6" "p.A382Cfs*6") ; not actual example (-)

        ;; Frame shift without termination site
        "chr17" 81537070 "G" "GTA" '("p.W514Cfs*?" "p.W490Cfs*?") ; not actual example (+)
        "chr17" 9771493 "CCT" "C" '("p.E310Gfs*?" "p.E311Gfs*?") ; not actual example (-)

        ;; Frame shift includes termination site deletion
        "chr17" 43045679 "CAGTAG" "C" '("p.H758Qfs*16" "p.H1815Qfs*16" "p.H1862Qfs*16" "p.H1883Qfs*16") ;; not actual example (-)
        "chr13" 24421120 "TAGCC" "T" '("p.G1724Kfs*8") ; not actual example (-)
        "chr13" 24421120 "TAGCCT" "T" '("p.G1724Sfs*36") ; not actual example (-)
        "chr13" 24421120 "TAGCCTTG" "T" '("p.Q1723Kfs*8") ; not actual example (-)
        "chr10" 87965466 "AGTCT" "A" '("p.V403Efs*12" "p.V576Efs*12" "p.V206Efs*12") ; not actual example (+)
        "chr10" 87965465 "AAGTCT" "A" '("p.V403Nfs*17" "p.V576Nfs*17" "p.V206Nfs*17") ; not actual example (+)
        "chr10" 87965463 "AAAAGTCT" "A" '("p.K402Efs*12" "p.K575Efs*12" "p.K205Efs*12") ; not actual example (+)
        "chr13" 24421117 "ACTTAGC" "A" '("p.G1724Vfs*3") ; not actual example (-)

        ;; Extension
        "chr2" 189011772 "T" "C" '("p.*1467Qext*45") ; cf. ClinVar 101338
        "chr11" 125655318 "TGA" "TAT" '("p.*477Yext*17" "p.*443Yext*17" "p.*477Yext*24")
        "chr10" 8074014 "C" "CATGGGTT" '("p.*445Yext*64" "p.*444Yext*64") ; not actual example (+)
        "chr10" 87965468 "TC" "T" '("p.*404Eext*11" "p.*577Eext*11" "p.*207Eext*11") ; not actual example (+)
        "chrX" 15823239 "ATAA" "A" '("p.*483Text*?") ; not actual example (+)
        "chr13" 24421118 "CTTA" "C" '("p.*1725Vext*2") ; not actual example (-)
        "chrX" 15823235 "CCAAA" "C" '("p.*483Lext*?") ; not actual example (+)
        ;; NOTE: There are very few correct examples...

        ;; Extension without termination site
        "chr17" 81537077 "CT" "C" '("p.*517Eext*?" "p.*493Eext*?") ; not actual example (+)
        "chr17" 9771487 "GT" "G" '("p.*312Yext*?" "p.*313Yext*?") ; not actual example (-)

        ;; Extension includes termination site deletion
        "chr13" 24421117 "ACTT" "A" '("p.*1725Fext*2") ; not actual example (-)
        "chr13" 24421116 "GACTT" "G" '("p.*1725Sext*6") ; not actual example (-)
        "chr13" 24421117 "ACT" "A" '("p.*1725Yext*35") ; not actual example (-)
        "chr13" 24421116 "GACT" "G" '("p.*1725Yext*2") ; not actual example (-)

        ;; no effect
        "chr7" 55181876 "A" "T" '("p.=") ; not actual example (+)
        "chr7" 55181874 "TGAT" "T" '("p.=") ; not actual example (+)
        "chr7" 55181876 "A" "AGGT" '("p.=") ; not actual example (+)
        "chr17" 80090386 "CAGCACGTGCATGAACAACACAGGACACACACAGCACGTGCATGAACAACACAGGACACACACA" "C" '("p.=") ; not actural example (+)
        "chr7" 55019277 "G" "GTC" '("p.=") ; not actual example (+)
        "chr17" 21042835 "T" "TG" '("p.=") ; not actual example (-)

        ;; unknown
        "chr12" 40393453 "G" "A" '("p.?") ; not actual example (+)

        ;; unknown because variant affects utr/initiation site boundary
        ;; e.g. NM_007298:c.-19_80del from BRCA Share
        "chr17" 43124016 "CCAGATGGGACACTCTAAGATTTTCTGCATAGCATTAATGACATTTTGTACTTCTTCAACGCGAAGAGCAGATAAATCCATTTCTTTCTGTTCCAATGAA" "C" '("p.?")
        "chr2" 197434979 "AGTCTTGGCGATCTTCGCCATTTT" "A" '("p.?")
        "chr3" 53495165 "GGAT" "G" '("p.?")

        ;; unknown because variant affects initiation site
        "chr2" 188974490 "A" "C" '("p.?")
        "chr11" 118436512 "CATG" "C" '("p.?")
        "chrX" 48791109 "CATG" "C" '("p.?")

        ;; unknwon because first amino acid of alt-prot-seq is M but variant affects utr/initiation site boundary in DNA level
        "chr7" 55019277 "GATGCGA" "ATG" '("p.?") ; not actual example (+)

        ;; unknown because variant includes termination site and alternative termination site is not found
        "chr17" 81537074 "GTACTGAGGC" "G" '("p.?") ; not actual example(+)
        )))

  (cavia-testing "options"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [chr pos ref alt opts e] (= (vcf-variant->protein-hgvs-texts
                                        {:chr chr, :pos pos, :ref ref, :alt alt}
                                        test-ref-seq-file rgidx opts)
                                       e)
        ;; prefer-deletion?, not actual example (+)
        "chr1" 47439008 "CCCGCAC" "C" {:prefer-deletion? false} '("p.P286_H287[3]")
        "chr1" 47439008 "CCCGCAC" "C" {:prefer-deletion? true} '("p.P292_H293del")
        "chr20" 58854572 "CCGCCCCAGCCGATCCCGACTCCGGGGCGGCCCCTGA" "C" {:prefer-deletion? true}
        '("p.P391_I402del" "p.S455_D466del")
        "chr7" 141764158 "CATG" "C" {:prefer-deletion? false} '("p.M1[1]") ; not actual example (+)
        "chr7" 141764158 "CATG" "C" {:prefer-deletion? true} '("p.M2del") ; not actual example (+)

        ;; prefer-insertion?, cf. rs3046924 (+)
        "chr1" 47438996 "T" "TCCGCAC" {:prefer-insertion? false} '("p.P286_H287[5]")
        "chr1" 47438996 "T" "TCCGCAC" {:prefer-insertion? true} '("p.H293_A294insPH")

        ;; prefer-insertion?
        "chr1" 26773690 "C" "CGCAGCA" {:prefer-insertion? true} '("p.Q1334_R1335insQQ") ;; not actual example (+)
        "chr10" 103089678 "T" "TTCC" {:prefer-insertion? true} '("p.*562Eext*1") ;; not actual example (-)

        ;; prefer-extension-for-initial-codon-alt?
        "chr10" 121593814 "CCATGGT" "C" {:prefer-extension-for-initial-codon-alt? true} '("p.M1Vext-17") ;; not actual example (-)
        "chr2" 197434979 "AGTCTTGGCGATCTTCGCCATTTT" "A" {:prefer-extension-for-initial-codon-alt? true} '("p.M1Sext-?") ;; not actual example (-)
        ;; prefer-extension-for-initial-codon-alt?, initiation codon is altered to termination codon
        "chr9" 27109592 "T" "TTTA" {:prefer-extension-for-initial-codon-alt? true} '("p.?") ;; not actual example (+)

        ;; three-prime-rule
        "chr13" 24421115 "TGACTTAGCCT" "T" {:three-prime-rule {:restrict-cds true}} '("p.G1724Nfs*6") ;; not actual example (-)
        "chr13" 24421115 "TGACTTAGCCT" "T" {:three-prime-rule {:restrict-cds false}} '("p.G1724delinsNETEF") ;; not actual example (-)
        "chr3" 53495165 "GGAT" "G" {:three-prime-rule {:restrict-cds true} :prefer-insertion? true :prefer-deletion? true} '("p.?") ;; not actual example (+)
        "chr3" 53495165 "GGAT" "G" {:three-prime-rule {:restrict-cds false} :prefer-insertion? true :prefer-deletion? true} '("p.M7del") ;; not actual example (+)
        "chr13" 24503771 "CACCATCC" "C" {:three-prime-rule {:restrict-cds true} :prefer-insertion? true :prefer-deletion? true} '() ; not actual example(-), affects exon/intron boundary
        "chr13" 24503771 "CACCATCC" "C" {:three-prime-rule {:restrict-cds false} :prefer-insertion? true :prefer-deletion? true} '("p.V2_M3del") ; not actual example(-)
        "chr6" 33086236 "TA" "T" {:three-prime-rule {:restrict-cds true} :prefer-insertion? true :prefer-deletion? true} '("p.*259=")
        "chr6" 33086236 "TA" "T" {:three-prime-rule {:restrict-cds false} :prefer-insertion? true :prefer-deletion? true} '("p.=") ; cf. rs67523850 (deletion in border of UTR)
        "chr13" 24421121 "A" "ATTA" {:three-prime-rule {:restrict-cds true} :prefer-insertion? true :prefer-deletion? true} '("p.*1725=")
        "chr13" 24421121 "A" "ATTA" {:three-prime-rule {:restrict-cds false} :prefer-insertion? true :prefer-deletion? true} '("p.=") ; not actual example (-)
        )))

  (cavia-testing "throws Exception if inputs are illegal"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (is (thrown? Exception
                   (vcf-variant->protein-hgvs {:chr "chr7", :pos 55191823, :ref "T", :alt "G"}
                                              test-ref-seq-file rgidx)))))

  (cavia-testing "case that inputs overlap exon/intron boundaries"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [chr pos ref alt]
           (= [] (vcf-variant->protein-hgvs {:chr chr, :pos pos, :ref ref, :alt alt}
                                            test-ref-seq-file rgidx))
        ;; Two variants at the each side of a GT dinucleotide splice donor site
        "chr1" 26773716 "CGGTGA" "CCAGGTGT"
        "chr1" 26773714 "AACGGTGAG" "AGCGGT"

        ;; p.=
        ;; Two synonymous snvs on the same strand.
        ;; The intron in between remains unchanged.
        "chr1" 26773714
        "AACGGTGAGTAAAGCCTGGTCTCGGTGCTGCTATGGATCAGGCTTCGCCACTGCCCACCCTAATCCTGTGTTTCTTTGCCTCCTATAGACAT"
        "AGCGGTGAGTAAAGCCTGGTCTCGGTGCTGCTATGGATCAGGCTTCGCCACTGCCCACCCTAATCCTGTGTTTCTTTGCCTCCTATAGACAC"
        ;; p.R1335delinsQRH
        ;; Two insertions (GCA at first, CAT at last) on the same strand.
        ;; The intron in between remains unchanged.
        "chr1" 26773714
        "AACGGTGAGTAAAGCCTGGTCTCGGTGCTGCTATGGATCAGGCTTCGCCACTGCCCACCCTAATCCTGTGTTTCTTTGCCTCCTATAGACAT"
        "AGCAACGGTGAGTAAAGCCTGGTCTCGGTGCTGCTATGGATCAGGCTTCGCCACTGCCCACCCTAATCCTGTGTTTCTTTGCCTCCTATAGACATCAT"

        ;; Variants include exon/intron boundary.
        ;; cf. ClinVar 965170
        "chr8" 42838217 "GAGATTAACAGGGGTCTGAAGAGGCGGCATTAGTAATCCAATAGCAGCATCAACCTGGGAAACAGGAGGCGGTAAAGGAGGTGGGGGAAGCTGTTCCTGTGGCTCCAGAAGATCTTCTTTCTAAAACAAAAATACAAAGTATGTTTGAATTTAGTAACTAAAAACAGTTTAAA" "G"
        ;; e.g. NM_001174080:c.3696-10_3710del
        "chr19" 13283381 "GCACAGGCGGCGAAGGCTGTTGGAGA" "G"))))

(deftest coding-dna-ref-gene?-test
  (testing "valid reference genes"
    (are [transcript] (true? (#'v2h/coding-dna-ref-gene? {:name transcript}))
      "NM_001005484.2"
      "NM_001005484"
      "ENST00000644969.2"
      "ENST00000644969"))
  (testing "invalid reference genes"
    (are [transcript] (false? (#'v2h/coding-dna-ref-gene? {:name transcript}))
      "XM_024451963.1"
      "XM_024451963"
      "NR_024540.1"
      "NR_024540"
      "ENSP00000496776.1"
      "ENSP00000496776")))

(defn- vcf-variant->hgvs-texts
  [variant seq-rdr rgidx & [options]]
  (map (fn [{:keys [coding-dna protein]}]
         {:coding-dna (hgvs/format coding-dna {:show-bases? true
                                               :range-format :coord})
          :protein (hgvs/format protein {:amino-acid-format :short
                                         :show-ter-site? true
                                         :ter-format :short})})
       (vcf-variant->hgvs variant seq-rdr rgidx options)))

(defslowtest vcf-variant->hgvs-test
  (cavia-testing "options"
                 (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
                   (are [chr pos ref alt opts e]
                        (= (vcf-variant->hgvs-texts {:chr chr, :pos pos, :ref ref, :alt alt}
                                                    test-ref-seq-file rgidx (merge {:prefer-insertion? true
                                                                                    :prefer-deletion? true}
                                                                                   opts)) e)
                     ;; three-prime-rule
                     "chr13" 24421115 "TGACTTAGCCT" "T" {:three-prime-rule {:coding-dna {:restrict-cds true} :protein {:restrict-cds true}}} '({:coding-dna "NM_006437:c.5169_*3delAGGCTAAGTC", :protein "p.G1724Nfs*6"}) ;; not actual example (-)
                     "chr13" 24421115 "TGACTTAGCCT" "T" {:three-prime-rule {:coding-dna {:restrict-cds false} :protein {:restrict-cds false}}} '({:coding-dna "NM_006437:c.5170_*4delGGCTAAGTCA", :protein "p.G1724delinsNETEF"}))))) ;; not actual example (-)
