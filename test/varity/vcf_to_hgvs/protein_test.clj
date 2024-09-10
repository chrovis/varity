(ns varity.vcf-to-hgvs.protein-test
  (:require [clojure.string :as string]
            [clojure.test :refer :all]
            [cljam.io.sequence :as cseq]
            [clj-hgvs.core :as hgvs]
            [varity.vcf-to-hgvs.protein :as prot]
            [varity.t-common :refer [test-ref-seq-file
                                     defslowtest]]))

(deftest overlap-exon-intron-boundary?-test
  (let [exon-ranges [[123 321] [456 654] [789 987]]]
    (are [p pos ref alt] (p (#'prot/overlap-exon-intron-boundary? exon-ranges pos ref alt))
      false? 454 "X" "Y"
      true? 454 "XXX" "X"
      false? 455 "XX" "X"
      false? 456 "XX" "X"
      false? 654 "XX" "X"
      true? 653 "XXX" "X"
      false? 653 "XX" "X")))

(deftest alt-exon-ranges-test
  ;; 1 [2 3 4] 5 6 7 [8 9 10 11] 12 13 14 15
  (are [p r a e] (= (#'prot/alt-exon-ranges [[2 4] [8 11]] p r a) e)
    3 "X" "XXX" [[2 6] [10 13]]
    6 "X" "XXX" [[2 4] [10 13]]
    2 "XX" "X" [[2 3] [7 10]]
    3 "XX" "X" [[2 3] [7 10]]
    6 "XX" "X" [[2 4] [7 10]]
    9 "XXX" "YYY" [[2 4] [8 11]])
  ;; Variants overlapping a boundary of exon/intron
  (are [p r a] (nil? (#'prot/alt-exon-ranges [[2 4] [8 11]] p r a))
    3 "XXX" "YYY"
    6 "XXX" "X"
    3 "XXX" "X"
    1 "XXXXX" "X"))

(deftest exon-sequence-test
  ;; A  C G T  A C G  T A C  G   T  A  C  G
  ;; 1 [2 3 4] 5 6 7 [8 9 10 11] 12 13 14 15
  (is (= (#'prot/exon-sequence "ACGTACGTACGTACG" 1 [[2 4] [8 11]]) "CGTTACG")))

(deftest is-insertion-variant?-test
  (are [p ref alt] (p (#'prot/is-insertion-variant? ref alt))
    false? "T" "A" ; substitution
    false? "TAGTCTA" "T" ; deletion
    true? "T" "TGTGATC" ; insertion
    false? "AT" "AGTCATCC" ; indel
    false? "GATC" "GCATGCAT" ; indel
    ))

(deftest cds-start-upstream-to-cds-variant?-test
  (are [p cds-start pos ref] (p (#'prot/cds-start-upstream-to-cds-variant? cds-start pos ref))
    true? 100 99 "TCGA"
    false? 100 100 "CGTA"
    false? 100 100 "G"
    false? 100 99 "A"
    false? 100 95 "GATGC"
    false? 100 101 "GTCAT"))

(deftest cds-to-cds-end-downstream-variant?-test
  (are [p cds-end pos ref] (p (#'prot/cds-to-cds-end-downstream-variant? cds-end pos ref))
    true? 100 99 "TCGA"
    true? 100 100 "CGTA"
    false? 100 100 "A"
    false? 100 95 "GATGC"
    false? 100 101 "GTCAT"))

(deftest make-alt-up-exon-seq-test
  (let [alt-up-seq "AATGCTTCTAGCTCCATGGCTATCGGC"
        start 3
        exon-ranges [[3 10] [13 18] [23 28] [35 45]]]
    (are [p end strand] (= p (#'prot/make-alt-up-exon-seq alt-up-seq
                                                          start
                                                          end
                                                          exon-ranges
                                                          strand))
      "TGCTTCGCTCCATAT"    25 :forward
      "AATGCTTCGCTCCATATC" 26 :forward
      "AATGCTTCGCTCCATAT"  25 :reverse
      "AATGCTTCGCTCCATATC" 26 :reverse)))

(deftest make-alt-down-exon-seq-test
  (let [alt-down-seq "AATGCTTCTAGCTCCATGGCTATCGGC"
        end 50
        exon-ranges [[3 10] [13 18] [23 28] [35 45]]]
    (are [p start strand] (= p (#'prot/make-alt-down-exon-seq alt-down-seq
                                                              start
                                                              end
                                                              exon-ranges
                                                              strand))
      "AATGCCTCCATGGCTA" 24 :forward
      "AATGGCTCCATGGCT"  25 :forward
      "AATGCCTCCATGGCT"  24 :reverse
      "AATGGCTCCATGGCT"  25 :reverse)))

(deftest make-ter-site-adjusted-alt-seq-test
  (let [alt-seq "XXXXXX"
        upstream-seq "YYYYYY"
        downstream-seq "ZZZZZZ"
        [cds-start cds-end] [7 12]]
    (are [p strand pos ref ref-include-ter-site] (= p (#'prot/make-ter-site-adjusted-alt-seq alt-seq
                                                                                             upstream-seq
                                                                                             downstream-seq
                                                                                             strand
                                                                                             cds-start
                                                                                             cds-end
                                                                                             pos
                                                                                             ref
                                                                                             ref-include-ter-site))
      "XXXXXX" :forward 8 "XX" false
      "XXXXXX" :forward 5 "YY" false
      "XXXXXX" :forward 5 "YYX" false
      "XXXXXX" :forward 13 "ZZ" false
      "XXXXXXZZZZZZ" :forward 10 "XX" true
      "XXXXXXZZZZZZ" :forward 12 "XZZ" true
      "YYYYYYXXXXXX" :reverse 8 "XX" true
      "YYYYYYXXXXXX" :reverse 9 "XX" true
      "XXXXXX" :reverse 5 "YY" false
      "YYYYYYXXXXXX" :reverse 5 "YYX" true
      "XXXXXX" :reverse 13 "ZZ" false
      "XXXXXX" :reverse 12 "XZZ" false)))

(deftest include-utr-ini-site-boundary?-test
  (let [forward-rg {:strand :forward
                    :cds-start 10
                    :cds-end 20}
        reverse-rg {:strand :reverse
                    :cds-start 10
                    :cds-end 20}]
    (are [p rg pos ref alt] (p (#'prot/include-utr-ini-site-boundary? rg pos ref alt))
      true? forward-rg 8 "CCAT" "C"
      true? forward-rg 8 "CCAT" "CGGG"
      false? forward-rg 9 "CAT" "C"
      true? reverse-rg 18 "CATGG" "C"
      true? reverse-rg 19 "ATGG" "AACC"
      false? reverse-rg 20 "TGG" "T")))

(deftest include-ter-site?-test
  (let [forward-rg {:strand :forward
                    :cds-start 10
                    :cds-end 20}
        reverse-rg {:strand :reverse
                    :cds-start 10
                    :cds-end 20}]
    (are [p rg pos ref alt] (p (#'prot/include-ter-site? rg pos ref alt))
      true? forward-rg 17 "ATG" "A"
      true? forward-rg 19 "GAG" "GCC"
      false? forward-rg 20 "AAT" "A"
      true? reverse-rg 9 "TT" "T"
      true? reverse-rg 9 "AACT" "AG"
      false? reverse-rg 8 "AG" "A")))

(deftest ter-site-same-pos?-test
  (are [p ref alt] (p (#'prot/ter-site-same-pos? ref alt))
    true? "MTGA*" "MTGA*"
    true? "MTGA*" "MTGA*CT"
    false? "MTGA*" "MTGAQCT*"
    false? "MTGA*" "MTGA"))

(deftest utr-variant?-test
  (let [cds-start 10
        cds-end 21]
    (are [p pos ref alt] (p (#'prot/utr-variant? cds-start cds-end pos ref alt))
      ;; cds-start upstream
      true? 9 "G" "T"
      true? 9 "G" "GA"
      true? 8 "GT" "G"
      true? 8 "TA" "TTCG"
      true? 7 "CGT" "CAGA"
      false? 10 "A" "T"
      false? 10 "A" "AT"
      false? 10 "ATG" "A"
      false? 9 "TA" "TTCG"
      false? 8 "CGAT" "CAGA"

      ;; cds-end downstream
      true? 22 "G" "T"
      true? 21 "A" "AT"
      true? 21 "AT" "A"
      true? 21 "AG" "AATC"
      true? 21 "CGTC" "CAGA"
      false? 21 "A" "T"
      false? 20 "A" "AG"
      false? 20 "AA" "A"
      false? 20 "GA" "GGCT"
      false? 20 "TATA" "TCG")))

(deftest get-alt-cds-start-pos-test
  (let [cds-start 40
        exon-ranges [[10 50] [80 120] [150 200]]
        pos* 40]
    (are [pos-start pos-end p] (= (#'prot/get-alt-cds-start-pos cds-start pos-start pos-end exon-ranges pos*) p)
      40 40 40
      35 45 46
      35 55 80)))

(deftest get-alt-cds-end-pos-test
  (let [cds-end 100
        exon-ranges [[10 50] [80 120] [150 200]]
        pos* 100]
    (are [pos-start pos-end p] (= (#'prot/get-alt-cds-end-pos cds-end pos-start pos-end exon-ranges pos*) p)
      100 100 100
      100 115 99
      75  110 50)))

(deftest apply-offset-test
  (testing "ref not include exon terminal"
    (let [ref "GCTGACC"
          alt "G"
          cds-start 40
          cds-end 110
          exon-ranges [[10 50] [80 120] [150 200]]]
      (are [pos pos* p] (= (#'prot/apply-offset pos ref alt cds-start cds-end exon-ranges pos*) p)
        ;; position is between cds-start and cds-end
        100 40 40
        100 110 104
        100 200 194
        ;; position is upstream of cds-start
        20 40 34
        20 110 104
        20 200 194
        ;; position is downstream of cds-end
        160 40 40
        160 110 110
        160 200 194)))
  (testing "ref includes exon terminal"
    (let [ref "GCTGACC"
          alt "G"
          exon-ranges [[10 50] [80 120] [150 200]]]
      (are [pos cds-start cds-end pos* p] (= (#'prot/apply-offset pos ref alt cds-start cds-end exon-ranges pos*) p)
        ;; variant around cds-start: 45_50del
        44 45 160 45 74
        ;; variant around cds-end: 150_155del
        149 45 155 155 120))))

(deftest get-first-diff-aa-info-test
  (let [ref-seq "ABCDEFGHIJKLMN"]
    (are [p alt-seq pos] (= (#'prot/get-first-diff-aa-info pos
                                                           ref-seq
                                                           alt-seq)
                            p)
      {:ppos 6 :pref "F" :palt "G"} "ABCDEG" 4
      {:ppos 13 :pref "M" :palt "K"} "ACBDEFGHIJKLK" 10)))

(deftest first-diff-aa-is-ter-site?-test
  (let [ref-seq "ABCDEFGHIJKLM*"]
    (are [pred alt-seq pos] (pred (#'prot/first-diff-aa-is-ter-site? pos
                                                                     ref-seq
                                                                     alt-seq))
      false? "ABCDEG" 1
      true? "ABCDEFGHIJKLMNO" 1
      true? "ABCDEFGHIJKLMNO*" 1)))

(def ref-gene-EGFR
  {:bin 125
   :name "NM_005228"
   :chr "chr7"
   :strand :forward
   :tx-start 55019032
   :tx-end 55207338
   :cds-start 55019278
   :cds-end 55205617
   :exon-count 28
   :exon-ranges [[55019032 55019365] [55142286 55142437] [55143305 55143488] [55146606 55146740] [55151294 55151362] [55152546 55152664] [55154011 55154152] [55155830 55155946] [55156533 55156659] [55156759 55156832] [55157663 55157753] [55160139 55160338] [55161499 55161631] [55163733 55163823] [55165280 55165437] [55171175 55171213] [55172983 55173124] [55173921 55174043] [55174722 55174820] [55181293 55181478] [55191719 55191874] [55192766 55192841] [55198717 55198863] [55200316 55200413] [55201188 55201355] [55201735 55201782] [55202517 55202625] [55205256 55207338]]
   :score 0
   :name2 "EGFR"
   :cds-start-stat :cmpl
   :cds-end-stat :cmpl,
   :exon-frames [0 1 0 1 1 1 0 1 1 2 1 2 1 2 0 2 2 0 0 0 0 0 1 1 0 0 0 1]})

(deftest protein-position-test
  (are [pos ppos] (= (#'prot/protein-position pos ref-gene-EGFR) ppos)
    ;; exon
    55191822 858
    55181378 790

    ;; intron
    55146585 142

    ;; CDS edges
    (:cds-start ref-gene-EGFR) 1
    (:cds-end ref-gene-EGFR) 1211

    ;; outside of CDS
    (dec (:cds-start ref-gene-EGFR)) 1
    (inc (:cds-end ref-gene-EGFR)) 1211
    (:tx-start ref-gene-EGFR) 1
    (:tx-end ref-gene-EGFR) 1211))

(deftest prot-seq-pstring-test
  (are [pref-seq palt-seq start end m e]
       (= (#'prot/prot-seq-pstring pref-seq palt-seq start end m) e)
    "LAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPI"
    "LAARNVLVKTPQHVKITDFGRAKLLGAEEKEYHAEGGKVPI"
    838 878 {:ppos 858, :pref "L", :palt "R"}
    (string/join \newline ["   841       851       861       871"
                           "LAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPI"
                           "LAARNVLVKTPQHVKITDFGRAKLLGAEEKEYHAEGGKVPI"
                           "                    ^"])

    "MTILTYPFKNLPTASKWALRFS"
    "MTTILTYPFKNLPTASKWALRFS"
    1 22 {:ppos 2, :pref "T", :palt "TT"}
    (string/join \newline ["1          11        21"
                           "MT ILTYPFKNLPTASKWALRFS"
                           "MTTILTYPFKNLPTASKWALRFS"
                           " ^^"])

    "KEGHQEGLVELPASFRELLTFFCTNATIHGAIRLVCSRGNR"
    "KEGHQEGLVELPASFRELLTFCTNATIHGAIRLVCSRGNRL"
    206 246 {:ppos 226, :pref "FF", :palt "F"}
    (string/join \newline ["     211       221       231       241"
                           "KEGHQEGLVELPASFRELLTFFCTNATIHGAIRLVCSRGNR"
                           "KEGHQEGLVELPASFRELLTF CTNATIHGAIRLVCSRGNR"
                           "                    ^^"])

    "YDIGGPDQEFGVDVGPVCFL*"
    "YDIGGPDQEFGVDVGPVCFLQ"
    1447 1467 {:ppos 1467, :pref "*", :palt "Q"}
    (string/join \newline ["    1451      1461"
                           "YDIGGPDQEFGVDVGPVCFL*"
                           "YDIGGPDQEFGVDVGPVCFLQ"
                           "                    ^"])))

(defslowtest protein-seq-3'-rule-test
  (let [tp53 {:name2 "TP53"
              :name "NM_000546.6"
              :chr "chr17"
              :tx-start 7668421
              :tx-end 7687490
              :cds-start 7669609
              :cds-end 7676594
              :strand :reverse
              :cds-start-stat :cmpl
              :cds-end-stat :cmpl
              :exon-ranges [[7668421 7669690] [7670609 7670715] [7673535 7673608] [7673701 7673837] [7674181 7674290] [7674859 7674971] [7675053 7675236] [7675994 7676272] [7676382 7676403] [7676521 7676622] [7687377 7687490]]
              :bin 643
              :exon-frames [2 0 1 2 0 1 0 0 2 0 -1]
              :exon-count 11}]
    (are [pos ref alt res]
         (= (with-open [seq-rdr (cseq/reader test-ref-seq-file)] (#'prot/mutation seq-rdr tp53 pos ref alt {})) res)
      7676197 "G" "GGTCTTGTCCCTTA" (:mutation (hgvs/parse "p.P58*"))
      7676202 "T" "TGTCCCTTAGTCTT" (:mutation (hgvs/parse "p.P58*")))))
