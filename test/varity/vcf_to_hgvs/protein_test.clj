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

(deftest is-deletion-variant?-test
  (are [p ref alt] (p (#'prot/is-deletion-variant? ref alt))
    false? "T" "A" ; substitution
    true? "TAGTCTA" "T" ; deletion
    false? "T" "TGTGATC" ; insertion
    true? "C" "GTCATCC" ; delins
    true? "ATC" "CATGCAT" ; delins
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
  (let [ref-up-exon-seq "AATGCTTCTAGCTCC"
        cds-start 100]
    (are [p pos ref alt] (= p (#'prot/make-alt-up-exon-seq ref-up-exon-seq
                                                           cds-start
                                                           pos
                                                           ref
                                                           alt))
      "AATGCTTCTAGCTCC" 102 "ATGTC" "A"
      "ATGCTTCTAGCT" 98 "CCTT" "C")))

(deftest make-alt-down-exon-seq-test
  (let [ref-up-exon-seq "CTTATAATAATAA"
        cds-end 1000]
    (are [p pos ref alt] (= p (#'prot/make-alt-down-exon-seq ref-up-exon-seq
                                                             cds-end
                                                             pos
                                                             ref
                                                             alt))
      "CTTATAATAATA" 1002 "TTATAA" "T"
      "TATAATAAT" 998 "GGCCT" "G")))

(deftest make-ter-site-adjusted-alt-seq-test
  (let [alt-seq "XXXXXX"
        upstream-seq "YYYYYY"
        downstream-seq "ZZZZZZ"
        [cds-start cds-end] [7 12]]
    (are [p strand pos ref] (#'prot/make-ter-site-adjusted-alt-seq alt-seq
                                                                   upstream-seq
                                                                   downstream-seq
                                                                   strand
                                                                   cds-start
                                                                   cds-end
                                                                   pos
                                                                   ref)
      "XXXXXX" :forward 8 "XX"
      "XXXXXX" :forward 5 "YY"
      "XXXXXX" :forward 5 "YYX"
      "XXXXXX" :forward 13 "ZZ"
      "XXXXXXZZZZZZ" :forward 12 "XZZ"
      "XXXXXX" :reverse 8 "XX"
      "XXXXXX" :reverse 5 "YY"
      "YYYYYYXXXXXX" :reverse 5 "YYX"
      "XXXXXX" :reverse 13 "ZZ"
      "XXXXXX" :reverse 12 "XZZ")))

(deftest get-pos-exon-end-tuple-test
  (let [exon-ranges [[1 10] [15 20] [25 40]]]
    (are [p pos] (= (#'prot/get-pos-exon-end-tuple pos exon-ranges) p)
      [15 20] 15
      [5 10] 5)))

(deftest include-utr-ini-site-boundary?-test
  (let [forward-rg {:strand :forward
                    :cds-start 10
                    :cds-end 20}
        reverse-rg {:strand :reverse
                    :cds-start 10
                    :cds-end 20}]
    (are [p rg pos ref alt] (p (#'prot/include-utr-ini-site-boundary? rg pos ref alt))
      true? forward-rg 8 "CCAT" "C"
      true? forward-rg 9 "CAT" "GGG"
      false? forward-rg 9 "CAT" "C"
      true? reverse-rg 18 "CATGG" "C"
      true? reverse-rg 20 "TGG" "ACC"
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
      true? forward-rg 20 "AG" "CC"
      false? forward-rg 20 "AAT" "A"
      true? reverse-rg 9 "TT" "T"
      true? reverse-rg 10 "ACT" "G"
      false? reverse-rg 8 "AG" "A")))

(deftest ter-site-same-pos?-test
  (are [p ref alt] (p (#'prot/ter-site-same-pos? ref alt))
    true? "MTGA*" "MTGA*"
    true? "MTGA*" "MTGA*CT"
    false? "MTGA*" "MTGAQCT*"
    false? "MTGA*" "MTGA"))

(deftest apply-offset-test
  (let [pos 100
        ref "GCTGACC"
        alt "G"
        exon-ranges [[10 50] [80 120] [150 200]]]
    (are [pos* ref-include-ter-site p] (= (#'prot/apply-offset pos ref alt exon-ranges ref-include-ter-site pos*) p)
      40 false 40
      110 false 104
      105 true 101
      112 true 106)))

(deftest get-first-diff-aa-info-test
  (let [ref-seq "ABCDEFGHIJKLMN"]
    (are [p alt-seq pos] (= (#'prot/get-first-diff-aa-info pos
                                                           ref-seq
                                                           alt-seq)
                            p)
      {:ppos 6 :pref "F" :palt "G"} "ABCDEG" 4
      {:ppos 13 :pref "M" :palt "K"} "ACBDEFGHIJKLK" 10)))

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
