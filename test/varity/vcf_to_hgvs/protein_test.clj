(ns varity.vcf-to-hgvs.protein-test
  (:require [clojure.string :as string]
            [clojure.test :refer :all]
            [cljam.io.sequence :as cseq]
            [clj-hgvs.core :as hgvs]
            [varity.vcf-to-hgvs.protein :as prot]
            [varity.t-common :refer [test-ref-seq-file
                                     defslowtest]]))

(deftest alt-exon-ranges-test
  ;; 1 [2 3 4] 5 6 7 [8 9 10 11] 12 13 14 15
  (are [p r a e] (= (#'prot/alt-exon-ranges [[2 4] [8 11]] p r a) e)
    3 "X" "XXX" [[2 6] [10 13]]
    6 "X" "XXX" [[2 4] [10 13]]
    2 "XX" "X" [[2 3] [7 10]]
    3 "XX" "X" [[2 3] [7 10]]
    6 "XX" "X" [[2 4] [7 10]]
    6 "XXX" "X" [[2 4] [7 9]]
    3 "XXX" "X" [[2 3] [6 9]]
    1 "XXXXX" "X" [[4 7]]
    9 "XXX" "XXX" [[2 4] [8 11]])
  ;; Can't determine whether the splice site is shifted or not
  (is (thrown-with-msg?
       Exception
       #"unsupported"
       (#'prot/alt-exon-ranges [[2 4] [8 11]] 3 "XXX" "XXX"))))

(deftest exon-sequence-test
  ;; A  C G T  A C G  T A C  G   T  A  C  G
  ;; 1 [2 3 4] 5 6 7 [8 9 10 11] 12 13 14 15
  (is (= (#'prot/exon-sequence "ACGTACGTACGTACG" 1 [[2 4] [8 11]]) "CGTTACG")))

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
