(ns varity.vcf-to-hgvs.protein-test
  (:require [clojure.test :refer :all]
            [varity.vcf-to-hgvs.protein :refer :all]))

(def alt-sequence #'varity.vcf-to-hgvs.protein/alt-sequence)

(deftest alt-sequence-test
  (are [st p r a e] (= (alt-sequence "ACGTACGTACGT" st p r a) e)
    5 8 "T" "A" "ACGAACGTACGT"
    5 8 "T" "TT" "ACGTTACGTACGT"
    5 8 "TA" "T" "ACGTCGTACGT"))

(def alt-exon-ranges #'varity.vcf-to-hgvs.protein/alt-exon-ranges)

(deftest alt-exon-ranges-test
  ;; 1 [2 3 4] 5 6 7 [8 9 10 11] 12 13 14 15
  (are [p r a e] (= (alt-exon-ranges [[2 4] [8 11]] p r a) e)
    3 "XXX" "XXX" [[2 4] [8 11]]
    3 "X" "XXX" [[2 6] [10 13]]
    6 "X" "XXX" [[2 4] [10 13]]
    3 "XX" "X" [[2 3] [7 10]]
    6 "XX" "X" [[2 4] [7 10]]
    6 "XXX" "X" [[2 4] [7 9]]
    3 "XXX" "X" [[2 3] [6 9]]
    1 "XXXXX" "X" [[4 7]]))

(def exon-sequence #'varity.vcf-to-hgvs.protein/exon-sequence)

(deftest exon-sequence-test
  ;; A  C G T  A C G  T A C  G   T  A  C  G
  ;; 1 [2 3 4] 5 6 7 [8 9 10 11] 12 13 14 15
  (is (= (exon-sequence "ACGTACGTACGTACG" 1 [[2 4] [8 11]]) "CGTTACG")))

(def protein-position #'varity.vcf-to-hgvs.protein/protein-position)

(def ref-gene-EGFR {:name2 "EGFR", :cds-start-stat :cmpl, :exon-ranges [[55019032 55019365] [55142286 55142437] [55143305 55143488] [55146606 55146740] [55151294 55151362] [55152546 55152664] [55154011 55154152] [55155830 55155946] [55156533 55156659] [55156759 55156832] [55157663 55157753] [55160139 55160338] [55161499 55161631] [55163733 55163823] [55165280 55165437] [55171175 55171213] [55172983 55173124] [55173921 55174043] [55174722 55174820] [55181293 55181478] [55191719 55191874] [55192766 55192841] [55198717 55198863] [55200316 55200413] [55201188 55201355] [55201735 55201782] [55202517 55202625] [55205256 55207338]], :cds-end-stat :cmpl, :tx-start 55019032, :name "NM_005228", :strand :forward, :cds-start 55019278, :score "0", :tx-end 55207338, :bin 125, :exon-frames "0,1,0,1,1,1,0,1,1,2,1,2,1,2,0,2,2,0,0,0,0,0,1,1,0,0,0,1,", :exon-count 28, :chr "chr7", :cds-end 55205617})

(deftest protein-position-test
  (are [pos ppos] (= (protein-position pos ref-gene-EGFR) ppos)
    55191822 858 ; exon
    55181378 790 ; exon
    55146585 nil ; intron
    (:tx-start ref-gene-EGFR) nil
    (:tx-end ref-gene-EGFR) nil
    (:cds-start ref-gene-EGFR) 1
    (:cds-end ref-gene-EGFR) 1211
    (dec (:cds-start ref-gene-EGFR)) nil
    (inc (:cds-end ref-gene-EGFR)) nil))
