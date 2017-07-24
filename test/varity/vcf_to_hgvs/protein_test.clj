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

;; TODO add test for #'varity.vcf-to-hgvs.protein/exon-sequence/protein-position
