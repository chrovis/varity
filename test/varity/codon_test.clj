(ns varity.codon-test
  (:require [clojure.test :refer :all]
            [varity.codon :refer :all]))

(deftest codon->amino-acid-test
  (are [s e] (= (codon->amino-acid s) e)
    "TCA" "S"
    "tca" "S"
    "TAA" "*")
  (are [s] (thrown? Error (codon->amino-acid s))
    "TC"
    "TCAG"
    ""
    nil))

(deftest amino-acid->codons-test
  (are [s e] (= (amino-acid->codons s) e)
    "L" ["CTG" "CTC" "CTT" "CTA" "TTG" "TTA"]
    \* ["TAG" "TAA" "TGA"])
  (are [s] (thrown? Error (amino-acid->codons s))
    "LS"
    ""
    nil))

(deftest amino-acid-sequence-test
  (are [s e] (= (amino-acid-sequence s) e)
    "ATGTCACTGTAA" "MSL*"
    "atgtcactgtaa" "MSL*"
    "" "")
  (is (thrown? Throwable (amino-acid-sequence nil))))
