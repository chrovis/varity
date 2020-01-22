(ns varity.vcf-to-hgvs.common-test
  (:require [clojure.test :refer :all]
            [cljam.io.sequence :as cseq]
            [varity.vcf-to-hgvs.common :refer :all]
            [varity.t-common :refer :all]
            [varity.vcf-to-hgvs.common :as common]))

(deftest diff-bases-test
  (testing "diff-bases returns a vector of diff info"
    (are [s1 s2 ret] (= (common/diff-bases s1 s2) ret)
      "C"   "G"       ["C" "G" 0 0]
      "TC"  "T"       ["C" "" 1 0]
      "C"   "CA"      ["" "A" 1 0]
      "TCG" "TCAG"    ["" "A" 2 1]
      "T"   "TTTTCTT" ["" "TTTCTT" 1 0]
      "CAG" "CTC"     ["AG" "TC" 1 0]
      "CA"  "CTC"     ["A"  "TC" 1 0]
      "CAG" "CT"      ["AG" "T" 1 0]
      "CAGTC" "CTGAC" ["AGT" "TGA" 1 1])))

(deftest repeat-units-test
  (are [s ret] (= (#'common/repeat-units s) ret)
    "AGT"          ["AGT"]
    "AGTAGT"       ["AGT" "AGTAGT"]
    "AGTAGTAGTAGT" ["AGT" "AGTAGT" "AGTAGTAGTAGT"]
    "AGTTGA"       ["AGTTGA"]
    "A"            ["A"]
    ""             []
    nil            [])
  (let [s (apply str (repeat (* common/max-repeat-unit-size 10) "A"))]
    (is (<= (count (first (take-last 2 (#'common/repeat-units s))))
            common/max-repeat-unit-size))))

(deftest forward-shift-test
  (testing "forward-shift calculates forward shifting steps"
    (are [s p b ret] (= (common/forward-shift s p b) ret)
      "AGGGGGGT" 2 "G" 5
      "AGGGGGGT" 2 "GG" 4
      "AGGGGGGT" 2 "GGG" 3
      "AAGTGTCC" 3 "GT" 2
      "AAGTGTCC" 5 "GT" 0))
  (testing "foward-shift throws exception when inputs are illegal"
    (are [s p b] (thrown? Exception (common/forward-shift s p b))
      "AGGGGGGT" 2 "T"
      "" 2 "T")))

(deftest backward-shift-test
  (testing "backward-shift calculates backward shifting steps"
    (are [s p b ret] (= (common/backward-shift s p b) ret)
      "AGGGGGGT" 7 "G" 5
      "AGGGGGGT" 6 "GG" 4
      "AGGGGGGT" 5 "GGG" 3
      "AAGTGTCC" 5 "GT" 2
      "AAGTGTCC" 3 "GT" 0))
  (testing "backward-shift throws exception when inputs are illegal"
    (are [s p b] (thrown? Exception (common/backward-shift s p b))
      "AGGGGGGT" 7 "A"
      "" 7 "A")))

(deftest repeat-info-test
  (are [s p a t e] (= (common/repeat-info s p a t) e)
    "XXXCAGTCXXX"       8  "AGT"    :ins ["AGT" 1 2]
    "XXXCAGTCXXX"       8  "AGTAGT" :ins ["AGT" 1 3]
    "XXXCAGTAGTCXXX"    11 "AGTAGT" :ins ["AGT" 2 4]
    "XXXCAGTAGTAGTCXXX" 11 "AGT"    :del ["AGT" 3 2]
    "XXXCAGTAGTAGTCXXX" 8  "AGTAGT" :del ["AGT" 3 1])
  (are [s p a t] (nil? (common/repeat-info s p a t))
    "XXXCAGTCXXX"       8  "ATG" :ins
    "XXXCAGTAGTAGTCXXX" 11 "ATG" :del
    "XXXCAGTCXXX"       5  "ATG" :del))

(deftest normalize-variant*-test
  (testing "normalize-variant* normalizes variant"
    (are [v st ret] (= (#'common/normalize-variant* v "NNNCAGTAGTAGTCNNN" st) ret)
      {:pos 7, :ref "T", :alt "TAGT"} :forward {:pos 13, :ref "T", :alt "TAGT"}
      {:pos 7, :ref "TAGT", :alt "T"} :forward {:pos 10, :ref "TAGT", :alt "T"}
      {:pos 7, :ref "T", :alt "TAGT"} :reverse {:pos 4, :ref "C", :alt "CAGT"}
      {:pos 7, :ref "TAGT", :alt "T"} :reverse {:pos 4, :ref "CAGT", :alt "C"})))

(defslowtest normalize-variant-test
  (cavia-testing "normalize without error"
    (with-open [seq-rdr (cseq/reader test-ref-seq-file)]
      (are [v rg] (map? (common/normalize-variant v seq-rdr rg))
        {:chr "chr17", :pos 43125270, :ref "CCTTTACCCAGAGCAGAGGGTGAAGGCCTCCTGAGCGCAGGGGCCCAGTTATCTGAGAAACCCCACAGCCTGTCCCCCGTCCAGGAAGTCTCAGCGAGCTCACGCCGCGCAGTCGCAGTTTTAATTTATCTGTAATTCCCGCGCTTTTCCGTTGCCACGGAAACCAAGGGGCTACCGCTAAG", :alt "C"}
        {:tx-start 43044295, :strand :reverse}))))

(deftest alt-sequence-test
  (are [st p r a e] (= (common/alt-sequence "ACGTACGTACGT" st p r a) e)
    5 8 "T" "A" "ACGAACGTACGT"
    5 8 "T" "TT" "ACGTTACGTACGT"
    5 8 "TA" "T" "ACGTCGTACGT"))
