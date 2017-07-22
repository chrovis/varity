(ns varity.vcf-to-hgvs.common-test
  (:require [clojure.test :refer :all]
            [clojure.java.io :as io]
            [cljam.io.sequence :as cseq]
            [varity.vcf-to-hgvs.common :refer :all]))

(deftest diff-bases-test
  (testing "diff-bases returns a vector of diff info"
    (are [s1 s2 ret] (= (diff-bases s1 s2) ret)
      "C"   "G"       ["C" "G" 0 0]
      "TC"  "T"       ["C" "" 1 0]
      "C"   "CA"      ["" "A" 1 0]
      "TCG" "TCAG"    ["" "A" 2 1]
      "T"   "TTTTCTT" ["" "TTTCTT" 1 0]
      "CAG" "CTC"     ["AG" "TC" 1 0]
      "CA"  "CTC"     ["A"  "TC" 1 0]
      "CAG" "CT"      ["AG" "T" 1 0]
      "CAGTC" "CTGAC" ["AGT" "TGA" 1 1])))

(deftest forward-shift-test
  (testing "forward-shift calculates forward shifting steps"
    (are [s p b ret] (= (forward-shift s p b) ret)
      "AGGGGGGT" 2 "G" 5
      "AGGGGGGT" 2 "GG" 4
      "AGGGGGGT" 2 "GGG" 3
      "AAGTGTCC" 3 "GT" 2
      "AAGTGTCC" 5 "GT" 0))
  (testing "foward-shift throws exception when inputs are illegal"
    (are [s p b] (thrown? Exception (forward-shift s p b))
      "AGGGGGGT" 2 "T"
      "" 2 "T")))

(deftest backward-shift-test
  (testing "backward-shift calculates backward shifting steps"
    (are [s p b ret] (= (backward-shift s p b) ret)
      "AGGGGGGT" 7 "G" 5
      "AGGGGGGT" 6 "GG" 4
      "AGGGGGGT" 5 "GGG" 3
      "AAGTGTCC" 5 "GT" 2
      "AAGTGTCC" 3 "GT" 0))
  (testing "backward-shift throws exception when inputs are illegal"
    (are [s p b] (thrown? Exception (backward-shift s p b))
      "AGGGGGGT" 7 "A"
      "" 7 "A")))

(deftest repeat-info-test
  (are [s p i e] (= (repeat-info s p i) e)
    "XXXCAGTCXXX" 8 "AGT" ["AGT" 1 1]
    "XXXCAGTCXXX" 8 "AGTAGT" ["AGT" 1 2]
    "XXXCAGTAGTCXXX" 11 "AGTAGT" ["AGT" 2 2]))

(def normalize-variant* #'varity.vcf-to-hgvs.common/normalize-variant*)

(deftest normalize-variant*-test
  (testing "normalize-variant* normalizes variant"
    (are [v st ret] (= (normalize-variant* v "NNNCAGTAGTAGTCNNN" st) ret)
      {:pos 7, :ref "T", :alt "TAGT"} "+" {:pos 13, :ref "T", :alt "TAGT"}
      {:pos 7, :ref "TAGT", :alt "T"} "+" {:pos 10, :ref "TAGT", :alt "T"}
      {:pos 7, :ref "T", :alt "TAGT"} "-" {:pos 4, :ref "C", :alt "CAGT"}
      {:pos 7, :ref "TAGT", :alt "T"} "-" {:pos 4, :ref "CAGT", :alt "C"})))


(def rg1 {:strand "+" :tx-start 5 :tx-end 42})
(def rg2 {:strand "-" :tx-start 3 :tx-end 32})

(deftest normalize-variant-test
  (testing "normalize-variant"
    (with-open [seq-rdr (cseq/reader (.getAbsolutePath (io/file (io/resource "test.fa"))))]
      (are [v rg ret] (= (normalize-variant v seq-rdr rg) ret)
        {:chr "ref1" :pos 38 :ref "AGCGC" :alt "A"} rg1 {:chr "ref1" :pos 38 :ref "AGCGC" :alt "A"}))))
