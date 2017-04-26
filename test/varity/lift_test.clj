(ns varity.lift-test
  (:require [clojure.test :refer :all]
            [varity.chain :as ch]
            [varity.lift :refer :all]
            [varity.t-common :refer :all]))

(deftest convert-coord-test
  (let [chidx (ch/index (ch/load-chain test-chain-file))]
    (testing "found"
      (are [c p r] (= (convert-coord {:chr c, :pos p} chidx) r)
        "chr1" 743267    {:chr "chr1", :pos 807887}
        "1"    743267    {:chr "chr1", :pos 807887}
        "chr1" 142613288 {:chr "chr4_GL000008v2_random", :pos 77854}
        "chrY" 1266342   {:chr "chrX", :pos 1198161}))
    (testing "found (strand -)"
      (are [c p r] (= (convert-coord {:chr c, :pos p} chidx) r)
        "chr1"  144274520 {:chr "chr1", :pos 149671178}
        "chr10" 47000000  {:chr "chr10", :pos 46549615}))
    (testing "not found"
      (are [c p] (nil? (convert-coord {:chr c, :pos p} chidx))
        "chr1" 300
        "chr10" 46478862))))
