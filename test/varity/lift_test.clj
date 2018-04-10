(ns varity.lift-test
  (:require [clojure.test :refer :all]
            [varity.chain :as ch]
            [varity.lift :refer :all]
            [varity.t-common :refer :all]))

(deftest convert-coord-test
  (let [chidx (ch/index (ch/load-chain test-chain-file))]
    (testing "found"
      (are [c p r] (= (convert-coord {:chr c, :pos p} chidx) r)
        "chr1"      10001 {:chr "chr1",     :pos 10001}
        "chr1"     177376 {:chr "chr1",    :pos 177376}
        "chr1"     177377 {:chr "chr19",   :pos 242824}
        "chr1"     177417 {:chr "chr19",   :pos 242864}
        "chr1"     227418 {:chr "chr1",    :pos 257667}
        "chr1"     267719 {:chr "chr1",    :pos 297968}
        "chr1"     743267 {:chr "chr1",    :pos 807887}
        "1"        743267 {:chr "chr1",    :pos 807887}
        "chr1"  142613288 {:chr "chr4_GL000008v2_random", :pos 77854}
        "chr1"  143342817 {:chr "chr14_GL000009v2_random", :pos 1}
        "chr1"  143544525 {:chr "chr14_GL000009v2_random", :pos 201709}
        "chr1"  249240620 {:chr "chr1", :pos 248946421}
        "chr1"  249240621 {:chr "chr1", :pos 248946422}
        "chrY"    1266342 {:chr "chrX",   :pos 1198161}))
    (testing "found (strand -)"
      (are [c p r] (= (convert-coord {:chr c, :pos p} chidx) r)
        "chr1"     317720 {:chr "chr1",    :pos 501617}
        "chr1"     471368 {:chr "chr1",    :pos 347969}
        "chr1"  144274520 {:chr "chr1", :pos 149671180}
        "chr1"  249060117 {:chr "chr3",  :pos 43941887}
        "chr1"  249060118 {:chr "chr3",  :pos 43941886}
        "chr1"  249060119 {:chr "chr13", :pos 81560696}
        "chr1"  249060188 {:chr "chr13", :pos 81560627}
        "chr10"  47000000 {:chr "chr10", :pos 46549617}))
    (testing "not found"
      (are [c p] (nil? (convert-coord {:chr c, :pos p} chidx))
        "chr1"        300
        "chr1"      10000
        "chr1"     177418
        "chr1"     227417
        "chr1"     267720
        "chr1"     317719
        "chr1"     471369
        "chr1"  143342816
        "chr1"  143544526
        "chr1"  249060189
        "chr1"  249240622
        "chr10"  46478862))))
