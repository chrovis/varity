(ns varity.lift-bench
  (:require [libra.bench :refer :all]
            [libra.criterium :as c]
            [varity.chain :as ch]
            [varity.lift :as lift]
            [varity.t-common :refer :all]))

(def src-coords
  [;; found (+)
   {:chr "chr1", :pos 10001}
   {:chr "chr1", :pos 177377}
   {:chr "chr1", :pos 227418}
   {:chr "chr1", :pos 743267}
   {:chr "chr1", :pos 142613288}
   {:chr "chr1", :pos 143544525}
   {:chr "chr1", :pos 249240620}
   {:chr "chrY", :pos 1266342}
   ;; found (-)
   {:chr "chr1", :pos 317720}
   {:chr "chr1", :pos 471368}
   {:chr "chr1", :pos 144274520}
   {:chr "chr1", :pos 249060117}
   {:chr "chr1", :pos 249060118}
   {:chr "chr1", :pos 249060119}
   {:chr "chr1", :pos 249060188}
   {:chr "chr10", :pos 47000000}
   ;; not found
   {:chr "chr1", :pos 300}
   {:chr "chr1", :pos 10000}
   {:chr "chr1", :pos 177418}
   {:chr "chr1", :pos 267720}
   {:chr "chr1", :pos 471369}
   {:chr "chr1", :pos 143342816}
   {:chr "chr1", :pos 249240622}
   {:chr "chr10", :pos 46478862}])

(defbench convert-coord-bench
  (let [chidx (ch/index (ch/load-chain test-chain-file))]
    (is (c/quick-bench (doseq [coord src-coords]
                         (lift/convert-coord coord chidx))))))
