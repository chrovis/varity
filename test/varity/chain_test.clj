(ns varity.chain-test
  (:require [clojure.test :refer :all]
            [cljam.util.intervals :as intervals]
            [varity.chain :as chain]))

(def ^:private sample-indexed-block
  (intervals/index-intervals
   (map (fn [{:keys [t-start t-end] :as m}]
          (assoc m :start t-start :end t-end))
        [{:t-start 2 :t-end 5}
         {:t-start 6 :t-end 10}
         {:t-start 8 :t-end 12}
         {:t-start 14 :t-end 16}
         {:t-start 18 :t-end 20}])
   {:structure :nclist}))

(def ^:private ^:const sample-indexed-chain
  {"chr1" [{:header {:t-name "chr1", :t-start 0, :t-end 100 :score 1}}
           {:header {:t-name "chr1", :t-start 100, :t-end 200 :score 2}}
           {:header {:t-name "chr1", :t-start 200, :t-end 300 :score 3}}
           {:header {:t-name "chr1", :t-start 250, :t-end 350 :score 4}}]
   "chr2" [{:header {:t-name "chr2", :t-start 0, :t-end 300 :score 1}}]
   "chr3" [{:header {:t-name "chr2", :t-start 0, :t-end 300 :score 1}}]})

(deftest search-overlap-blocks-test
  (are [start end ans] (= (map #(select-keys % [:t-start :t-end])
                               (chain/search-overlap-blocks start end sample-indexed-block))
                          ans)
    7 12 [{:t-start 6, :t-end 10} {:t-start 8, :t-end 12}]
    15 17 [{:t-start 14, :t-end 16}]
    6 7 [{:t-start 6, :t-end 10}]
    4 4 [{:t-start 2, :t-end 5}]
    5 5 []
    6 6 [{:t-start 6, :t-end 10}]
    1 1 []
    30 31 []))

(deftest search-containing-chains-test
  (is (= (chain/search-containing-chains "chr1" 10 20 sample-indexed-chain)
         [{:header {:t-name "chr1", :t-start 0, :t-end 100 :score 1}}]))
  (is (= (chain/search-containing-chains "chr1" 260 270 sample-indexed-chain)
         [{:header {:t-name "chr1", :t-start 250, :t-end 350 :score 4}}
          {:header {:t-name "chr1", :t-start 200, :t-end 300 :score 3}}]))
  (is (= (chain/search-containing-chains "chr1" 10 120 sample-indexed-chain)
         []))
  (is (= (chain/search-containing-chains "chr2" 10 20 sample-indexed-chain)
         [{:header {:t-name "chr2", :t-start 0, :t-end 300 :score 1}}])))
