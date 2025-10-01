(ns varity.chain-test
  (:require [clojure.test :refer :all]
            [varity.chain :as chain]))

(def ^:private sample-indexed-block
  (sorted-map
   2 {:t-start 2 :t-end 5}
   6 {:t-start 6 :t-end 10}
   8 {:t-start 8 :t-end 12}
   14 {:t-start 14 :t-end 16}
   18 {:t-start 18 :t-end 20}))

(def ^:private ^:const sample-indexed-chain
  {"chr1" [{:header {:t-name "chr1", :t-start 0, :t-end 100 :score 1}}
           {:header {:t-name "chr1", :t-start 100, :t-end 200 :score 2}}
           {:header {:t-name "chr1", :t-start 200, :t-end 300 :score 3}}
           {:header {:t-name "chr1", :t-start 250, :t-end 350 :score 4}}]
   "chr2" [{:header {:t-name "chr2", :t-start 0, :t-end 300 :score 1}}]
   "chr3" [{:header {:t-name "chr2", :t-start 0, :t-end 300 :score 1}}]})

(deftest search-overlap-blocks-test
  (is (= (chain/search-overlap-blocks 7 12 sample-indexed-block)
         [{:t-start 6, :t-end 10} {:t-start 8, :t-end 12}]))
  (is (= (chain/search-overlap-blocks 15 17 sample-indexed-block)
         [{:t-start 14, :t-end 16}]))
  (is (= (chain/search-overlap-blocks 6 7 sample-indexed-block)
         [{:t-start 6, :t-end 10}]))
  (is (= (chain/search-overlap-blocks 4 4 sample-indexed-block)
         [{:t-start 2, :t-end 5}]))
  (is (= (chain/search-overlap-blocks 5 5 sample-indexed-block) []))
  (is (= (chain/search-overlap-blocks 6 6 sample-indexed-block)
         [{:t-start 6, :t-end 10}]))
  (is (= (chain/search-overlap-blocks 1 1 sample-indexed-block) []))
  (is (= (chain/search-overlap-blocks 30 31 sample-indexed-block) [])))

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
