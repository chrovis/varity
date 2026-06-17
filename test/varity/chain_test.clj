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

(deftest search-chains-test
  (testing "simple test"
    (let [chain [{:header
                  {:t-name "chr1" :q-name "chr1"
                   :t-start 0 :q-start 0
                   :t-end 26 :q-end 18
                   :tsize 26 :q-size 18}
                  :data
                  [{:size 10 :dt 2 :dq 0}
                   {:size 1 :dt 0 :dq 2}
                   {:size 1 :dt 4 :dq 0}
                   {:size 1 :dt 4 :dq 0}
                   {:size 3}]}]
          chain-idx (chain/index chain)]
      (are [pos ans] (= (-> (chain/search-chains "chr1" pos chain-idx)
                            first
                            :in-block
                            (select-keys [:q-start :t-start])
                            not-empty)
                        ans)
        9 {:t-start 1 :q-start 1}
        10 {:t-start 1 :q-start 1}
        11 nil
        12 nil
        13 {:t-start 13 :q-start 11} ;;gap=+2
        14 {:t-start 14 :q-start 14} ;;gap=+2-2
        15  nil
        16  nil
        17  nil
        18  nil
        19 {:t-start 19 :q-start 15} ;;gap=+2-2+4
        20  nil
        21  nil
        22  nil
        23  nil
        24  {:t-start 24 :q-start 16} ;;gap-+2-2+4+4
        25  {:t-start 24 :q-start 16}
        26  {:t-start 24 :q-start 16})))
  (testing "size=0 test"
    (let [chain [{:header
                  {:t-name "chr1" :q-name "chr1"
                   :t-start 0 :q-start 0
                   :t-end 21 :q-end 12
                   :tsize 21 :q-size 12}
                  :data
                  [{:size 0 :dt 4 :dq 0} ;; [1,4]
                   {:size 6 :dt 0 :dq 0} ;; [5,10]
                   {:size 0 :dt 4 :dq 0} ;; [11,14]
                   {:size 6}]}] ;;[15,20]
          chain-idx (chain/index chain)]
      (is
       (= (map #(-> (chain/search-chains "chr1" % chain-idx)
                    first
                    :in-block
                    some?)
               (range 1 21))
          [false false false false ;; [1,4]
           true true true true true true ;; [5,10]
           false false false false ;; [11, 14]
           true true true true true true]))))) ;; [15,20]

(deftest search-overlap-blocks-test
  (are [start end ans] (= (map #(select-keys % [:t-start :t-end])
                               (chain/search-overlap-blocks start end sample-indexed-block))
                          ans)
    7 12 [{:t-start 6, :t-end 10} {:t-start 8, :t-end 12}]
    15 17 [{:t-start 14, :t-end 16}]
    6 7 [{:t-start 6, :t-end 10}]
    4 4 [{:t-start 2, :t-end 5}]
    5 5 [{:t-start 2, :t-end 5}]
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
