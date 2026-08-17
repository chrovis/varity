(ns varity.chain-test
  (:require [clojure.test :refer [deftest is are testing]]
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

(def ^:private sample-indexed-chain
  (chain/index
   [{:header {:t-name "chr1", :t-start 0, :t-end 100, :score 1, :id 1,
              :q-name "1", :q-start 0, :q-end 100}, :data [{:size 100}]}
    {:header {:t-name "chr1", :t-start 100, :t-end 200, :score 2, :id 2,
              :q-name "1", :q-start 100, :q-end 200}, :data [{:size 100}]}
    {:header {:t-name "chr1", :t-start 200, :t-end 300, :score 3, :id 3,
              :q-name "1", :q-start 200, :q-end 300}, :data [{:size 100}]}
    {:header {:t-name "chr1", :t-start 250, :t-end 350, :score 4, :id 4,
              :q-name "1", :q-start 250, :q-end 350}, :data [{:size 100}]}
    {:header {:t-name "chr2", :t-start 0, :t-end 300, :score 1, :id 5,
              :q-name "2", :q-start 0, :q-end 300}, :data [{:size 300}]}
    {:header {:t-name "chr3", :t-start 0, :t-end 300, :score 1, :id 6,
              :q-name "3", :q-start 0, :q-end 300}, :data [{:size 300}]}]))

(deftest search-chains-test
  (testing "single chain"
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
  (testing "multiple chains"
    (letfn [(overlapping-ids [chr pos]
              (->> sample-indexed-chain
                   (chain/search-chains chr pos)
                   (map (comp :id :header))
                   sort))]
      (is (= (overlapping-ids "chr1"   1) [1]))
      (is (= (overlapping-ids "chr1" 100) [1]))
      (is (= (overlapping-ids "chr1" 101) [2]))
      (is (= (overlapping-ids "chr1" 101) [2]))
      (is (= (overlapping-ids "chr1" 200) [2]))
      (is (= (overlapping-ids "chr1" 201) [3]))
      (is (= (overlapping-ids "chr1" 250) [3]))
      (is (= (overlapping-ids "chr1" 251) [3 4]))
      (is (= (overlapping-ids "chr1" 300) [3 4]))
      (is (= (overlapping-ids "chr1" 301) [4]))
      (is (= (overlapping-ids "chr1" 350) [4]))
      (is (= (overlapping-ids "chr1" 351) [])))))

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
  (letfn [(containing-ids [chr start end]
            (->> sample-indexed-chain
                 (chain/search-containing-chains chr start end)
                 (map (comp :id :header))))]
    (is (= (containing-ids "chr1"  10  20) [1]))
    (is (= (containing-ids "chr1"  99  99) [1]))
    (is (= (containing-ids "chr1" 100 100) [1]))
    (is (= (containing-ids "chr1" 100 101) []))
    (is (= (containing-ids "chr1" 101 101) [2]))
    (is (= (containing-ids "chr1" 200 200) [2]))
    (is (= (containing-ids "chr1" 200 201) []))
    (is (= (containing-ids "chr1" 201 201) [3]))
    (is (= (containing-ids "chr1" 250 250) [3]))
    (is (= (containing-ids "chr1" 250 251) [3]))
    (is (= (containing-ids "chr1" 251 251) [4 3]))
    (is (= (containing-ids "chr1" 260 270) [4 3]))
    (is (= (containing-ids "chr1"  10 120) []))
    (is (= (containing-ids "chr2"  10  20) [5]))))
