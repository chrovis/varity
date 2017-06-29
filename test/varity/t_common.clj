(ns varity.t-common
  (:require [clojure.java.io :as io]
            [clojure.test :refer [deftest testing]]
            [cavia.core :as cavia :refer [defprofile with-profile]]))

(defn- in-cloverage? []
  (some? (resolve 'cloverage.coverage/*covered*)))

(defmacro defslowtest
  [name & body]
  (if-not (in-cloverage?)
    `(deftest ~(vary-meta name assoc :slow true)
       ~@body)))

(defprofile prof
  {:resources [{:id "test.2bit"
                :url "https://test.chrov.is/data/varity/hg38.2bit"
                :sha1 "6fb20ba4de0b49247b78e08c2394d0c4f8594148"}
               {:id "test-refGene.txt.gz"
                :url "https://test.chrov.is/data/varity/hg38-refGene.txt.gz"
                :sha1 "941d514e57f4e842743f5c9269a0279906a072a0"}]})

(defn prepare-cavia! []
  (with-profile prof
    (cavia/without-print (cavia/get!))))

(defmacro cavia-testing
  [s & body]
  `(testing ~s
     (prepare-cavia!)
     ~@body))

(def test-ref-seq-file (cavia/resource prof "test.2bit"))
(def test-ref-gene-file (cavia/resource prof "test-refGene.txt.gz"))

(def test-chain-file (io/resource "hg19ToHg38.over.chain.gz"))
