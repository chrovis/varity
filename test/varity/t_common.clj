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
  {:resources [{:id "test.fa"
                :url "https://test.chrov.is/data/varity/hg38.fa.gz"
                :md5 "b2aee9f885accc00531e59c4736bee63"
                :packed :gzip}
               {:id "test.fa.fai"
                :url "https://test.chrov.is/data/varity/hg38.fa.fai"
                :sha1 "cfd5f7326d2e23e4e13fc2e2da21f7086b47e545"}
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

(def test-fa-file (cavia/resource prof "test.fa"))
(def test-ref-gene-file (cavia/resource prof "test-refGene.txt.gz"))

(def test-chain-file (io/resource "hg19ToHg38.over.chain.gz"))
