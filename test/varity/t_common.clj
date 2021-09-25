(ns varity.t-common
  (:require [clojure.java.io :as io]
            [clojure.tools.logging :refer [*logger-factory*]]
            [clojure.tools.logging.impl :refer [disabled-logger-factory]]
            [clojure.test :refer [assert-expr deftest do-report testing]]
            [cavia.core :as cavia :refer [defprofile with-profile]]))

(defmethod assert-expr 'thrown-with-error-type? [msg form]
  ;; (is (thrown-with-error-type? key expr))
  ;; Asserts that evaluating expr throws an exception of error type key.
  (let [key (nth form 1)
        body (nthnext form 2)]
    `(try ~@body
          (do-report {:type :fail, :message ~msg, :expected '~form, :actual nil})
          (catch Exception e#
            (let [m# (ex-data e#)]
              (if (= ~key (:type m#))
                (do-report {:type :pass, :message ~msg,
                            :expected '~form, :actual e#})
                (do-report {:type :fail, :message ~msg,
                            :expected '~form, :actual e#})))
            e#))))

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
  {:style/indent 1}
  [s & body]
  `(testing ~s
     (prepare-cavia!)
     ~@body))

(defn disable-log-fixture [f]
  (binding [*logger-factory* disabled-logger-factory]
    (f)))

(def test-ref-seq-file (cavia/resource prof "test.2bit"))
(def test-ref-gene-file (cavia/resource prof "test-refGene.txt.gz"))

(def test-chain-file (io/resource "hg19ToHg38.over.chain.gz"))

(def test-gtf-file "./test-resources/gtf_parse_test.gtf")

(def test-gff3-file "./test-resources/gff3_parse_test.gff3")
