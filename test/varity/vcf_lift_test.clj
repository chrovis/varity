(ns varity.vcf-lift-test
  (:require [clojure.test :refer [deftest is testing]]
            [cljam.io.vcf :as vcf]
            [cljam.io.sequence :as io-seq]
            [varity.chain :as ch]
            [varity.vcf-lift :as vcf-lift]))

(def small-reverse-chains
  (ch/index
   [{:header {:t-name "chr1", :t-start 0, :t-end 540, :t-size 540,
              :q-name "chr1", :q-start 0, :q-end 540, :q-size 540,
              :q-strand :reverse}
     :data [{:chr "chr1", :size 540, :dt nil, :dq nil}]}]))

(def small-forward-chains
  (ch/index
   [{:header {:t-name "chr1", :t-start 0, :t-end 540, :t-size 540,
              :q-name "chr1", :q-start 0, :q-end 540, :q-size 540,
              :q-strand :forward}
     :data [{:chr "chr1", :size 540, :dt nil, :dq nil}]}]))

(deftest small-new-interval-test
  (testing "reverse"
    (is (= (@#'vcf-lift/calc-new-interval "chr1" 1 1 small-reverse-chains)
           {:chr "chr1", :start 540, :end 540, :strand :reverse}))
    (is (= (@#'vcf-lift/calc-new-interval "chr1" 1000 1000
                                          small-reverse-chains)
           nil)))
  (testing "forward"
    (is (= (@#'vcf-lift/calc-new-interval "chr1" 1 1 small-forward-chains)
           {:chr "chr1", :start 1, :end 1, :strand :forward}))
    (is (= (@#'vcf-lift/calc-new-interval "chr1" 1000 1000
                                          small-reverse-chains)
           nil))))

(deftest vcf-reverse-lift-test
  (with-open [vcf-success-rdr
              (vcf/reader "./test-resources/vcf-lift/sample-rev.vcf")
              vcf-reject-rdr
              (vcf/reader "./test-resources/vcf-lift/sample-reject.vcf")
              seq-reader
              (io-seq/reader "./test-resources/vcf-lift/sample.fasta")]
    (is (= (vcf-lift/liftover-variants
            seq-reader
            (ch/index
             (ch/load-chain "./test-resources/vcf-lift/sample-rev.chain"))
            (vcf/read-variants vcf-success-rdr))
           {:success
            [{:alt ["TAT" "TT"], :ref "GCA", :pos 51, :filter [:PASS],
              :id "8", :info nil, :qual 100.0, :chr "chr1"}
             {:alt ["CT" "CA"], :ref "C", :pos 52, :filter [:PASS],
              :id "6", :info nil, :qual 100.0, :chr "chr1"}
             {:alt ["TAT" "T"], :ref "A", :pos 53, :filter [:PASS],
              :id "7", :info nil, :qual 100.0, :chr "chr1"}
             {:alt ["A" "ATCT"], :ref "AT", :pos 53, :filter [:PASS],
              :id "5", :info nil, :qual 100.0, :chr "chr1"}
             {:alt ["T"], :ref "TG", :pos 54, :filter [:PASS],
              :id "2", :info {:END 55}, :qual 100.0, :chr "chr1"}
             {:alt ["T"], :ref "TG", :pos 54, :filter [:PASS],
              :id "4", :info {:END 55}, :qual 100.0, :chr "chr1"}
             {:alt ["GA"], :ref "G", :pos 55, :filter [:PASS],
              :id "1", :info {:END 55}, :qual 100.0, :chr "chr1"}
             {:alt ["T"], :ref "G", :pos 55, :filter [:PASS],
              :id "3", :info {:END 55}, :qual 100.0, :chr "chr1"}]
            :failure nil}))
    (is (= (vcf-lift/liftover-variants
            seq-reader
            (ch/index
             (ch/load-chain "./test-resources/vcf-lift/sample-rev.chain"))
            (vcf/read-variants vcf-reject-rdr))
           {:failure [{:alt ["A"], :ref "AT", :pos 5, :filter [:PASS], :id nil,
                       :info {:END 6}, :qual 100.0, :chr "chr1"}]
            :success nil}))))

(deftest vcf-forward-lift-test
  (with-open [vcf-forward-rdr
              (vcf/reader "./test-resources/vcf-lift/sample-forward.vcf")
              seq-reader
              (io-seq/reader "./test-resources/vcf-lift/sample.fasta")]
    (is (= (vcf-lift/liftover-variants
            seq-reader
            (ch/index (ch/load-chain
                       "./test-resources/vcf-lift/sample-forward.chain"))
            (vcf/read-variants vcf-forward-rdr))
           {:success [{:alt ["AT"], :ref "A", :pos 5, :filter [:PASS],
                       :id "1", :info {:END 5}, :qual 100.0, :chr "chr1"}
                      {:alt ["ATA" "G"], :ref "AT", :pos 5, :filter [:PASS],
                       :id "2", :info nil, :qual 100.0, :chr "chr1"}
                      {:alt ["CG"], :ref "C", :pos 8, :filter [:PASS],
                       :id "3", :info {:END 8}, :qual 100.0, :chr "chr1"}]
            :failure nil}))))

(deftest vcf-multi-lift-test
  (with-open [vcf-forward-rdr
              (vcf/reader "./test-resources/vcf-lift/sample-multi.vcf")
              seq-reader
              (io-seq/reader "./test-resources/vcf-lift/sample.fasta")]
    (is (= (vcf-lift/liftover-variants
            seq-reader
            (ch/index (ch/load-chain
                       "./test-resources/vcf-lift/sample-multi.chain"))
            (vcf/read-variants vcf-forward-rdr))
           {:failure nil,
            :success [{:alt ["A"], :ref "ATGCATGC", :pos 9, :filter [:PASS],
                       :id "1", :info {:END 16}, :qual 100.0, :chr "chr1"}
                      {:alt ["A"], :ref "ATGCATGC", :pos 9, :filter [:PASS],
                       :id "3", :info {:END 16}, :qual 100.0, :chr "chr1"}
                      {:alt ["A"], :ref "ATGCATGCATGCATGCATGC", :pos 9,
                       :filter [:PASS], :id "2", :info {:END 28}, :qual 100.0,
                       :chr "chr1"}
                      {:alt ["T"], :ref "TG", :pos 54, :filter [:PASS], :id "4",
                       :info {:END 55}, :qual 100.0, :chr "chr1"}]}))))

(deftest multi-hit-chain
  (let [chain-index (ch/index (ch/load-chain "./test-resources/vcf-lift/sample-multi.chain"))]
    (testing "descending score chain"
      (is (= (#'vcf-lift/calc-new-interval "chr1" 178 179 chain-index)
             {:chr "chr1", :start 58, :end 59, :strand :forward}))
      (is (= (#'vcf-lift/calc-new-interval "chr1" 182 183 chain-index)
             {:chr "chr1", :start 52, :end 53, :strand :forward}))
      (is (= (#'vcf-lift/calc-new-interval "chr1" 191 192 chain-index)
             {:chr "chr1", :start 51, :end 52, :strand :forward})))
    (testing "ascending score chain"
      (is (= (#'vcf-lift/calc-new-interval "chr4" 131 191 chain-index)
             {:chr "chr1", :start 101, :end 160, :strand :forward})))))

(deftest several-gaps-chain
  (let [chain-index (->> "./test-resources/vcf-lift/sample-several-gaps.chain"
                         ch/load-chain
                         ch/index)]
    ;; 12345  67890  12
    ;; ATGCA--ATTGG--CC
    ;; AAG-ACGTTT-GACCC
    ;; 123 456789 01234
    (is (= {:chr "1", :pos 10, :ref "G", :alt ["T"]}
           (vcf-lift/liftover-variant*
            chain-index {:chr "1", :pos 10, :ref "G", :alt ["T"]})))
    (is (= {:chr "1", :pos 13, :ref "C", :alt ["T"]}
           (vcf-lift/liftover-variant*
            chain-index {:chr "1", :pos 11, :ref "C", :alt ["T"]})))
    (is (= {:chr "1", :pos 14, :ref "C", :alt ["T"]}
           (vcf-lift/liftover-variant*
            chain-index {:chr "1", :pos 12, :ref "C", :alt ["T"]})))))
