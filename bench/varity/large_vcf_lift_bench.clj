(ns varity.large-vcf-lift-bench
  (:require [libra.bench :refer [defbench is]]
            [libra.criterium :as c]
            [cljam.io.sequence :as cseq]
            [cljam.io.vcf :as vcf]
            [varity.chain :as ch]
            [varity.vcf-lift :as vcf-lift]
            [varity.t-common
             :refer [test-ref-seq-file
                     test-large-vcf-file
                     test-large-chain-file
                     prepare-cavia!]]))

(defbench large-vcf-lift-bench
  (prepare-cavia!)
  (with-open [seq-rdr (cseq/reader test-ref-seq-file)
              vcf-rdr (vcf/reader test-large-vcf-file)]
    (let [chidx (ch/index (ch/load-chain test-large-chain-file))]
      (is (c/quick-bench
           (doall (vcf-lift/liftover-variants
                   seq-rdr chidx (vcf/read-variants vcf-rdr))))))))
