(ns varity.hgvs-test
  (:require [clojure.test :refer [are testing]]
            [clj-hgvs.core :as hgvs]
            [cljam.io.sequence :as cseq]
            [varity.hgvs :as vhgvs]
            [varity.t-common :refer [defslowtest cavia-testing
                                     test-ref-gene-file
                                     test-ncbi-ref-seq-file
                                     test-ref-seq-file]]
            [varity.ref-gene :as rg]))

(defslowtest find-aliases-test
  (cavia-testing "find HGVS aliases"
    (with-open [rdr (cseq/reader test-ref-seq-file)]
      (let [ref-genes (rg/load-ref-genes test-ref-gene-file)
            ref-seqs (rg/load-ref-seqs test-ncbi-ref-seq-file)
            ref-gene-idx (rg/index ref-genes)
            ref-seq-idx (rg/index ref-seqs)]
        (are [idx hgvs-str expected]
            (= (->> (vhgvs/find-aliases (hgvs/parse hgvs-str) rdr idx)
                    (map hgvs/format)
                    set)
               expected)
          ;; cf. rs11571587 (+)
          ref-gene-idx "NM_000059:c.162CAA[1]" #{"NM_000059:c.162CAA[1]"
                                                 "NM_000059:c.165_167del"}
          ref-gene-idx "NM_000059:c.165_167del" #{"NM_000059:c.162CAA[1]"
                                                  "NM_000059:c.165_167del"}

          ;; cf. rs727502907 (-)
          ref-gene-idx "NM_004333:c.-95GCCTCC[3]" #{"NM_004333:c.-95GCCTCC[3]"
                                                    "NM_004333:c.-77_-72del"}
          ref-gene-idx "NM_004333:c.-77_-72del" #{"NM_004333:c.-95GCCTCC[3]"
                                                  "NM_004333:c.-77_-72del"}

          ;; cf. rs2307882 (-)
          ref-gene-idx "NM_144639:c.1510-122AG[3]" #{"NM_144639:c.1510-122AG[3]"
                                                     "NM_144639:c.1510-121_1510-120insAGAG"}
          ref-gene-idx "NM_144639:c.1510-121_1510-120insAGAG" #{"NM_144639:c.1510-122AG[3]"
                                                                "NM_144639:c.1510-121_1510-120insAGAG"}

          ref-seq-idx "NM_000059.4:c.162CAA[1]" #{"NM_000059.4:c.162CAA[1]"
                                                 "NM_000059.4:c.165_167del"}
          ref-seq-idx "NM_000059.4:c.165_167del" #{"NM_000059.4:c.162CAA[1]"
                                                  "NM_000059.4:c.165_167del"}

          ref-seq-idx "NM_004333.6:c.-95GCCTCC[3]" #{"NM_004333.6:c.-95GCCTCC[3]"
                                                    "NM_004333.6:c.-77_-72del"}
          ref-seq-idx "NM_004333.6:c.-77_-72del" #{"NM_004333.6:c.-95GCCTCC[3]"
                                                  "NM_004333.6:c.-77_-72del"}

          ref-seq-idx "NM_144639.3:c.1510-122AG[3]" #{"NM_144639.3:c.1510-122AG[3]"
                                                     "NM_144639.3:c.1510-121_1510-120insAGAG"}
          ref-seq-idx "NM_144639.3:c.1510-121_1510-120insAGAG" #{"NM_144639.3:c.1510-122AG[3]"
                                                                 "NM_144639.3:c.1510-121_1510-120insAGAG"})
        (testing "false coordinate"
          (are [s] (thrown-with-error-type?
                    ::vhgvs/invalid-variant
                    (vhgvs/find-aliases (hgvs/parse s) rdr ref-gene-idx))
            "NM_001276760:c.559+1G>T" ; actually "NM_001276760:c.560G>T"
            "NM_005228:c.2574-1T>G" ; actually "NM_005228:c.2573T>G"
            ))
        (let [find-rg (fn [refs transcript]
                        (first (filter #(= (:name %) transcript) refs)))]
          (are [refs s transcript expected] (= (->> (vhgvs/find-aliases (hgvs/parse s) rdr (find-rg refs transcript))
                                                    (map hgvs/format)
                                                    set)
                                               expected)
            ;; cf. rs11571587 (+)
            ref-genes "NM_000059:c.162CAA[1]" "NM_000059" #{"NM_000059:c.162CAA[1]"
                                                            "NM_000059:c.165_167del"}
            ref-genes "NM_000059:c.165_167del" "NM_000059" #{"NM_000059:c.162CAA[1]"
                                                             "NM_000059:c.165_167del"}

            ;; cf. rs727502907 (-)
            ref-genes "NM_004333:c.-95GCCTCC[3]" "NM_004333" #{"NM_004333:c.-95GCCTCC[3]"
                                                               "NM_004333:c.-77_-72del"}
            ref-genes "NM_004333:c.-77_-72del" "NM_004333" #{"NM_004333:c.-95GCCTCC[3]"
                                                             "NM_004333:c.-77_-72del"}

            ;; cf. rs2307882 (-)
            ref-genes "NM_144639:c.1510-122AG[3]" "NM_144639" #{"NM_144639:c.1510-122AG[3]"
                                                                "NM_144639:c.1510-121_1510-120insAGAG"}
            ref-genes "NM_144639:c.1510-121_1510-120insAGAG" "NM_144639" #{"NM_144639:c.1510-122AG[3]"
                                                                           "NM_144639:c.1510-121_1510-120insAGAG"}

            ref-seqs "NM_000059.4:c.162CAA[1]" "NM_000059.4" #{"NM_000059.4:c.162CAA[1]"
                                                               "NM_000059.4:c.165_167del"}
            ref-seqs "NM_000059.4:c.165_167del" "NM_000059.4" #{"NM_000059.4:c.162CAA[1]"
                                                                "NM_000059.4:c.165_167del"}

            ;; cf. rs727502907 (-)
            ref-seqs "NM_004333.6:c.-95GCCTCC[3]" "NM_004333.6" #{"NM_004333.6:c.-95GCCTCC[3]"
                                                                  "NM_004333.6:c.-77_-72del"}
            ref-seqs "NM_004333.6:c.-77_-72del" "NM_004333.6" #{"NM_004333.6:c.-95GCCTCC[3]"
                                                                "NM_004333.6:c.-77_-72del"}

            ;; cf. rs2307882 (-)
            ref-seqs "NM_144639.3:c.1510-122AG[3]" "NM_144639.3" #{"NM_144639.3:c.1510-122AG[3]"
                                                                   "NM_144639.3:c.1510-121_1510-120insAGAG"}
            ref-seqs "NM_144639.3:c.1510-121_1510-120insAGAG" "NM_144639.3" #{"NM_144639.3:c.1510-122AG[3]"
                                                                              "NM_144639.3:c.1510-121_1510-120insAGAG"}))))))
