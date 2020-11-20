(ns varity.hgvs-test
  (:require [clojure.test :refer [are testing]]
            [clj-hgvs.core :as hgvs]
            [cljam.io.sequence :as cseq]
            [varity.hgvs :as vhgvs]
            [varity.t-common :refer [defslowtest cavia-testing test-ref-gene-file
                                     test-ref-seq-file]]
            [varity.ref-gene :as rg]))

(defslowtest find-aliases-test
  (cavia-testing "find HGVS aliases"
    (with-open [rdr (cseq/reader test-ref-seq-file)]
      (let [rgs (rg/load-ref-genes test-ref-gene-file)
            rgidx (rg/index rgs)
            find-rg (fn [transcript]
                      (first (filter #(= (:name %) transcript) rgs)))]
        (are [s xs] (= (->> (vhgvs/find-aliases (hgvs/parse s) rdr rgidx)
                            (map hgvs/format)
                            set)
                       xs)
          ;; cf. rs11571587 (+)
          "NM_000059:c.162CAA[1]" #{"NM_000059:c.162CAA[1]"
                                    "NM_000059:c.165_167del"}
          "NM_000059.3:c.165_167del" #{"NM_000059:c.162CAA[1]"
                                       "NM_000059:c.165_167del"}

          ;; cf. rs727502907 (-)
          "NM_004333.6:c.-95GCCTCC[3]" #{"NM_004333:c.-95GCCTCC[3]"
                                         "NM_004333:c.-77_-72del"}
          "NM_004333:c.-77_-72del" #{"NM_004333:c.-95GCCTCC[3]"
                                     "NM_004333:c.-77_-72del"}

          ;; cf. rs2307882 (-)
          "NM_144639:c.1510-122AG[3]" #{"NM_144639:c.1510-122AG[3]"
                                        "NM_144639:c.1510-121_1510-120insAGAG"}
          "NM_144639:c.1510-121_1510-120insAGAG" #{"NM_144639:c.1510-122AG[3]"
                                                   "NM_144639:c.1510-121_1510-120insAGAG"})

        (are [s rg xs] (= (->> (vhgvs/find-aliases (hgvs/parse s) rdr rg)
                               (map hgvs/format)
                               set)
                          xs)
          ;; cf. rs11571587 (+)
          "NM_000059:c.162CAA[1]" (find-rg "NM_000059") #{"NM_000059:c.162CAA[1]"
                                                          "NM_000059:c.165_167del"}
          "NM_000059.3:c.165_167del" (find-rg "NM_000059") #{"NM_000059:c.162CAA[1]"
                                                             "NM_000059:c.165_167del"}

          ;; cf. rs727502907 (-)
          "NM_004333.6:c.-95GCCTCC[3]" (find-rg "NM_004333") #{"NM_004333:c.-95GCCTCC[3]"
                                                               "NM_004333:c.-77_-72del"}
          "NM_004333:c.-77_-72del" (find-rg "NM_004333") #{"NM_004333:c.-95GCCTCC[3]"
                                                           "NM_004333:c.-77_-72del"}

          ;; cf. rs2307882 (-)
          "NM_144639:c.1510-122AG[3]" (find-rg "NM_144639") #{"NM_144639:c.1510-122AG[3]"
                                                              "NM_144639:c.1510-121_1510-120insAGAG"}
          "NM_144639:c.1510-121_1510-120insAGAG" (find-rg "NM_144639") #{"NM_144639:c.1510-122AG[3]"
                                                                         "NM_144639:c.1510-121_1510-120insAGAG"})

        (testing "false coordinate"
          (are [s] (thrown-with-error-type?
                    ::vhgvs/invalid-variant
                    (vhgvs/find-aliases (hgvs/parse s) rdr rgidx))
            "NM_001276760:c.559+1G>T" ; actually "NM_001276760:c.560G>T"
            "NM_005228:c.2574-1T>G" ; actually "NM_005228:c.2573T>G"
            ))))))
