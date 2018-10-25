(ns varity.vcf-to-hgvs.genome-test
  (:require [clojure.string :as string]
            [clojure.test :refer :all]
            [varity.vcf-to-hgvs.genome :as genome]))

(deftest sequence-pstring-test
  (are [seq* start end m e] (= (#'genome/sequence-pstring seq* start end m) e)
    "CAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGG"
    2 42 {:pos 22, :ref "T", :alt "G"}
    (string/join \newline ["         11        21        31        41"
                           "CAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGG"
                           "CAAGATCACAGATTTTGGGCGGGCCAAACTGCTGGGTGCGG"
                           "                    ^"])

    "CAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGG"
    2 42 {:pos 22, :ref "T", :alt "TCC"}
    (string/join \newline ["         11        21          31        41"
                           "CAAGATCACAGATTTTGGGCT  GGCCAAACTGCTGGGTGCGG"
                           "CAAGATCACAGATTTTGGGCTCCGGCCAAACTGCTGGGTGCGG"
                           "                    ^^^"])

    "CAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGG"
    2 42 {:pos 22, :ref "TGG", :alt "T"}
    (string/join \newline ["         11        21        31        41"
                           "CAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGG"
                           "CAAGATCACAGATTTTGGGCT  CCAAACTGCTGGGTGCGG"
                           "                    ^^^"])

    "TCAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCG"
    1 41 {:pos 21, :ref "C", :alt "A"}
    (string/join \newline ["1         11        21        31        41"
                           "TCAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCG"
                           "TCAAGATCACAGATTTTGGGATGGCCAAACTGCTGGGTGCG"
                           "                    ^"])))
