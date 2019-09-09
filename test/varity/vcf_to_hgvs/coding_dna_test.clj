(ns varity.vcf-to-hgvs.coding-dna-test
  (:require [clojure.string :as string]
            [clojure.test :refer :all]
            [varity.vcf-to-hgvs.coding-dna :as coding-dna]))

(deftest sequence-pstring-test
  (are [ref-seq alt-seq start end m rg e]
      (= (#'coding-dna/sequence-pstring ref-seq alt-seq start end m rg) e)
    ;; NM_005228:c.2573T>G
    "CAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGGA"
    "CAAGATCACAGATTTTGGGCGGGCCAAACTGCTGGGTGCGGA"
    55191802 55191843 {:pos 55191822, :ref "T", :alt "G"}
    {:strand :forward, :tx-start 55019032, :tx-end 55207338, :cds-start 55019278, :cds-end 55205617,
     :exon-ranges [[55019032 55019365] [55142286 55142437] [55143305 55143488] [55146606 55146740] [55151294 55151362] [55152546 55152664] [55154011 55154152] [55155830 55155946] [55156533 55156659] [55156759 55156832] [55157663 55157753] [55160139 55160338] [55161499 55161631] [55163733 55163823] [55165280 55165437] [55171175 55171213] [55172983 55173124] [55173921 55174043] [55174722 55174820] [55181293 55181478] [55191719 55191874]]}
    (string/join \newline ["         2562      2572      2582      2592"
                           "CAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGGA"
                           "CAAGATCACAGATTTTGGGCGGGCCAAACTGCTGGGTGCGGA"
                           "                    ^"])

    ;; NM_004369:c.6063+6[9]
    "ctggtctgtgtataaagacataaaaaaaactcacAAGCTGCT"
    "ctggtctgtgtataaagacataaaaaaaaactcacAAGCTGCT"
    237363219 237363260 {:pos 237363239, :ref "t", :alt "ta"}
    {:strand :reverse, :tx-start 237324012, :tx-end 237414207, :cds-start 237324774, :cds-end 237396817,
     :exon-ranges [[237363253 237363398] [237364350 237364428] [237365698 237366035] [237366687 237367286] [237368563 237369177] [237371732 237372337] [237374412 237375020] [237376772 237377344] [237378636 237379235] [237380915 237381499] [237387582 237388184] [237394587 237395204] [237396727 237396847] [237413953 237414207]]}
    (string/join \newline ["         6063+2    6063+12  6063+22   6063+32"
                           "AGCAGCTTgtgagtttttttt atgtctttatacacagaccag"
                           "AGCAGCTTgtgagtttttttttatgtctttatacacagaccag"
                           "                    ^^"])))
