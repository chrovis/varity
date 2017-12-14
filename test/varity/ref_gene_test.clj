(ns varity.ref-gene-test
  (:require [clojure.test :refer :all]
            [clj-hgvs.coordinate :as coord]
            [varity.ref-gene :as rg]
            [varity.t-common :refer :all]))

(def parse-ref-gene-line #'varity.ref-gene/parse-ref-gene-line)

(deftest parse-ref-gene-line-test
  (is (= (parse-ref-gene-line "592\tNM_001291366\tchr1\t-\t975198\t982117\t976171\t981029\t4\t975198,976498,978880,982064,\t976269,976624,981047,982117,\t0\tPERM1\tcmpl\tcmpl\t1,1,0,-1,")
         {:bin 592
          :name "NM_001291366"
          :chr "chr1"
          :strand "-"
          :tx-start 975199
          :tx-end 982117
          :cds-start 976172
          :cds-end 981029
          :exon-count 4
          :exon-ranges [[975199 976269] [976499 976624] [978881 981047] [982065 982117]]
          :score 0
          :name2 "PERM1"
          :cds-start-stat :cmpl
          :cds-end-stat :cmpl
          :exon-frames '(1 1 0 -1)})))

(defslowtest in-any-exon?-test
  (cavia-testing "in-any-exon? (slow)"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (is (true? (rg/in-any-exon? "chr7" 55191822 rgidx)))
      (is (false? (rg/in-any-exon? "chr7" 55019367 rgidx))))))

(deftest cds-coord-test
  ;; 1 [2 3 4] 5 6 7 [8 9 10 11] 12 13 14 15
  (testing "strand +"
    (are [p s e r] (= (coord/format
                       (rg/cds-coord p {:strand "+"
                                        :cds-start s
                                        :cds-end e
                                        :exon-ranges [[2 4] [8 11]]}))
                      r)
      2  2 11 "1"
      9  2 11 "5"
      9  3 11 "4"
      9  2 10 "5"
      5  2 11 "3+1"
      6  2 11 "3+2"
      7  2 11 "4-1"
      3  9 11 "-3"
      9  2 3  "*3"
      5  9 11 "-2+1"
      13 2 3  "*5+2"))
  (testing "strand -"
    (are [p s e r] (= (coord/format
                       (rg/cds-coord p {:strand "-"
                                        :cds-start s
                                        :cds-end e
                                        :exon-ranges [[2 4] [8 11]]}))
                      r)
      9  2 11 "3"
      3  2 11 "6"
      9  3 11 "3"
      9  2 10 "2"
      7  2 11 "4+1"
      6  2 11 "4+2"
      5  2 11 "5-1"
      9  2 3  "-3"
      3  9 11 "*3"
      12 2 3  "-5-1"
      5  9 11 "*2-1")))

(defn- cds-coord
  [chr pos rgidx]
  (->> (rg/ref-genes chr pos rgidx)
       (map #(rg/cds-coord pos %))
       (map coord/format)))

(defslowtest cds-coord-slow-test
  (cavia-testing "cds-coord (slow)"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [c p r] (= (cds-coord c p rgidx) r)
        "chr7"  55191822  '("2573")
        "chr19" 1220596   '("613")
        "chr1"  948129    '("1659+2")
        "chr7"  140753336 '("1799")))))

(deftest cds-coord->genomic-pos-test
  ;; 1 [2 3 4] 5 6 7 [8 9 10 11] 12 13 14 15
  (testing "strand +"
    (are [c s e r] (= (rg/cds-coord->genomic-pos (coord/parse-cdna-coordinate c)
                                                 {:strand "+"
                                                  :cds-start s
                                                  :cds-end e
                                                  :exon-ranges [[2 4] [8 11]]})
                      r)
      "1"    2 11 2
      "5"    2 11 9
      "4"    3 11 9
      "5"    2 10 9
      "3+1"  2 11 5
      "3+2"  2 11 6
      "4-1"  2 11 7
      "-3"   9 11 3
      "*3"   2 3  9
      "-2+1" 9 11 5
      "*5+2" 2 3  13))
  (testing "strand -"
    (are [c s e r] (= (rg/cds-coord->genomic-pos (coord/parse-cdna-coordinate c)
                                                 {:strand "-"
                                                  :cds-start s
                                                  :cds-end e
                                                  :exon-ranges [[2 4] [8 11]]})
                      r)
      "3"    2 11 9
      "6"    2 11 3
      "3"    3 11 9
      "2"    2 10 9
      "4+1"  2 11 7
      "4+2"  2 11 6
      "5-1"  2 11 5
      "-3"   2 3  9
      "*3"   9 11 3
      "-5-1" 2 3  12
      "*2-1" 9 11 5)))
