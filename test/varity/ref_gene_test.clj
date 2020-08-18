(ns varity.ref-gene-test
  (:require [clojure.test :refer :all]
            [clojure.string :as string]
            [cljam.io.sequence :as cseq]
            [cljam.io.protocols :as p]
            [clj-hgvs.coordinate :as coord]
            [varity.ref-gene :as rg]
            [varity.t-common :refer :all]))

(def parse-ref-gene-line #'varity.ref-gene/parse-ref-gene-line)

(def ^:private ^:const test-ref-gene
  {:bin 592
   :name "NM_001291366"
   :chr "chr1"
   :strand :reverse
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
   :exon-frames [1 1 0 -1]})

(deftest parse-ref-gene-line-test
  (is (= (parse-ref-gene-line "592\tNM_001291366\tchr1\t-\t975198\t982117\t976171\t981029\t4\t975198,976498,978880,982064,\t976269,976624,981047,982117,\t0\tPERM1\tcmpl\tcmpl\t1,1,0,-1,")
         test-ref-gene)))

(defslowtest in-any-exon?-test
  (cavia-testing "in-any-exon? (slow)"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (is (true? (rg/in-any-exon? "chr7" 55191822 rgidx)))
      (is (false? (rg/in-any-exon? "chr7" 55019367 rgidx))))))

(deftest regions
  (testing "tx-region"
    (are [?ref-gene ?output]
        (= ?output (rg/tx-region ?ref-gene))
      test-ref-gene {:chr "chr1", :start 975199, :end 982117, :strand :reverse}
      (merge test-ref-gene {:cds-start 982118 :cds-end 982117}) {:chr "chr1", :start 975199, :end 982117, :strand :reverse}))
  (testing "cds-region"
    (are [?ref-gene ?output]
        (= ?output (rg/cds-region ?ref-gene))
      test-ref-gene {:chr "chr1", :start 976172, :end 981029, :strand :reverse}
      (merge test-ref-gene {:cds-start 982118 :cds-end 982117}) nil))
  (testing "exon-seq"
    (are [?ref-gene ?output]
        (= ?output (rg/exon-seq ?ref-gene))
      test-ref-gene [{:exon-index 1, :exon-count 4, :chr "chr1", :start 982065, :end 982117, :strand :reverse}
                     {:exon-index 2, :exon-count 4, :chr "chr1", :start 978881, :end 981047, :strand :reverse}
                     {:exon-index 3, :exon-count 4, :chr "chr1", :start 976499, :end 976624, :strand :reverse}
                     {:exon-index 4, :exon-count 4, :chr "chr1", :start 975199, :end 976269, :strand :reverse}]))
  (testing "cds-seq"
    (are [?ref-gene ?output]
        (= ?output (rg/cds-seq ?ref-gene))
      test-ref-gene [{:exon-index 2, :exon-count 4, :chr "chr1", :start 978881, :end 981029, :strand :reverse}
                     {:exon-index 3, :exon-count 4, :chr "chr1", :start 976499, :end 976624, :strand :reverse}
                     {:exon-index 4, :exon-count 4, :chr "chr1", :start 976172, :end 976269, :strand :reverse}])))

(defslowtest regions-slow
  (cavia-testing "exon-seq"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [?nm ?output]
          (= (some->> ?output (map (partial zipmap [:exon-index :exon-count :chr :start :end :strand])))
             (rg/exon-seq (first (rg/ref-genes ?nm rgidx))))
        "NM_000024" [[1 1 "chr5" 148826593 148828634 :forward]]
        "NM_000361" [[1 1 "chr20" 23045633 23049664 :reverse]]
        "NM_000015" [[1 2 "chr8" 18391245 18391345 :forward]
                     [2 2 "chr8" 18399998 18401213 :forward]]
        "NM_000025" [[1 2 "chr8" 37965265 37966666 :reverse]
                     [2 2 "chr8" 37962996 37964239 :reverse]]
        "NR_024420" [[1 2 "chr12" 8390270 8390752 :reverse]
                     [2 2 "chr12" 8356964 8357141 :reverse]]
        "NR_037934" [[1 2 "chr4" 25160672 25160753 :forward]
                     [2 2 "chr4" 25197468 25198505 :forward]]
        "NM_020403" [[1 4 "chr13" 67229780 67230336 :reverse]
                     [2 4 "chr13" 67225405 67228575 :reverse]
                     [3 4 "chr13" 66631210 66631411 :reverse]
                     [4 4 "chr13" 66302834 66305028 :reverse]])))
  (cavia-testing "cds-seq"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [?nm ?output]
          (= (some->> ?output (map (partial zipmap [:exon-index :exon-count :chr :start :end :strand])))
             (rg/cds-seq (first (rg/ref-genes ?nm rgidx))))
        "NM_000024" [[1 1 "chr5" 148826832 148828073 :forward]]
        "NM_000361" [[1 1 "chr20" 23047777 23049504 :reverse]]
        "NM_000015" [[2 2 "chr8" 18400004 18400876 :forward]]
        "NM_000025" [[1 2 "chr8" 37965265 37966469 :reverse]
                     [2 2 "chr8" 37964218 37964239 :reverse]]
        "NR_024420" nil
        "NR_037934" nil
        "NM_020403" [[2 4 "chr13" 67225405 67228440 :reverse]
                     [3 4 "chr13" 66631210 66631411 :reverse]
                     [4 4 "chr13" 66304655 66305028 :reverse]]))))

(deftest read-sequences
  (testing "read-transcript-sequence"
    (let [r (reify p/ISequenceReader
              (p/read-sequence [this {:keys [start end]}]
                (subs "AAAATTTTGGGGCCCCAAAATTTTGGGGCCCC" (dec start) end)))]
      (are [?ref-gene ?seq]
          (= ?seq (rg/read-transcript-sequence r ?ref-gene))
        {:strand :forward :exon-ranges [[2 5]]} "AAAT"
        {:strand :reverse :exon-ranges [[2 5]]} "ATTT"
        {:strand :forward :exon-ranges [[2 5] [8 10]]} "AAATTGG"
        {:strand :reverse :exon-ranges [[2 5] [8 10]]} "CCAATTT")))
  (testing "read-coding-sequence"
    (let [r (reify p/ISequenceReader
              (p/read-sequence [this {:keys [start end]}]
                (subs "AAAATTTTGGGGCCCCAAAATTTTGGGGCCCC" (dec start) end)))]
      (are [?ref-gene ?seq]
          (= ?seq (rg/read-coding-sequence r ?ref-gene))
        {:strand :forward :cds-start 3 :cds-end 9 :exon-ranges [[2 5]]} "AAT"
        {:strand :reverse :cds-start 3 :cds-end 9 :exon-ranges [[2 5]]} "ATT"
        {:strand :forward :cds-start 3 :cds-end 9 :exon-ranges [[2 5] [8 10]]} "AATTG"
        {:strand :reverse :cds-start 3 :cds-end 9 :exon-ranges [[2 5] [8 10]]} "CAATT"
        {:strand :forward :cds-start 8 :cds-end 9 :exon-ranges [[2 5] [8 10]]} "TG"
        {:strand :reverse :cds-start 2 :cds-end 5 :exon-ranges [[2 5] [8 10]]} "ATTT"))))

(defslowtest read-sequences-slow
  (cavia-testing "read-transcript-sequence"
                 (with-open [r (cseq/reader test-ref-seq-file)]
                   (are [?ref-gene ?prefix ?suffix]
                       (let [s (rg/read-transcript-sequence r ?ref-gene)]
                         (and (string/starts-with? s ?prefix)
                              (string/ends-with? s ?suffix)))
                     test-ref-gene "ATGCCG" "CTTTTT")))
  (cavia-testing "read-coding-sequence"
                 (with-open [r (cseq/reader test-ref-seq-file)]
                   (are [?ref-gene ?prefix ?suffix]
                       (let [s (rg/read-coding-sequence r ?ref-gene)]
                         (and (string/starts-with? s ?prefix)
                              (string/ends-with? s ?suffix)))
                     test-ref-gene "ATGGAA" "TCCTAG"))))

(defslowtest read-sequences-rg-slow
  (cavia-testing "read-transcript-sequence"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (with-open [r (cseq/reader test-ref-seq-file)]
        (are [?nm ?prefix ?suffix ?length]
            (= {:prefix ?prefix :suffix ?suffix :length ?length}
               (->> rgidx
                    (rg/ref-genes ?nm)
                    first
                    (rg/read-transcript-sequence r)
                    ((juxt #(string/join (take 6 %))
                           #(string/join (take-last 6 %))
                           count))
                    (zipmap [:prefix :suffix :length])))
          "NM_000024" "GCACAT" "ATTGCA" 2042
          "NM_000361" "GGCTGC" "ATCCCA" 4032
          "NM_000015" "TGAGAT" "TTGTGG" 1317
          "NM_000025" "GCTACT" "TTACAA" 2646
          "NR_024420" "CAGGCA" "TATCAA" 661
          "NR_037934" "AGTTAA" "TAATAA" 1120
          "NM_020403" "AGTTCA" "TCCTGA" 6125
          "NM_201282" "CCCCGG" "ATTTGA" 2239
          "NM_004333" "CGCCTC" "TTATAA" 2946
          "NM_000314" "CCTCCC" "TGACTA" 8702
          "NM_000546" "GATGGG" "GGGGTG" 2591
          "NM_019063" "GGGGCG" "TTCTAA" 5561
          "NM_004304" "AGCTGC" "GACTAA" 6265))))
  (cavia-testing "read-coding-sequence"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (with-open [r (cseq/reader test-ref-seq-file)]
        (are [?nm ?prefix ?suffix ?length]
            (= {:prefix ?prefix :suffix ?suffix :length ?length}
               (->> rgidx
                    (rg/ref-genes ?nm)
                    first
                    (rg/read-coding-sequence r)
                    ((juxt #(string/join (take 6 %))
                           #(string/join (take-last 6 %))
                           count))
                    (zipmap [:prefix :suffix :length])))
          "NM_000024" "ATGGGG" "CTGTAA" 1242
          "NM_000361" "ATGCTT" "CTCTGA" 1728
          "NM_000015" "ATGGAC" "ATTTAG" 873
          "NM_000025" "ATGGCT" "TCTTAG" 1227
          "NR_024420" "" "" 0
          "NR_037934" "" "" 0
          "NM_020403" "ATGGAC" "CTCTAA" 3612
          "NM_201282" "ATGCGA" "TCCTAA" 1887
          "NM_004333" "ATGGCG" "CACTGA" 2301
          "NM_000314" "ATGACA" "GTCTGA" 1212
          "NM_000546" "ATGGAG" "GACTGA" 1182
          "NM_019063" "ATGGAC" "TCCTAA" 2946
          "NM_004304" "ATGGGA" "CCCTGA" 4863)))))

(deftest exon-ranges->intron-ranges-test
  (testing "exon-ranges->intron-ranges"
    (are [exon-ranges r] (= r (rg/exon-ranges->intron-ranges exon-ranges))
             [] []
             [[1 10] [20 30]] [[11 19]]
             [[1000 2000] [3000 3500] [5000 6000]] [[2001 2999] [3501 4999]]
             [[100 300]] [])))

(defslowtest seek-gene-region-test
  (cavia-testing "seek-gene-region (slow)"
    (let [rgidx (rg/index (rg/load-ref-genes test-ref-gene-file))]
      (are [c p tn exs] (= exs
                           (->> (seek-gene-region c p rgidx tn)
                                (map :region-types)
                                (mapv (fn [rt]
                                        (mapv #(vector (:region %) (:index %) (:count %))
                                              rt)))))
        "chr4" 54736520 nil [[["exon" 18 21]] [["exon" 18 21]]]
        "chr7" 116771976 "NM_000245" [[["exon" 14 21]]]
        "chrX" 61197987 nil []
        "chr3" 41224090 "NM_001904" [[["intron" 2 14]]]
        "chr5" 12575053 nil [[["UTR-5" 1 1] ["intron" 1 3]]]
        "chr10" 79512600 "NM_001099692" [[["UTR-5" 1 1]]]
        "chr12" 101128642 "NM_001286615" [[["UTR-3" 1 1]] [["UTR-5" 1 1]]]))))

(deftest cds-coord-test
  ;; 1 [2 3 4] 5 6 7 [8 9 10 11] 12 13 14 15
  (testing "strand +"
    (are [p s e r] (= (coord/format
                       (rg/cds-coord p {:strand :forward
                                        :tx-start 2
                                        :tx-end 11
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
      7  2 3  "*2-1"
      1  9 11 "-5"
      13 2 3  "*7"))
  (testing "strand -"
    (are [p s e r] (= (coord/format
                       (rg/cds-coord p {:strand :reverse
                                        :tx-start 2
                                        :tx-end 11
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
      7  2 3  "-2+1"
      5  9 11 "*2-1"
      13 2 3  "-7"
      1  9 11 "*5")))

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
    (are [c s e r] (= (rg/cds-coord->genomic-pos (coord/parse-coding-dna-coordinate c)
                                                 {:strand :forward
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
      "*2-1" 2 3  7
      "-5"   9 11 1
      "*7"   2 3  13))
  (testing "strand -"
    (are [c s e r] (= (rg/cds-coord->genomic-pos (coord/parse-coding-dna-coordinate c)
                                                 {:strand :reverse
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
      "-2+1" 2 3  7
      "*2-1" 9 11 5
      "-7"   2 3  13
      "*5"   9 11 1)))
