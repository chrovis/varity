(ns varity.fusion-test
  (:require [clojure.test :refer [deftest is are testing]]
            [cljam.io.sequence :as cseq]
            [varity.ref-gene :as rg]
            [varity.fusion :as fusion]
            [varity.t-common :as t-common
             :refer [defslowtest cavia-testing]]))

(def eml4-alk
  {:chr "chr2",
   :pos 29223826,
   :ref "G",
   :alt ["G]chr2:42299120]"],
   :info {:SVTYPE "INV"},
   :id "EML4-ALK",
   :pattern ["EML4" :intron 13 nil "ALK" :intron 19 nil]
   :gene1 "NM_019063.5",
   :gene2 "NM_004304.5",
   :aa (str
        "MDGFAGSLDDSISAASTSDVQDRLSALESRVQQQEDEITVLKAALADVLR"
        "RLAISEDHVASVKKSVSSKGQPSPRAVIPMSCITNGSGANRKPSHTSAVS"
        "IAGKETLSSAAKSGTEKKKEKPQGQREKKEESHSNDQSPQIRASPSPQPS"
        "SQPLQIHRQTPESKNATPTKSIKRPSPAEKSHNSWENSDDSRNKLSKIPS"
        "TPKLIPKVTKTADKHKDVIINQEGEYIKMFMRGRPITMFIPSDVDNYDDI"
        "RTELPPEKLKLEWAYGYRGKDCRANVYLLPTGKIVYFIASVVVLFNYEER"
        "TQRHYLGHTDCVKCLAIHPDKIRIATGQIAGVDKDGRPLQPHVRVWDSVT"
        "LSTLQIIGLGTFERGVGCLDFSKADSGVHLCIIDDSNEHMLTVWDWQKKA"
        "KGAEIKTTNEVVLAVEFHPTDANTIITCGKSHIFFWTWSGNSLTRKQGIF"
        "GKYEKPKFVQCLAFLGNGDVLTGDSGGVMLIWSKTTVEPTPGKGPKVYRR"
        "KHQELQAMQMELQSPEYKLSKLRTSTIMTDYNPNYCFAGKTSSISDLKEV"
        "PRKNITLIRGLGHGAFGEVYEGQVSGMPNDPSPLQVAVKTLPEVCSEQDE"
        "LDFLMEALIISKFNHQNIVRCIGVSLQSLPRFILLELMAGGDLKSFLRET"
        "RPRPSQPSSLAMLDLLHVARDIACGCQYLEENHFIHRDIAARNCLLTCPG"
        "PGRVAKIGDFGMARDIYRASYYRKGGCAMLPVKWMPPEAFMEGIFTSKTD"
        "TWSFGVLLWEIFSLGYMPYPSKSNQEVLEFVTSGGRMDPPKNCPGPVYRI"
        "MTQCWQHQPEDRPNFAIILERIEYCTQDPDVINTALPIEYGPLVEEEEKV"
        "PVRPKDPEGVPPLLVSQQAKREEERSPAAPPPLPTTSSGKAAKKPTAAEI"
        "SVRVPRGPAVEGGHVNMAFSQSNPPSELHKVHGSRNKPTSLWNPTYGSWF"
        "TEKPTKKNNPIAKKEPHDRGNLGLEGSCTVPPNVATGRLPGASLLLEPSS"
        "LTANMKEVPLFRLRHFPCGNVNYGYQQQGLPLEAATAPGAGHYEDTILKS"
        "KNSMNQPGP*")})

(def fus-ddit3
  {:chr "chr12",
   :pos 57520469,
   :ref "C",
   :alt ["C]chr16:31187788]"],
   :info {:SVTYPE "TRA"},
   :id "FUS-DDIT3",
   :pattern ["FUS" :intron 7 nil "DDIT3" :exon 1 :5'-utr]
   :gene1 "NM_004960.4"
   :gene2 "NM_004083.6",
   :aa (str
        "MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGY"
        "GQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGY"
        "GQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYN"
        "PPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQS"
        "GGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGG"
        "RGGMGGSDRGGFNKFGVFKKEVYLHTSPHLKADVLFQTDPTAEMAAESLP"
        "FSFGTLSSWELEAWYEDLQEVLSSDENGGTYVSPPGNEEEESKIFTTLDP"
        "ASLAWLTEEEPEPAEVTSTSQSPHSPDSSQSSLAQEEEEEDQGRTRKRKQ"
        "SGHSPARAGKQRMKEKEQENERKVAQLAEENERLKQEIERLTREVEATRR"
        "ALIDRMVNLHQA*")})

(def kif5b-ret
  {:chr "chr10",
   :pos 32027449,
   :ref "A",
   :alt ["[chr10:43115670[A"],
   :info {:SVTYPE "INV"},
   :id "KIF5B-RET",
   :pattern ["KIF5B" :intron 15 nil "RET" :intron 11 nil]
   :gene1 "NM_004521.3",
   :gene2 "NM_020630.6",
   :aa (str
        "MADLAECNIKVMCRFRPLNESEVNRGDKYIAKFQGEDTVVIASKPYAFDR"
        "VFQSSTSQEQVYNDCAKKIVKDVLEGYNGTIFAYGQTSSGKTHTMEGKLH"
        "DPEGMGIIPRIVQDIFNYIYSMDENLEFHIKVSYFEIYLDKIRDLLDVSK"
        "TNLSVHEDKNRVPYVKGCTERFVCSPDEVMDTIDEGKSNRHVAVTNMNEH"
        "SSRSHSIFLINVKQENTQTEQKLSGKLYLVDLAGSEKVSKTGAEGAVLDE"
        "AKNINKSLSALGNVISALAEGSTYVPYRDSKMTRILQDSLGGNCRTTIVI"
        "CCSPSSYNESETKSTLLFGQRAKTIKNTVCVNVELTAEQWKKKYEKEKEK"
        "NKILRNTIQWLENELNRWRNGETVPIDEQFDKEKANLEAFTVDKDITLTN"
        "DKPATAIGVIGNFTDAERRKCEEEIAKLYKQLDDKDEEINQQSQLVEKLK"
        "TQMLDQEELLASTRRDQDNMQAELNRLQAENDASKEEVKEVLQALEELAV"
        "NYDQKSQEVEDKTKEYELLSDELNQKSATLASIDAELQKLKEMTNHQKKR"
        "AAEMMASLLKDLAEIGIAVGNNDVKEDPKWEFPRKNLVLGKTLGEGEFGK"
        "VVKATAFHLKGRAGYTTVAVKMLKENASPSELRDLLSEFNVLKQVNHPHV"
        "IKLYGACSQDGPLLLIVEYAKYGSLRGFLRESRKVGPGYLGSGGSRNSSS"
        "LDHPDERALTMGDLISFAWQISQGMQYLAEMKLVHRDLAARNILVAEGRK"
        "MKISDFGLSRDVYEEDSYVKRSQGRIPVKWMAIESLFDHIYTTQSDVWSF"
        "GVLLWEIVTLGGNPYPGIPPERLFNLLKTGHRMERPDNCSEEMYRLMLQC"
        "WKQEPDKRPVFADISKDLEKMMVKRRDYLDLAASTPSDSLIYDDGLSEEE"
        "TPLVDCNNAPLPRALPSTWIENKLYGRISHAFTRF*")})

(def ewsr1-fli1
  {:chr "chr11",
   :pos 128799605,
   :ref "C", :alt ["]chr22:29287545]C"],
   :info {:SVTYPE "TRA"},
   :id "EWSR1-FLI1",
   :pattern ["EWSR1" :intron 8 nil "FLI1" :intron 5 nil]
   :gene1 "NM_013986.4",
   :gene2 "NM_002017.5",
   :aa (str
        "MASTDYSTYSQAAAQQGYSAYTAQPTQGYAQTTQAYGQQSYGTYGQPTDV"
        "SYTQAQTTATYGQTAYATSYGQPPTVEGTSTGYTTPTAPQAYSQPVQGYG"
        "TGAYDTTTATVTTTQASYAAQSAYGTQPAYPAYGQQPAATAPTRPQDGNK"
        "PTETSQPQSSTGGYNQPSLGYGQSNYSYPQVPGSYPMQPVTAPPSYPPTS"
        "YSSTQPTSYDQSSYSQQNTYGQPSSYGQQSSYGQQSSYGQQPPTSYPPQT"
        "GSYSQAPSQYSQQSSSYGQQNPSYDSVRRGAWGNNMNSGLNKSPPLGGAQ"
        "TISKNTEQRPQPDPYQILGPTSSRLANPGSGQIQLWQFLLELLSDSANAS"
        "CITWEGTNGEFKMTDPDEVARRWGERKSKPNMNYDKLSRALRYYYDKNIM"
        "TKVHGKRYAYKFDFHGIAQALQPHPTESSMYKYPSDISYMPSYHAHQQKV"
        "NFVPPHPSSMPVTSSSFFGAASQYWTSPTGGIYPNPNVPRHPNTHVPSHL"
        "GSYY*")})

(def ndrg1-erg
  {:chr "chr8",
   :pos 133265103,
   :ref "T",
   :alt ["]chr21:38453057]T"],
   :info {:SVTYPE "TRA"},
   :id "NDRG1-ERG",
   :pattern ["NDRG1" :intron 3 nil "ERG" :intron 3 nil]
   :gene1 "NM_006096.4",
   :gene2 "NM_001243428.1",
   :aa (str
        "MSREMQDVDLAEVKPLVEKGETITGLLQEFDVQEALSVVSEDQSLFECAY"
        "GTPHLAKTEMTASSSSDYGQTSKMSPRVPQQDWLSQPPARVTIKMECNPS"
        "QVNGSRNSPDECSVAKGGKMVGSPDTVGMNYGSYMEEKHMPPPNMTTNER"
        "RVIVPADPTLWSTDHVRQWLEWAVKEYGLPDVNILLFQNIDGKELCKMTK"
        "DDFQRLTPSYNADILLSHLHYLRETPLPHLTSDDVDKALQNSPRLMHARN"
        "TGGAAFIFPNTSVYPEATQRITTRPDLPYEPPRRSAWTGHGHPTPQSKAA"
        "QPSPSTVPKTEDQRPQLDPYQILGPTSSRLANPGSGQIQLWQFLLELLSD"
        "SSNSSCITWEGTNGEFKMTDPDEVARRWGERKSKPNMNYDKLSRALRYYY"
        "DKNIMTKVHGKRYAYKFDFHGIAQALQPHPPESSLYKYPSDLPYMGSYHA"
        "HPQKMNFVAPHPPALPVTSSSFFAAPNPYWNSPTGGIYPNTRLPTSHMPS"
        "HLGTYY*")})

(def rg (->> "test-resources/test-refseq-fusion.txt.gz"
             rg/load-ref-seqs
             rg/index))

(deftest aligned-genes-test
  (are [?r1 ?s1 ?r2 ?s2 ?result-id ?inserted-seq]
       (= [?result-id ?inserted-seq]
          ((juxt (comp :id first) last)
           (#'fusion/aligned-genes
            {:retained ?r1, :id 1}
            {:retained ?r2, :id 2}
            {:strand (case ?s1 :f :forward :r :reverse)}
            {:strand (case ?s2 :f :forward :r :reverse)}
            "A")))       ;;        this              mate
    :L :f :L :f nil nil
    :L :f :L :r 1   "A"  ;;   o----->|A         <------|
    :L :f :R :f 1   "A"  ;;   o----->|A                |------>
    :L :f :R :r nil nil
    :L :r :L :f 2   "T"  ;;   <------|A         o----->|
    :L :r :L :r nil nil
    :L :r :R :f nil nil
    :L :r :R :r 2   "T"  ;;   <------|A                |<-----o
    :R :f :L :f 2   "A"  ;;         A|------>   o----->|
    :R :f :L :r nil nil
    :R :f :R :f nil nil
    :R :f :R :r 2   "A"  ;;         A|------>          |<-----o
    :R :r :L :f nil nil
    :R :r :L :r 1   "T"  ;;         A|<-----o   <------|
    :R :r :R :f 1   "T"  ;;         A|<-----o          |------>
    :R :r :R :r nil nil))

(deftest exon-intron-seq-test
  (are [?gene ?result]
       (= ?result
          (map #(dissoc % :gene :chr :strand :exon-count :exon-index)
               (#'fusion/exon-intron-seq ?gene)))
    {:strand :forward,
     :exon-ranges [[1 2] [5 7] [9 10]]}
    [{:start 1, :end 2, :index 1, :type :exon}
     {:start 3, :end 4, :index 1, :type :intron}
     {:start 5, :end 7, :index 2, :type :exon}
     {:start 8, :end 8, :index 2, :type :intron}
     {:start 9, :end 10, :index 3, :type :exon}]

    {:strand :reverse,
     :exon-ranges [[1 2] [5 7] [9 10]]}
    [{:start 9, :end 10, :index 1, :type :exon}
     {:start 8, :end 8, :index 1, :type :intron}
     {:start 5, :end 7, :index 2, :type :exon}
     {:start 3, :end 4, :index 2, :type :intron}
     {:start 1, :end 2, :index 3, :type :exon}]))

(deftest transcript-regions-test
  (are [?strand1 ?cds1 ?exon1 ?bp1 ?strand2 ?cds2 ?exon2 ?bp2 ?in-seq ?result]
      (= ?result
         (let [gene1 {:strand ?strand1,
                      :tx-start (ffirst ?exon1), :tx-end (last (last ?exon1)),
                      :cds-start (first ?cds1), :cds-end (last ?cds1),
                      :exon-ranges ?exon1}
               gene2 {:strand ?strand2,
                      :tx-start (ffirst ?exon2), :tx-end (last (last ?exon2)),
                      :cds-start (first ?cds2), :cds-end (last ?cds2),
                      :exon-ranges ?exon2}
               split #(#'fusion/split-at-breakpoint
                       %1 (#'fusion/exon-intron-seq %2))
               [xs1 [x1 & _]] (split ?bp1 gene1)
               [_ [x2 & xs2]] (split ?bp2 gene2)]
           (-> (#'fusion/transcript-regions
                xs1 ?bp1 x1 ?in-seq x2 ?bp2 xs2)
               (update :breakpoint-regions
                       (comp seq (partial map (juxt :start :end :type :utr))))
               (update :transcript-regions
                       (comp seq (partial map (juxt :start :end :strand
                                                    :inserted-seq)))))))
    ;; exon-exon, with an inserted seq => fusion
    :forward [2 9] [[1 2] [5 7] [9 10]] 6
    :reverse [4 8] [[1 5] [8 10] [14 20] [22 28]] 9
    "ATT"
    {:fusion? true,
     :breakpoint-regions [[5 7 :exon nil]
                          [8 10 :exon :5'-utr]],
     :transcript-regions [[2 2 :forward nil]
                          [5 6 :forward nil]
                          [nil nil nil "ATT"]
                          [8 9 :reverse nil]
                          [4 5 :reverse nil]]}
    ;; 5'UTR-exon => no fusion
    :forward [5 9] [[1 2] [5 7] [9 10]] 1
    :reverse [4 8] [[1 5] [8 10] [14 20] [22 28]] 5
    nil
    {:fusion? false,
     :breakpoint-regions [[1 2 :exon :5'-utr]
                          [1 5 :exon nil]],
     :transcript-regions nil}
    ;; 3'UTR-exon => no fusion
    :forward [2 9] [[1 2] [5 7] [9 10]] 10
    :reverse [4 8] [[1 5] [8 10] [14 20] [22 28]] 5
    nil
    {:fusion? false,
     :breakpoint-regions [[9 10 :exon :3'-utr]
                          [1 5 :exon nil]],
     :transcript-regions [[2 2 :forward nil]
                          [5 7 :forward nil]
                          [9 9 :forward nil]]}
    ;; intron in 3'UTR-exon => no fusion
    :forward [2 7] [[1 2] [5 7] [9 10]] 8
    :reverse [4 8] [[1 5] [8 10] [14 20] [22 28]] 5
    nil
    {:fusion? false,
     :breakpoint-regions [[8 8 :intron :3'-utr]
                          [1 5 :exon nil]],
     :transcript-regions [[2 2 :forward nil]
                          [5 7 :forward nil]]}
    ;; exon-3'UTR => no fusion
    :forward [2 9] [[1 2] [5 7] [9 10]] 6
    :reverse [4 8] [[1 5] [8 10] [14 20] [22 28]] 3
    nil
    {:fusion? false,
     :breakpoint-regions [[5 7 :exon nil]
                          [1 5 :exon :3'-utr]],
     :transcript-regions [[2 2 :forward nil]
                          [5 6 :forward nil]]}
    ;; intron-3'UTR => no fusion
    :forward [2 9] [[1 2] [5 7] [9 10]] 8
    :reverse [4 8] [[1 5] [8 10] [14 20] [22 28]] 3
    nil
    {:fusion? false,
     :breakpoint-regions [[8 8 :intron nil]
                          [1 5 :exon :3'-utr]],
     :transcript-regions [[2 2 :forward nil]
                          [5 7 :forward nil]]}
    ;; exon-intron => no fusion
    :forward [2 9] [[1 2] [5 7] [9 10]] 7
    :reverse [4 8] [[1 5] [8 10] [14 20] [22 28]] 6
    nil
    {:fusion? false,
     :breakpoint-regions [[5 7 :exon nil]
                          [6 7 :intron nil]],
     :transcript-regions [[2 2 :forward nil]
                          [5 7 :forward nil]]}
    ;; intron-intron => fusion
    :forward [2 9] [[1 2] [5 7] [9 10]] 8
    :reverse [4 8] [[1 5] [8 10] [14 20] [22 28]] 6
    "ATC"
    {:fusion? true,
     :breakpoint-regions [[8 8 :intron nil]
                          [6 7 :intron nil]],
     :transcript-regions [[2 2 :forward nil]
                          [5 7 :forward nil]
                          [4 5 :reverse nil]]}
    ;; intron-exon => fusion
    :forward [2 9] [[1 2] [5 7] [9 10]] 8
    :reverse [4 8] [[1 5] [8 10] [14 20] [22 28]] 8
    "ATG"
    {:fusion? true,
     :breakpoint-regions [[8 8 :intron nil]
                          [8 10 :exon nil]],
     :transcript-regions [[2 2 :forward nil]
                          [5 7 :forward nil]
                          [4 5 :reverse nil]]}))

(defslowtest fusion-transcripts-test
  (cavia-testing "transcripts of fusion genes"
    (with-open [r (cseq/reader t-common/test-ref-seq-file)]
      (doseq [{:keys [chr pos ref id pattern aa gene1 gene2]
               [alt] :alt} [eml4-alk fus-ddit3 kif5b-ret
                            ewsr1-fli1 ndrg1-erg]
              :let [[{:keys [transcript]
                      [{t1 :type i1 :index u1 :utr}
                       {t2 :type i2 :index u2 :utr}] :breakpoint-regions
                      [g1 g2] :genes}] (fusion/fusion-transcripts
                                        r
                                        rg
                                        (fusion/parse-vcf-breakpoints
                                         chr pos ref alt))]]
        (testing id
          (is (= gene1 (:name g1)))
          (is (= gene2 (:name g2)))
          (is (= pattern [(:name2 g1) t1 i1 u1
                          (:name2 g2) t2 i2 u2]))
          (is (= aa transcript)))))))

(deftest parse-vcf-breakpoints
  ;; These test cases are taken from the VCF specification
  (testing "bnd_W"
    (is (= {:breakpoints [{:chr  "2", :pos 321681, :retained :L}
                          {:chr "17", :pos 198982, :retained :L}],
            :inserted-seq nil}
           (fusion/parse-vcf-breakpoints "2" 321681 "G" "G]17:198982]"))))
  (testing "bnd_V"
    (is (= {:breakpoints [{:chr  "2", :pos 321682, :retained :R}
                          {:chr "13", :pos 123456, :retained :L}],
            :inserted-seq nil}
           (fusion/parse-vcf-breakpoints "2" 321682 "T" "]13:123456]T"))))
  (testing "bnd_U"
    (is (= {:breakpoints [{:chr "13", :pos 123456, :retained :L}
                          {:chr "2", :pos 321682, :retained :R}],
            :inserted-seq nil}
           (fusion/parse-vcf-breakpoints "13" 123456 "C" "C[2:321682["))))
  (testing "bnd_X"
    (is (= {:breakpoints [{:chr "13", :pos 123457, :retained :R}
                          {:chr "17", :pos 198983, :retained :R}],
            :inserted-seq nil}
           (fusion/parse-vcf-breakpoints "13" 123457 "A" "[17:198983[A"))))
  (testing "bnd_Y"
    (is (= {:breakpoints [{:chr "17", :pos 198982, :retained :L}
                          {:chr  "2", :pos 321681, :retained :L}],
            :inserted-seq nil}
           (fusion/parse-vcf-breakpoints "17" 198982 "A" "A]2:321681]"))))
  (testing "bnd_Z"
    (is (= {:breakpoints [{:chr "17", :pos 198983, :retained :R}
                          {:chr "13", :pos 123457, :retained :R}],
            :inserted-seq nil}
           (fusion/parse-vcf-breakpoints "17" 198983 "C" "[13:123457[C"))))
  (testing "bnd_V ins"
    (is (= {:breakpoints [{:chr  "2", :pos 321682, :retained :R}
                          {:chr "13", :pos 123456, :retained :L}],
            :inserted-seq "AGTNNNNNCA"}
           (fusion/parse-vcf-breakpoints "2" 321682 "T" "]13:123456]AGTNNNNNCAT"))))
  (testing "bnd_U ins"
    (is (= {:breakpoints [{:chr "13", :pos 123456, :retained :L}
                          {:chr  "2", :pos 321682, :retained :R}],
            :inserted-seq "AGTNNNNNCA"}
           (fusion/parse-vcf-breakpoints "13" 123456 "C" "CAGTNNNNNCA[2:321682[")))))

(deftest variant->breakpoints-test
  (testing "breakends"
    (is (= [{:breakpoints [{:chr  "2", :pos 321681, :retained :L}
                           {:chr "17", :pos 198982, :retained :L}],
             :inserted-seq nil}]
           (fusion/variant->breakpoints
            {:chr "2", :pos 321681, :ref "G", :alt ["G]17:198982]"],
             :info {:SVTYPE "BND", :MATEID "bnd_Y"}}))))
  (testing "duplication in alt"
    (is (= [{:breakpoints [{:chr "2", :pos 321681, :retained :R}
                           {:chr "2", :pos 321800, :retained :L}],
             :inserted-seq nil}]
           (fusion/variant->breakpoints
            {:chr "2", :pos 321681, :ref "G", :alt ["<DUP>"],
             :info {:SVTYPE "CNV", :END 321800}}))))
  (testing "duplication in info"
    (is (= [{:breakpoints [{:chr "2", :pos 321681, :retained :R}
                           {:chr "2", :pos 321800, :retained :L}],
             :inserted-seq nil}]
           (fusion/variant->breakpoints
            {:chr "2", :pos 321681, :ref "G", :alt [nil],
             :info {:SVTYPE "DUP", :END 321800}}))))
  (testing "deletion in alt"
    (is (= [{:breakpoints [{:chr "2", :pos 321681, :retained :L}
                           {:chr "2", :pos 321800, :retained :R}],
             :inserted-seq nil}]
           (fusion/variant->breakpoints
            {:chr "2", :pos 321681, :ref "G", :alt ["<DEL>"],
             :info {:SVTYPE "CNV", :END 321800}}))))
  (testing "deleteion in info"
    (is (= [{:breakpoints [{:chr "2", :pos 321681, :retained :L}
                           {:chr "2", :pos 321800, :retained :R}],
             :inserted-seq nil}]
           (fusion/variant->breakpoints
            {:chr "2", :pos 321681, :ref "G", :alt [nil],
             :info {:SVTYPE "DEL", :END 321800}}))))
  (testing "nil for CNVs"
    (is (= [nil]
           (fusion/variant->breakpoints
            {:chr "2", :pos 321681, :ref "G", :alt ["<CNV>"],
             :info {:SVTYPE "CNV", :END 321800}}))))
  (testing "not supporting inversion"
    (is (= [nil]
           (fusion/variant->breakpoints
            {:chr "2", :pos 321681, :ref "G", :alt ["<INV>"],
             :info {:SVTYPE "INV", :END 321800}})))))
