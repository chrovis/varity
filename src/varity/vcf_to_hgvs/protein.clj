(ns varity.vcf-to-hgvs.protein
  (:require [clojure.core.memoize :as memo]
            [clojure.pprint :as pp]
            [clojure.set :as s]
            [clojure.string :as string]
            [clojure.tools.logging :as log]
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.core :as hgvs]
            [clj-hgvs.mutation :as mut]
            [cljam.io.sequence :as cseq]
            [cljam.util.sequence :as util-seq]
            [proton.core :as proton]
            [proton.string :as pstring]
            [varity.codon :as codon]
            [varity.ref-gene :as rg]
            [varity.vcf-to-hgvs.common :refer [diff-bases] :as common]))

(defn- overlap-exon-intron-boundary?*
  [exon-ranges pos ref alt]
  (let [nref (count ref)
        nalt (count alt)
        [_ _ d _] (diff-bases ref alt)
        pos (+ pos d)
        nref (- nref d)]
    (boolean
     (and (not (= 1 nref nalt))
          (not= 1 (count exon-ranges))
          (some (fn [[s e]]
                  (and (not= s e)
                       (or (and (< pos s) (<= s (+ pos nref -1)))
                           (and (<= pos e) (< e (+ pos nref -1))))))
                exon-ranges)))))

(def ^:private default-memoize-threshold 1024)

(defn ^:private default-memoize-factory
  [f]
  (memo/lru f :lru/threshold default-memoize-threshold))

(def ^:private overlap-exon-intron-boundary?-impl (atom (default-memoize-factory overlap-exon-intron-boundary?*)))

(defn overlap-exon-intron-boundary?
  "Returns a boolean whether a variant overlaps boundary of exon/intron or not."
  [exon-ranges pos ref alt]
  (@overlap-exon-intron-boundary?-impl exon-ranges pos ref alt))

(defn reset-overlap-exon-intron-boundary?!
  "Resets the memoized function to change memoization strategy."
  [memoize-factory]
  (reset! overlap-exon-intron-boundary?-impl (memoize-factory overlap-exon-intron-boundary?*)))

(defn alt-exon-ranges
  "Returns exon ranges a variant applied."
  [exon-ranges pos ref alt]
  (let [nref (count ref)
        nalt (count alt)
        typ (cond
              (< nref nalt) :ins
              (> nref nalt) :del
              :else :same)
        tpos (+ pos (min nref nalt))
        d (Math/abs (- nref nalt))]
    (if (overlap-exon-intron-boundary? exon-ranges pos ref alt)
      (do (log/warn "Variants overlapping a boundary of exon/intron are unsupported")
          nil)
      (->> exon-ranges
           (keep (fn [[s e]]
                   (case typ
                     :ins (cond
                            (< tpos s) [(+ s d) (+ e d)]
                            (<= s tpos e) [s (+ e d)]
                            :else [s e])
                     :del (let [dels tpos
                                dele (dec (+ tpos d))]
                            (cond
                              (< dele s) [(- s d) (- e d)]
                              (= s e dels dele) [dels dele]
                              (<= dels s) (when (< dele e) [dels (- e d)])
                              (<= dels e) (if (< dele e)
                                            [s (- e d)]
                                            [s (dec dels)])
                              :else [s e]))
                     :same [s e])))
           vec))))

(defn exon-sequence
  "Extracts bases in exon from supplied sequence, returning the sequence of
  bases as string."
  ([sequence* start exon-ranges]
   (exon-sequence sequence* start (+ start (count sequence*)) exon-ranges))
  ([sequence* start end exon-ranges]
   (let [limit (+ (count sequence*) start -1)]
     (->> exon-ranges
          (keep (fn [[s e]]
                  (cond
                    (< limit s) nil
                    (< e start) nil
                    (< end s) nil
                    (and (<= start s) (<= e end)) [s (min e limit)]
                    (and (< s start) (< end e)) [start (min end limit)]
                    (<= s start) [start (min e limit)]
                    (<= end e) [s (min end limit)])))
          (map (fn [[s e]]
                 (subs sequence* (- s start) (inc (- e start)))))
          (apply str)))))

(defn- is-insertion-variant?
  [ref alt]
  (let [[del ins offset _] (diff-bases ref alt)
        ndel (count del)
        nins (count ins)]
    (and (= ndel 0) (<= 1 nins) (= offset 1))))

(defn- cds-start-upstream-to-cds-variant?
  [cds-start pos ref]
  (and (< pos cds-start)
       (<= cds-start (dec (+ pos (count ref))))))

(defn- cds-to-cds-end-downstream-variant?
  [cds-end pos ref]
  (and (<= pos cds-end)
       (< cds-end (dec (+ pos (count ref))))))

(defn- make-alt-up-exon-seq
  [alt-up-exon-seq start end alt-exon-ranges strand]
  (let [alt-up-exon-seq* (exon-sequence alt-up-exon-seq start end alt-exon-ranges)]
    (if (= strand :forward)
      (subs alt-up-exon-seq* (mod (count alt-up-exon-seq*) 3))
      alt-up-exon-seq*)))

(defn- make-alt-down-exon-seq
  [alt-down-exon-seq start end alt-exon-ranges strand]
  (let [alt-down-exon-seq* (exon-sequence alt-down-exon-seq start end alt-exon-ranges)
        nalt-down-exon-seq (count alt-down-exon-seq*)]
    (if (= strand :forward)
      alt-down-exon-seq*
      (subs alt-down-exon-seq* 0 (- nalt-down-exon-seq (mod nalt-down-exon-seq 3))))))

(defn- make-ter-site-adjusted-alt-seq
  [alt-exon-seq alt-up-exon-seq alt-down-exon-seq strand cds-start cds-end pos ref ref-include-ter-site]
  (cond
    (and (= strand :forward)
         (or (cds-to-cds-end-downstream-variant? cds-end pos ref)
             ref-include-ter-site))
    (str alt-exon-seq alt-down-exon-seq)

    (and (= strand :reverse)
         (or (cds-start-upstream-to-cds-variant? cds-start pos ref)
             ref-include-ter-site))
    (str alt-up-exon-seq alt-exon-seq)

    :else
    alt-exon-seq))

(defn- include-utr-ini-site-boundary?
  [{:keys [strand cds-start cds-end]} pos ref alt]
  (let [ini-site-boundary (set (cond
                                 (= strand :forward) (range (- cds-start 1) (inc cds-start))
                                 (= strand :reverse) (range cds-end (inc (+ cds-end 1)))))
        nref (count ref)
        alt-positions (set (if (= (first ref) (first alt))
                             (range (inc pos) (+ (inc pos) (dec nref)))
                             (range pos (+ pos nref))))]
    (= 2 (count (s/intersection ini-site-boundary alt-positions)))))

(defn- ini-site-affected?
  [ref-cds-seq alt-cds-seq]
  (not= (take 3 ref-cds-seq) (take 3 alt-cds-seq)))

(defn- include-ter-site?
  [{:keys [strand cds-start cds-end]} pos ref alt]
  (let [ter-site-positions (set (cond
                                  (= strand :forward) (range (- cds-end 2) (inc cds-end))
                                  (= strand :reverse) (range cds-start (inc (+ cds-start 2)))))
        [del _ offset _] (diff-bases ref alt)
        pos (+ pos offset)
        ndel (count del)
        end-pos (+ pos (if (is-insertion-variant? ref alt) 0 (dec ndel)))
        alt-positions (set (range pos (inc end-pos)))]
    (boolean (seq (s/intersection ter-site-positions alt-positions)))))

(defn- ref-include-from-ter-start-and-over-ter-end?
  [{:keys [strand cds-start cds-end]} pos ref alt]
  (let [[del _ offset _] (diff-bases ref alt)
        pos (+ pos offset)
        ndel (count del)
        ter-start-pos (if (= strand :forward)
                        (- cds-end 2)
                        (+ cds-start 2))
        ter-end-pos (if (= strand :forward)
                      cds-end
                      cds-start)
        pos-end (+ pos (if (= ndel 0) 0 (dec ndel)))]
    (if (= strand :forward)
      (and (= pos ter-start-pos) (<= ter-end-pos pos-end))
      (and (= pos-end ter-start-pos) (<= pos ter-end-pos)))))

(defn- ref-include-from-ter-upstream-and-over-ter-end?
  [{:keys [strand cds-start cds-end]} pos ref alt]
  (let [[del _ offset _] (diff-bases ref alt)
        pos (+ pos offset)
        ndel (count del)
        ter-start-pos (if (= strand :forward)
                        (- cds-end 2)
                        (+ cds-start 2))
        ter-end-pos (if (= strand :forward)
                      cds-end
                      cds-start)
        pos-end (+ pos (if (= ndel 0) 0 (dec ndel)))]
    (if (= strand :forward)
      (and (< pos ter-start-pos) (< ter-end-pos pos-end))
      (and (< ter-start-pos pos-end) (< pos ter-end-pos)))))

(defn- ter-site-same-pos?
  [ref-prot-seq alt-prot-seq]
  (and (string/includes? ref-prot-seq "*")
       (string/includes? alt-prot-seq "*")
       (let [ref-ter-pos (count ref-prot-seq)
             alt-ter-pos (inc (count (first (string/split alt-prot-seq #"\*"))))]
         (= ref-ter-pos alt-ter-pos))))

(defn- cds-start-upstream?
  [cds-start pos ref alt]
  (let [[del _ offset _] (diff-bases ref alt)
        ndel (count del)
        pos (+ pos offset)
        pos-end (+ pos (if (= ndel 0) 0 (dec ndel)))]
    (if (is-insertion-variant? ref alt)
      (<= pos cds-start)
      (<= pos pos-end (dec cds-start)))))

(defn- cds-end-downstream?
  [cds-end pos ref alt]
  (let [[del _ offset _] (diff-bases ref alt)
        ndel (count del)
        pos (+ pos offset)
        pos-end (+ pos (if (= ndel 0) 0 (dec ndel)))]
    (if (is-insertion-variant? ref alt)
      (< cds-end pos)
      (<= (inc cds-end) pos pos-end))))

(defn- cds-variant?
  [cds-start cds-end pos ref alt]
  (let [[del _ offset _] (diff-bases ref alt)
        ndel (count del)
        pos (+ pos offset)
        pos-end (+ pos (if (= ndel 0) 0 (dec ndel)))]
    (if (is-insertion-variant? ref alt)
      (and (< cds-start pos) (<= pos cds-end))
      (<= cds-start pos pos-end cds-end))))

(defn- utr-variant?
  [cds-start cds-end pos ref alt]
  (or (cds-start-upstream? cds-start pos ref alt)
      (cds-end-downstream? cds-end pos ref alt)))

(defn- frameshift-within-cds?
  [{:keys [cds-start cds-end]} pos ref alt]
  (let [[del ins _ _] (diff-bases ref alt)
        ndel (count del)
        nins (count ins)]
    (and (cds-variant? cds-start cds-end pos ref alt)
         (not= 0 (rem (Math/abs (- ndel nins)) 3)))))

(defn- get-alt-cds-start-pos
  [cds-start pos-start pos-end exon-ranges pos*]
  (let [[_ exon-end] (first (filter (fn [[s e]] (<= s pos* e)) exon-ranges))
        pos-end* (inc pos-end)]
    (cond
      (= cds-start pos-start pos-end)
      pos*
      (< exon-end pos-end*)
      (let [[exon-start* _] (first (filter (fn [[s e]] (or (<= pos-end* s)
                                                           (<= s pos-end* e))) exon-ranges))]
        (if (< pos-end* exon-start*)
          exon-start*
          pos-end*))
      :else
      pos-end*)))

(defn- get-alt-cds-end-pos
  [cds-end pos-start pos-end exon-ranges pos*]
  (let [[exon-start _] (first (filter (fn [[s e]] (<= s pos* e)) exon-ranges))
        pos-start* (dec pos-start)]
    (cond
      (= cds-end pos-start pos-end)
      pos*
      (< pos-start* exon-start)
      (let [[_ exon-end*] (first (filter (fn [[s e]] (or (<= e pos-start*)
                                                         (<= s pos-start* e))) (reverse exon-ranges)))]
        (if (< exon-end* pos-start*)
          exon-end*
          pos-start*))
      :else
      pos-start*)))

(defn- in-frame?
  [pos ref alt {:keys [cds-start cds-end] :as _rg}]
  (let [[del ins offset] (diff-bases ref alt)
        start-pos* (+ pos offset)
        end-pos* (+ pos offset (dec (count del)))
        start-pos (if (< start-pos* cds-start)
                    cds-start
                    start-pos*)
        end-pos (if (< cds-end end-pos*)
                  cds-end
                  end-pos*)
        ndel (inc (- end-pos start-pos))
        nins (count ins)
        multiple-of-3? #(= 0 (rem % 3))]
    (if (or (< start-pos* end-pos* cds-start)
            (< cds-end start-pos* end-pos*))
      true
      (multiple-of-3? (- ndel nins)))))

(defn- apply-offset
  [pos ref alt cds-start cds-end exon-ranges pos*]
  (let [[del _ offset _] (diff-bases ref alt)
        ndel (count del)
        pos-start (+ pos offset)
        pos-end (+ pos-start (if (zero? ndel) 0 (dec ndel)))
        apply-offset* (fn [exon-ranges*]
                        (ffirst (alt-exon-ranges exon-ranges* pos ref alt)))
        pos** (cond (and (= pos* cds-start) (<= pos-start cds-start pos-end))
                    (get-alt-cds-start-pos cds-start pos-start pos-end exon-ranges pos*)
                    (and (= pos* cds-end) (<= pos-start cds-end pos-end))
                    (get-alt-cds-end-pos cds-end pos-start pos-end exon-ranges pos*)
                    :else
                    pos*)]
    (apply-offset* [[pos** pos**]])))

(defn- read-sequence-info
  [seq-rdr rg pos ref alt]
  (if-let [overlap-exon-intron-boundary (overlap-exon-intron-boundary? (:exon-ranges rg) pos ref alt)]
    {:overlap-exon-intron-boundary overlap-exon-intron-boundary}
    (let [{:keys [chr tx-start tx-end cds-start cds-end exon-ranges strand]} rg
          ref-seq (cseq/read-sequence seq-rdr {:chr chr, :start tx-start, :end tx-end})
          ref-cds-seq (cseq/read-sequence seq-rdr {:chr chr, :start cds-start, :end cds-end})
          ref-cds-exon-seq (exon-sequence ref-cds-seq cds-start cds-end exon-ranges)
          ref-include-utr-ini-site-boundary (include-utr-ini-site-boundary? rg pos ref alt)
          ref-include-ter-site (include-ter-site? rg pos ref alt)
          ref-include-from-ter-start-and-over-ter-end (ref-include-from-ter-start-and-over-ter-end? rg pos ref alt)
          ref-include-from-ter-upstream-and-over-ter-end (ref-include-from-ter-upstream-and-over-ter-end? rg pos ref alt)
          frameshift-within-cds (frameshift-within-cds? rg pos ref alt)
          alt-seq (common/alt-sequence ref-seq tx-start pos ref alt)
          alt-exon-ranges* (alt-exon-ranges exon-ranges pos ref alt)
          apply-offset* (partial apply-offset pos ref alt cds-start cds-end exon-ranges)
          alt-cds-start (apply-offset* cds-start)
          alt-cds-end (apply-offset* cds-end)
          alt-tx-end (apply-offset* tx-end)
          [alt-up-exon-seq alt-cds-exon-seq alt-down-exon-seq] (when (and alt-cds-start alt-cds-end)
                                                                 (pstring/split-at alt-seq [(- alt-cds-start tx-start) (inc (- alt-cds-end tx-start))]))
          alt-cds-exon-seq (exon-sequence alt-cds-exon-seq alt-cds-start alt-cds-end alt-exon-ranges*)
          alt-up-exon-seq (make-alt-up-exon-seq alt-up-exon-seq tx-start (dec alt-cds-start) alt-exon-ranges* strand)
          alt-down-exon-seq (make-alt-down-exon-seq alt-down-exon-seq (inc alt-cds-end) alt-tx-end alt-exon-ranges* strand)
          ter-site-adjusted-alt-seq (make-ter-site-adjusted-alt-seq alt-cds-exon-seq alt-up-exon-seq alt-down-exon-seq
                                                                    strand cds-start cds-end pos ref ref-include-ter-site)
          in-frame (in-frame? pos ref alt rg)
          ref-cds-exon-seq* (cond-> ref-cds-exon-seq
                              (= strand :reverse) util-seq/revcomp)
          alt-cds-exon-seq* (cond-> alt-cds-exon-seq
                              (= strand :reverse) util-seq/revcomp)]
      {:ref-exon-seq ref-cds-exon-seq*
       :ref-prot-seq (codon/amino-acid-sequence ref-cds-exon-seq*)
       :alt-exon-seq alt-cds-exon-seq*
       :alt-prot-seq (codon/amino-acid-sequence alt-cds-exon-seq*)
       :alt-tx-prot-seq (codon/amino-acid-sequence
                         (cond-> (str alt-up-exon-seq alt-cds-exon-seq alt-down-exon-seq)
                           (= strand :reverse) util-seq/revcomp))
       :ini-offset (quot (count (case strand
                                  :forward alt-up-exon-seq
                                  :reverse alt-down-exon-seq))
                         3)
       :c-ter-adjusted-alt-prot-seq (codon/amino-acid-sequence
                                     (cond-> ter-site-adjusted-alt-seq
                                       (= strand :reverse) util-seq/revcomp))
       :alt-rg (when alt-exon-ranges*
                 (assoc rg
                        :exon-ranges alt-exon-ranges*
                        :cds-start alt-cds-start
                        :cds-end alt-cds-end
                        :tx-end alt-tx-end))
       :ref-include-utr-ini-site-boundary ref-include-utr-ini-site-boundary
       :ref-include-ter-site ref-include-ter-site
       :ref-include-from-ter-start-and-over-ter-end ref-include-from-ter-start-and-over-ter-end
       :ref-include-from-ter-upstream-and-over-ter-end ref-include-from-ter-upstream-and-over-ter-end
       :frameshift-within-cds frameshift-within-cds
       :utr-variant (utr-variant? cds-start cds-end pos ref alt)
       :in-frame in-frame
       :ini-site-affected (ini-site-affected? ref-cds-exon-seq* alt-cds-exon-seq*)})))

(defn- protein-position
  "Converts genomic position to protein position. If pos is outside of CDS,
  returns protein position of CDS start or CDS end. If pos is in intron, returns
  the nearest next protein position."
  [pos rg]
  (let [->protein-coord #(inc (quot (dec %) 3))
        pos (proton/clip pos (:cds-start rg) (:cds-end rg))
        pos (->> (:exon-ranges rg)
                 (some (fn [[s e]]
                         (cond
                           (<= s pos e) pos
                           (< pos s) s))))
        cds-coord (rg/cds-coord pos rg)]
    (->protein-coord (:position cds-coord))))

(defn format-alt-prot-seq
  [{:keys [alt-prot-seq alt-tx-prot-seq ini-offset]}]
  (if (= \* (last alt-prot-seq))
    alt-prot-seq
    (let [s (subs alt-tx-prot-seq
                  ini-offset)
          [s-head s-tail] (pstring/split-at s (count alt-prot-seq))]
      (if-let [end (string/index-of s-tail \*)]
        (str s-head
             (subs s-tail 0 (inc end)))
        s))))

(defn- repeat-info*
  [seq* pos ref-only alt-only]
  (when-let [[alt type] (cond
                          (and (empty? ref-only) (seq alt-only)) [alt-only :ins]
                          (and (seq ref-only) (empty? alt-only)) [ref-only :del])]
    (common/repeat-info seq* pos alt type)))

(defn- get-first-diff-aa-info
  [pos ref-seq alt-seq]
  (let [alt-seq* (subs alt-seq 0 (min (count ref-seq) (count alt-seq)))
        [ref-only alt-only offset _] (diff-bases ref-seq alt-seq* (dec pos))]
    (when (and (seq ref-only)
               (seq alt-only))
      {:ppos (+ pos offset)
       :pref (str (first ref-only))
       :palt (str (first alt-only))})))

(defn- first-diff-aa-is-ter-site?
  [pos ref-seq alt-seq]
  (= "*" (:pref (get-first-diff-aa-info pos ref-seq alt-seq))))

(defn- ->protein-variant
  [{:keys [strand] :as rg} pos ref alt
   {:keys [ref-exon-seq ref-prot-seq alt-exon-seq alt-rg ref-include-ter-site
           ref-include-from-ter-start-and-over-ter-end utr-variant ini-site-affected] :as seq-info}
   {:keys [prefer-deletion? prefer-insertion? prefer-deletion-insertion?
           prefer-extension-for-initial-codon-alt?]}]
  (cond
    (:overlap-exon-intron-boundary seq-info)
    {:type :overlap-exon-intron-boundary, :pos nil, :ref nil, :alt nil}

    (or (= ref-exon-seq alt-exon-seq) utr-variant)
    {:type :no-effect, :pos 1, :ref nil, :alt nil}

    (pos? (mod (count ref-exon-seq) 3))
    (do (log/warnf "CDS length is indivisible by 3: %d (%s, %s)"
                   (count ref-exon-seq) (:name rg) (:name2 rg))
        {:type :unknown, :pos nil, :ref nil, :alt nil})

    (not= (last ref-prot-seq) \*)
    (do (log/warnf "Last codon is not stop codon: %s (%s, %s)"
                   (str (last ref-prot-seq)) (:name rg) (:name2 rg))
        {:type :unknown, :pos nil, :ref nil, :alt nil})

    (and (not prefer-extension-for-initial-codon-alt?)
         ini-site-affected)
    {:type :unknown, :pos nil, :ref nil, :alt nil}

    :else
    (let [alt-prot-seq* (format-alt-prot-seq seq-info)
          ppos (protein-position pos rg)
          base-ppos (case strand
                      :forward ppos
                      :reverse (protein-position (+ pos (count ref) -1) rg))
          [_ pref ref-prot-rest] (pstring/split-at
                                  ref-prot-seq
                                  [(dec base-ppos)
                                   (case strand
                                     :forward (protein-position (+ pos (count ref) -1) rg)
                                     :reverse ppos)])
          [_ palt alt-prot-rest] (pstring/split-at
                                  alt-prot-seq*
                                  [(min (dec base-ppos) (count alt-prot-seq*))
                                   (min (case strand
                                          :forward (protein-position (+ pos (count alt) -1) alt-rg)
                                          :reverse (protein-position pos alt-rg))
                                        (count alt-prot-seq*))])
          [pref-only palt-only offset _] (diff-bases pref palt)
          npref (count pref)
          nprefo (count pref-only)
          npalto (count palt-only)
          [unit ref-repeat alt-repeat] (repeat-info* ref-prot-seq
                                                     (+ base-ppos offset)
                                                     pref-only
                                                     palt-only)
          first-diff-aa-is-ter-site (first-diff-aa-is-ter-site? base-ppos
                                                                ref-prot-seq
                                                                alt-prot-seq*)
          t (cond
              (and base-ppos (= (first pref) (first palt) \*)) :substitution
              (and ref-include-from-ter-start-and-over-ter-end
                   (not first-diff-aa-is-ter-site)) :frame-shift
              (= (+ base-ppos offset) (count ref-prot-seq)) (if (ter-site-same-pos? ref-prot-seq alt-prot-seq*)
                                                              :substitution
                                                              :extension)
              (= (+ base-ppos offset) 1) (if (or (= ref-prot-rest alt-prot-rest)
                                                 (and prefer-extension-for-initial-codon-alt?
                                                      (not= (first ref-prot-seq) (first alt-prot-seq*))))
                                           :extension
                                           :frame-shift)
              (and (pos? npref) (= (first palt-only) \*)) :substitution
              (not= ref-prot-rest alt-prot-rest) (cond
                                                   (or (and (= (first alt-prot-rest) \*)
                                                            (>= nprefo npalto)
                                                            (= palt (subs pref 0 (count palt))))
                                                       (= (first palt-only) \*)) :fs-ter-substitution
                                                   ref-include-ter-site :indel
                                                   first-diff-aa-is-ter-site :extension
                                                   :else :frame-shift)
              (or (and (zero? nprefo) (zero? npalto))
                  (and (= nprefo 1) (= npalto 1))) :substitution
              (and ref-include-ter-site (pos? (count pref)) (= (first palt) \*)) :ter-substitution
              (and prefer-deletion? (pos? nprefo) (zero? npalto)) :deletion
              (and prefer-insertion? (zero? nprefo) (pos? npalto)) (if first-diff-aa-is-ter-site
                                                                     :extension
                                                                     :insertion)
              (and (some? unit) (= ref-repeat 1) (= alt-repeat 2)) :duplication
              (and (some? unit) (pos? alt-repeat)) (if (and prefer-deletion-insertion? (pos? npref) (pos? npalto))
                                                     :indel
                                                     :repeated-seqs)
              (and (pos? nprefo) (zero? npalto)) :deletion
              (and (pos? nprefo) (pos? npalto)) (if (= base-ppos 1)
                                                  :extension
                                                  :indel)
              (and (zero? nprefo) (pos? npalto)) (if first-diff-aa-is-ter-site
                                                   :extension
                                                   :insertion)
              :else (throw (ex-info "Unsupported variant" {:type ::unsupported-variant})))]
      {:type (if (#{:fs-ter-substitution :ter-substitution} t) :substitution t)
       :pos base-ppos
       :ref (cond
              (= t :fs-ter-substitution)
              (let [pref-len (count pref)
                    palt-len (count palt)
                    palt-ter-len (inc palt-len)]
                (if (<= pref-len palt-ter-len)
                  (str pref (subs ref-prot-rest 0 (inc (- palt-len pref-len))))
                  (subs pref 0 palt-ter-len)))

              (= t :ter-substitution)
              (first (string/split pref #"\*"))

              :else pref)
       :alt (if (#{:fs-ter-substitution :ter-substitution} t)
              (str palt \*)
              palt)})))

(defn- protein-substitution
  [ppos pref palt {:keys [ref-prot-seq alt-prot-seq]}]
  (let [[s-ref s-alt offset _] (if (and ppos (= (first pref) (first palt) \*))
                                 ["*" "*" 0 nil]
                                 (diff-bases pref palt))]
    (cond
      (and (empty? s-ref) (empty? s-alt))
      (mut/protein-substitution (mut/->long-amino-acid (last pref))
                                (coord/protein-coordinate ppos)
                                (mut/->long-amino-acid (last palt)))
      (empty? s-ref)
      (let [{:keys [ppos pref palt]} (get-first-diff-aa-info ppos ref-prot-seq alt-prot-seq)]
        (mut/protein-substitution (mut/->long-amino-acid pref)
                                  (coord/protein-coordinate ppos)
                                  (mut/->long-amino-acid palt)))
      :else
      (mut/protein-substitution (mut/->long-amino-acid (first s-ref))
                                (coord/protein-coordinate (+ ppos offset))
                                (mut/->long-amino-acid (first s-alt))))))

(defn- protein-deletion
  [ppos pref palt]
  (let [[del _ offset _] (diff-bases pref palt)
        ndel (count del)]
    (mut/protein-deletion (mut/->long-amino-acid (first del))
                          (coord/protein-coordinate (+ ppos offset))
                          (when (> ndel 1)
                            (mut/->long-amino-acid (last del)))
                          (when (> ndel 1)
                            (coord/protein-coordinate (+ ppos offset ndel -1))))))

(defn- protein-duplication
  [ppos pref palt]
  (let [[_ ins offset _] (diff-bases pref palt)
        nins (count ins)]
    (mut/protein-duplication (mut/->long-amino-acid (first ins))
                             (coord/protein-coordinate (- (+ ppos offset) nins))
                             (when (> nins 1) (mut/->long-amino-acid (last ins)))
                             (when (> nins 1) (coord/protein-coordinate (dec (+ ppos offset)))))))

(defn- protein-insertion
  [ppos pref palt seq-info]
  (let [[_ ins offset _] (diff-bases pref palt)
        start (dec (+ ppos offset))
        end (+ ppos offset)]
    (mut/protein-insertion (mut/->long-amino-acid (nth (:ref-prot-seq seq-info) (dec start)))
                           (coord/protein-coordinate start)
                           (mut/->long-amino-acid (nth (:ref-prot-seq seq-info) (dec end)))
                           (coord/protein-coordinate end)
                           (->> (seq ins) (map mut/->long-amino-acid)))))

(defn- protein-repeated-seqs
  [ppos pref palt seq-info]
  (let [[pref-only palt-only offset _] (diff-bases pref palt)
        [unit ref-repeat alt-repeat] (repeat-info* (:ref-prot-seq seq-info)
                                                   (+ ppos offset)
                                                   pref-only
                                                   palt-only)
        nunit (count unit)
        start (+ ppos offset (- (* nunit (min ref-repeat alt-repeat))))
        end (dec (+ start nunit))]
    (mut/protein-repeated-seqs (mut/->long-amino-acid (first unit))
                               (coord/protein-coordinate start)
                               (when (< start end) (mut/->long-amino-acid (last unit)))
                               (when (< start end) (coord/protein-coordinate end))
                               alt-repeat)))

(defn- protein-frame-shift
  [ppos {:keys [ref-prot-seq ref-include-utr-ini-site-boundary
                ref-include-from-ter-start-and-over-ter-end] :as seq-info}]
  (let [alt-prot-seq* (format-alt-prot-seq seq-info)]
    (if (or ref-include-utr-ini-site-boundary
            (not= (first ref-prot-seq) (first alt-prot-seq*)))
      (mut/protein-unknown-mutation)
      (let [{:keys [ppos pref palt]} (get-first-diff-aa-info ppos ref-prot-seq alt-prot-seq*)
            [_ _ offset _] (diff-bases (or pref "") (or palt ""))
            ppos (if ref-include-from-ter-start-and-over-ter-end (count ref-prot-seq) ppos)
            alt-prot-seq* (format-alt-prot-seq seq-info)
            ref (nth ref-prot-seq (dec (+ ppos offset)))
            alt (nth alt-prot-seq* (dec (+ ppos offset)))
            ter-site (-> seq-info
                         format-alt-prot-seq
                         (subs (dec (+ ppos offset)))
                         (string/index-of "*"))]
        (if (= alt \*)
          (protein-substitution (+ ppos offset) (str ref) (str alt) seq-info) ; eventually fs-ter-substitution
          (mut/protein-frame-shift (mut/->long-amino-acid ref)
                                   (coord/protein-coordinate (+ ppos offset))
                                   (mut/->long-amino-acid alt)
                                   (if (and ter-site (pos? ter-site))
                                     (coord/protein-coordinate (inc ter-site))
                                     (coord/unknown-coordinate))))))))

(defn- protein-extension
  [ppos pref palt {:keys [ref-prot-seq alt-tx-prot-seq ini-offset prefer-extension-for-initial-codon-alt?] :as seq-info}]
  (let [alt-prot-seq* (format-alt-prot-seq seq-info)]
    (cond
      (and (not= ppos 1)
           (ter-site-same-pos? ref-prot-seq alt-prot-seq*))
      (protein-substitution (count ref-prot-seq) "*" "*" seq-info)

      (and (= ppos 1) (not prefer-extension-for-initial-codon-alt?))
      (mut/protein-unknown-mutation)

      :else
      (let [[_ ins _ _] (diff-bases (or pref "") (or palt ""))
            ini-site ((comp str first) ref-prot-seq)
            first-diff-aa-info (if (= ppos 1)
                                 {:ppos 1
                                  :pref ini-site}
                                 (get-first-diff-aa-info ppos
                                                         ref-prot-seq
                                                         alt-prot-seq*))
            rest-seq (if (= ppos 1)
                       (-> alt-tx-prot-seq
                           (subs 0 ini-offset)
                           reverse
                           (#(apply str %)))
                       (subs alt-prot-seq* (:ppos first-diff-aa-info)))
            alt-aa (mut/->long-amino-acid (if (= ppos 1)
                                            (or (last ins) (first rest-seq))
                                            (:palt first-diff-aa-info)))
            alt-aa-offset (if (and (= ppos 1) (nil? (last ins))) -1 0)
            new-aa-pos (some-> (string/index-of rest-seq (:pref first-diff-aa-info)) inc (+ alt-aa-offset))]
        (if (and (= ppos 1) (= alt-aa "Ter"))
          (mut/protein-unknown-mutation)
          (mut/protein-extension (if (= ppos 1) (mut/->long-amino-acid ini-site) "Ter")
                                 (coord/protein-coordinate (if (= ppos 1) 1 (:ppos first-diff-aa-info)))
                                 alt-aa
                                 (if (= ppos 1) :upstream :downstream)
                                 (if new-aa-pos
                                   (coord/protein-coordinate new-aa-pos)
                                   (coord/unknown-coordinate))))))))

(defn- protein-indel
  [ppos pref palt {:keys [ref-prot-seq c-ter-adjusted-alt-prot-seq
                          ref-include-ter-site frameshift-within-cds
                          ref-include-from-ter-upstream-and-over-ter-end in-frame] :as seq-info}]
  (let [[pref* palt* ppos*] (if ref-include-ter-site
                              (let [{adjusted-ppos :ppos} (get-first-diff-aa-info ppos ref-prot-seq c-ter-adjusted-alt-prot-seq)
                                    ppos (or adjusted-ppos ppos)
                                    get-seq-between-pos-ter-site (fn [seq pos]
                                                                   (-> (subs seq (dec pos))
                                                                       (string/split #"\*")
                                                                       first))
                                    pref (get-seq-between-pos-ter-site ref-prot-seq ppos)
                                    palt (get-seq-between-pos-ter-site c-ter-adjusted-alt-prot-seq ppos)]
                                [pref palt ppos])
                              [pref palt ppos])
        [del ins offset _] (diff-bases (or pref* "") (or palt* ""))
        ndel (count del)
        alt-retain-ter-site? (if ref-include-ter-site
                               (string/includes? (subs c-ter-adjusted-alt-prot-seq (dec ppos*)) "*")
                               true)]
    (cond
      (first-diff-aa-is-ter-site? ppos*
                                  ref-prot-seq
                                  c-ter-adjusted-alt-prot-seq)
      (protein-extension ppos* pref* palt* seq-info)

      frameshift-within-cds
      (protein-frame-shift ppos* seq-info)

      (every? empty? [del ins])
      (mut/protein-no-effect)

      (= 1 (count del) (count ins))
      (protein-substitution ppos* pref* palt* seq-info)

      (empty? del)
      (protein-insertion ppos* pref* palt* seq-info)

      (empty? ins)
      (protein-deletion ppos* pref* palt*)

      (and ref-include-from-ter-upstream-and-over-ter-end
           (not in-frame))
      (protein-frame-shift ppos* seq-info)

      alt-retain-ter-site?
      (mut/protein-indel (mut/->long-amino-acid (first del))
                         (coord/protein-coordinate (+ ppos* offset))
                         (when (> ndel 1)
                           (mut/->long-amino-acid (last del)))
                         (when (> ndel 1)
                           (coord/protein-coordinate (+ ppos* offset ndel -1)))
                         (->> (seq ins) (map mut/->long-amino-acid)))
      :else
      (mut/protein-unknown-mutation))))

(defn- mutation
  [seq-rdr rg pos ref alt options]
  (let [seq-info (read-sequence-info seq-rdr rg pos ref alt)]
    (when-let [pvariant (->protein-variant rg pos ref alt seq-info options)]
      (let [{ppos :pos, pref :ref, palt :alt}
            (if-not (#{:no-effect :unknown :overlap-exon-intron-boundary} (:type pvariant))
              (common/apply-3'-rule pvariant (:ref-prot-seq seq-info))
              pvariant)
            seq-info (merge seq-info options)]
        (case (:type pvariant)
          :substitution (protein-substitution ppos pref palt seq-info)
          :deletion (protein-deletion ppos pref palt)
          :duplication (protein-duplication ppos pref palt)
          :insertion (protein-insertion ppos pref palt seq-info)
          :indel (protein-indel ppos pref palt seq-info)
          :repeated-seqs (protein-repeated-seqs ppos pref palt seq-info)
          :frame-shift (protein-frame-shift ppos seq-info)
          :extension (protein-extension ppos pref palt seq-info)
          :no-effect (mut/protein-no-effect)
          :unknown (mut/protein-unknown-mutation)
          :overlap-exon-intron-boundary nil)))))

(defn ->hgvs
  ([variant seq-rdr rg]
   (->hgvs variant seq-rdr rg {}))
  ([{:keys [pos ref alt]} seq-rdr rg options]
   (when-let [mutation (mutation seq-rdr rg pos ref alt options)]
     (hgvs/hgvs nil :protein mutation))))

(defn- prot-seq-pstring
  [pref-seq palt-seq start end {:keys [ppos pref palt]}]
  (let [[pref-up _ pref-down] (pstring/split-at pref-seq
                                                [(- ppos start)
                                                 (min (count pref-seq)
                                                      (+ (- ppos start) (count pref)))])
        [palt-up _ palt-down] (pstring/split-at palt-seq
                                                [(- ppos start)
                                                 (min (count palt-seq)
                                                      (+ (- ppos start) (count palt)))])
        frame-shift? (if (<= (count pref) (count palt))
                       (or (empty? palt-down)
                           (nil? (string/index-of pref-down palt-down)))
                       (or (empty? pref-down)
                           (nil? (string/index-of palt-down pref-down))))
        palt-down (if-not frame-shift?
                    (subs palt-down 0 (min (count pref-down) (count palt-down)))
                    palt-down)
        nmut (max (count pref) (count palt))
        nmutseq (if-not frame-shift? nmut 0)
        ticks (->> (iterate #(+ % 10) (inc (- start (mod start 10))))
                   (take-while #(<= % end))
                   (remove #(< % start)))
        tick-intervals (conj
                        (->> ticks
                             (map #(+ % (if (< ppos %) (max 0 (- (count palt) (count pref))) 0)))
                             (map #(- % start))
                             (partition 2 1)
                             (mapv #(- (second %) (first %))))
                        0)]
    (string/join
     \newline
     [(apply pp/cl-format nil
             (apply str "~V@T" (repeat (count ticks) "~VA"))
             (- (first ticks) start) (interleave tick-intervals ticks))
      (pp/cl-format nil "~A~VA~A" pref-up nmutseq pref pref-down)
      (pp/cl-format nil "~A~VA~A" palt-up nmutseq palt palt-down)
      (pp/cl-format nil "~V@T~V@{~A~:*~}" (- ppos start) nmut "^")])))

(defn debug-string
  [{:keys [pos ref alt]} seq-rdr rg]
  (let [seq-info (read-sequence-info seq-rdr rg pos ref alt)
        alt-prot-seq* (format-alt-prot-seq seq-info)
        ppos (protein-position pos rg)
        base-ppos (case (:strand rg)
                    :forward ppos
                    :reverse (protein-position (+ pos (count ref) -1) rg))
        pref (subs (:ref-prot-seq seq-info)
                   (dec base-ppos)
                   (case (:strand rg)
                     :forward (protein-position (+ pos (count ref) -1) rg)
                     :reverse ppos))
        palt (subs alt-prot-seq*
                   (dec base-ppos)
                   (case (:strand rg)
                     :forward (protein-position (+ pos (count alt) -1) rg)
                     :reverse (protein-position (- pos (- (count alt) (count ref))) rg)))
        {base-ppos :pos, pref :ref, palt :alt} (common/apply-3'-rule
                                                {:pos base-ppos, :ref pref, :alt palt}
                                                (:ref-prot-seq seq-info))
        start (max 1 (- base-ppos 20))
        end (min (+ base-ppos 20) (count (:ref-prot-seq seq-info)))
        pref-seq (subs (:ref-prot-seq seq-info) (dec start) end)
        palt-seq (subs alt-prot-seq* (dec start) (+ end (max 0 (- (count palt) (count pref)))))
        [ticks refs alts muts] (string/split-lines
                                (prot-seq-pstring pref-seq palt-seq start end
                                                  {:ppos base-ppos, :pref pref, :palt palt}))]
    (string/join \newline [(str "         " ticks)
                           (str "  p.ref: " refs)
                           (str "  p.alt: " alts)
                           (str "         " muts)])))
