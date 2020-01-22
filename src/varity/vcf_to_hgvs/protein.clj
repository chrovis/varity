(ns varity.vcf-to-hgvs.protein
  (:require [clojure.pprint :as pp]
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
    (->> exon-ranges
         (keep (fn [[s e]]
                 (case typ
                   :ins (cond
                          (< tpos s) [(+ s d) (+ e d)]
                          (<= s tpos e) [s (+ e d)]
                          :else [s e])
                   :del (let [dels tpos
                              dele (+ tpos d)]
                          (cond
                            (< dele s) [(- s d) (- e d)]
                            (<= dels s) (if (< dele e) [dels (- e d)])
                            (<= dels e) (if (< dele e)
                                          [s (- e d)]
                                          [s (dec dels)])
                            :else [s e]))
                   :same [s e])))
         vec)))

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

(defn- read-exon-sequence
  [seq-rdr chr start end exon-ranges]
  (exon-sequence (cseq/read-sequence seq-rdr {:chr chr, :start start, :end end})
                 start end exon-ranges))

(defn- read-sequence-info
  [seq-rdr rg pos ref alt]
  (let [{:keys [chr tx-start tx-end cds-start cds-end exon-ranges strand]} rg
        ref-seq (cseq/read-sequence seq-rdr {:chr chr, :start cds-start, :end cds-end})
        alt-seq (common/alt-sequence ref-seq cds-start pos ref alt)
        alt-exon-ranges* (alt-exon-ranges exon-ranges pos ref alt)
        ref-exon-seq1 (exon-sequence ref-seq cds-start exon-ranges)
        ref-up-exon-seq1 (read-exon-sequence seq-rdr chr tx-start (dec cds-start) exon-ranges)
        ref-up-exon-seq1 (subs ref-up-exon-seq1 (mod (count ref-up-exon-seq1) 3))
        ref-down-exon-seq1 (read-exon-sequence seq-rdr chr (inc cds-end) tx-end exon-ranges)
        nref-down-exon-seq1 (count ref-down-exon-seq1)
        ref-down-exon-seq1 (subs ref-down-exon-seq1 0 (- nref-down-exon-seq1 (mod nref-down-exon-seq1 3)))
        alt-exon-seq1 (exon-sequence alt-seq cds-start alt-exon-ranges*)]
    {:ref-exon-seq ref-exon-seq1
     :ref-prot-seq (codon/amino-acid-sequence (cond-> ref-exon-seq1
                                                (= strand :reverse) util-seq/revcomp))
     :alt-exon-seq alt-exon-seq1
     :alt-prot-seq (codon/amino-acid-sequence (cond-> alt-exon-seq1
                                                (= strand :reverse) util-seq/revcomp))
     :alt-tx-prot-seq (codon/amino-acid-sequence
                       (cond-> (str ref-up-exon-seq1 alt-exon-seq1 ref-down-exon-seq1)
                         (= strand :reverse) util-seq/revcomp))
     :ini-offset (quot (:position (rg/cds-coord (case strand
                                                  :forward tx-start
                                                  :reverse tx-end) rg))
                       3)}))

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

(defn- ->protein-variant
  [{:keys [strand] :as rg} pos ref alt
   {:keys [ref-exon-seq ref-prot-seq alt-exon-seq] :as seq-info}
   {:keys [prefer-deletion?]}]
  (cond
    (= ref-exon-seq alt-exon-seq)
    {:type :no-effect, :pos 1, :ref nil, :alt nil}

    (pos? (mod (count ref-exon-seq) 3))
    (do (log/warnf "CDS length is indivisible by 3: %d (%s, %s)"
                   (count ref-exon-seq) (:name rg) (:name2 rg))
        {:type :unknown, :pos nil, :ref nil, :alt nil})

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
                                          :forward (protein-position (+ pos (count alt) -1) rg)
                                          :reverse (protein-position (- pos (- (count alt) (count ref))) rg))
                                        (count alt-prot-seq*))])
          [pref-only palt-only offset _] (diff-bases pref palt)
          nprefo (count pref-only)
          npalto (count palt-only)
          [unit ref-repeat alt-repeat] (repeat-info* ref-prot-seq
                                                     (+ base-ppos offset)
                                                     pref-only
                                                     palt-only)
          t (cond
              (= (+ base-ppos offset) (count ref-prot-seq)) :extension
              (= (+ base-ppos offset) 1) (if (= ref-prot-rest alt-prot-rest)
                                           :extension
                                           :frame-shift)
              (and (pos? nprefo) (= (first palt-only) \*)) :substitution
              (not= ref-prot-rest alt-prot-rest) (if (or (and (empty? palt-only)
                                                              (= (first alt-prot-rest) \*))
                                                         (= (last palt-only) \*)) :fs-ter-substitution
                                                     :frame-shift)
              (or (and (zero? nprefo) (zero? npalto))
                  (and (= nprefo 1) (= npalto 1))) :substitution
              (and (some? unit) (= ref-repeat 1) (= alt-repeat 2)) :duplication
              (and prefer-deletion? (pos? nprefo) (zero? npalto)) :deletion
              (and (some? unit) (pos? alt-repeat)) :repeated-seqs
              (and (pos? nprefo) (zero? npalto)) :deletion
              (and (pos? nprefo) (pos? npalto)) (if (= base-ppos 1)
                                                  :extension
                                                  :indel)
              (and (zero? nprefo) (pos? npalto)) :insertion
              :else (throw (ex-info "Unsupported variant" {:type ::unsupported-variant})))]
      {:type (if (= t :fs-ter-substitution) :substitution t)
       :pos base-ppos
       :ref (if (= t :fs-ter-substitution)
              (str pref (subs ref-prot-rest 0 (max 0 (inc (- (count palt) (count pref))))))
              pref)
       :alt (if (= t :fs-ter-substitution)
              (str palt \*)
              palt)})))

(defn- protein-substitution
  [ppos pref palt]
  (let [[s-ref s-alt offset _] (diff-bases pref palt)]
    (if (and (empty? s-ref) (empty? s-alt))
      (mut/protein-substitution (mut/->long-amino-acid (last pref))
                                (coord/protein-coordinate ppos)
                                (mut/->long-amino-acid (last palt)))
      (mut/protein-substitution (mut/->long-amino-acid (first s-ref))
                                (coord/protein-coordinate (+ ppos offset))
                                (mut/->long-amino-acid (first s-alt))))))

(defn- protein-deletion
  [ppos pref palt]
  (let [[del _ offset _] (diff-bases pref palt)
        ndel (count del)]
    (mut/protein-deletion (mut/->long-amino-acid (first del))
                          (coord/protein-coordinate (+ ppos offset))
                          (if (> ndel 1)
                            (mut/->long-amino-acid (last del)))
                          (if (> ndel 1)
                            (coord/protein-coordinate (+ ppos offset ndel -1))))))

(defn- protein-duplication
  [ppos pref palt]
  (let [[_ ins offset _] (diff-bases pref palt)
        nins (count ins)]
    (mut/protein-duplication (mut/->long-amino-acid (first ins))
                             (coord/protein-coordinate (- (+ ppos offset) nins))
                             (if (> nins 1) (mut/->long-amino-acid (last ins)))
                             (if (> nins 1) (coord/protein-coordinate (dec (+ ppos offset)))))))

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

(defn- protein-indel
  [ppos pref palt]
  (let [[del ins offset _] (diff-bases pref palt)
        ndel (count del)]
    (mut/protein-indel (mut/->long-amino-acid (first del))
                       (coord/protein-coordinate (+ ppos offset))
                       (if (> ndel 1)
                         (mut/->long-amino-acid (last del)))
                       (if (> ndel 1)
                         (coord/protein-coordinate (+ ppos offset ndel -1)))
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
                               (if (< start end) (mut/->long-amino-acid (last unit)))
                               (if (< start end) (coord/protein-coordinate end))
                               alt-repeat)))

(defn- protein-frame-shift
  [ppos pref palt seq-info]
  (let [[ppos pref palt] (if (= pref palt)
                           (->> (map vector (:ref-prot-seq seq-info) (:alt-prot-seq seq-info))
                                (drop (dec ppos))
                                (map-indexed vector)
                                (drop-while (fn [[_ [r a]]] (= r a)))
                                first
                                ((fn [[i [r a]]]
                                   [(+ ppos i) (str r) (str a)])))
                           [ppos pref palt])
        [_ _ offset _] (diff-bases pref palt)
        alt-prot-seq (format-alt-prot-seq seq-info)
        ref (nth (:ref-prot-seq seq-info) (dec (+ ppos offset)))
        alt (nth alt-prot-seq (dec (+ ppos offset)))
        ter-site (-> seq-info
                     format-alt-prot-seq
                     (subs (dec (+ ppos offset)))
                     (string/index-of "*"))]
    (mut/protein-frame-shift (mut/->long-amino-acid ref)
                             (coord/protein-coordinate (+ ppos offset))
                             (mut/->long-amino-acid alt)
                             (if (and ter-site (pos? ter-site))
                               (coord/protein-coordinate (inc ter-site))
                               (coord/unknown-coordinate)))))

(defn- protein-extension
  [ppos pref palt seq-info]
  (let [{:keys [alt-prot-seq alt-tx-prot-seq ini-offset]} seq-info
        [_ ins offset _] (diff-bases pref palt)
        rest-seq (if (= ppos 1)
                   (-> alt-tx-prot-seq
                       (subs 0 ini-offset)
                       reverse
                       (#(apply str %)))
                   (-> alt-tx-prot-seq
                       (subs (+ ini-offset (count alt-prot-seq)))))
        ter-site (some-> (string/index-of rest-seq (if (= ppos 1) "M" "*")) inc)]
    (mut/protein-extension (if (= ppos 1) "Met" "Ter")
                           (coord/protein-coordinate (if (= ppos 1) 1 (+ ppos offset)))
                           (mut/->long-amino-acid (last ins))
                           (if (= ppos 1) :upstream :downstream)
                           (if ter-site
                             (coord/protein-coordinate ter-site)
                             (coord/unknown-coordinate)))))

(defn- mutation
  [seq-rdr rg pos ref alt options]
  (let [seq-info (read-sequence-info seq-rdr rg pos ref alt)]
    (if-let [{mut-type :type, ppos :pos, pref :ref, palt :alt}
             (->protein-variant rg pos ref alt seq-info options)]
      (case mut-type
        :substitution (protein-substitution ppos pref palt)
        :deletion (protein-deletion ppos pref palt)
        :duplication (protein-duplication ppos pref palt)
        :insertion (protein-insertion ppos pref palt seq-info)
        :indel (protein-indel ppos pref palt)
        :repeated-seqs (protein-repeated-seqs ppos pref palt seq-info)
        :frame-shift (protein-frame-shift ppos pref palt seq-info)
        :extension (protein-extension ppos pref palt seq-info)
        :no-effect (mut/protein-no-effect)
        :unknown (mut/protein-unknown-mutation)))))

(defn ->hgvs
  ([variant seq-rdr rg]
   (->hgvs variant seq-rdr rg {}))
  ([{:keys [pos ref alt]} seq-rdr rg options]
   (if-let [mutation (mutation seq-rdr rg pos ref alt options)]
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
