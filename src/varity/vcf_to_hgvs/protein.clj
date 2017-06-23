(ns varity.vcf-to-hgvs.protein
  (:require [clojure.string :as string]
            [clj-hgvs.coordinate :as coord]
            [clj-hgvs.core :as hgvs]
            [clj-hgvs.mutation :as mut]
            [cljam.fasta :as fa]
            [varity.codon :as codon]
            [varity.ref-gene :as rg]
            [varity.util :refer [revcomp-bases]]
            [varity.vcf-to-hgvs.common :refer [diff-bases] :as common]))

(defn- alt-sequence
  "Returns sequence a variant applied."
  [ref-seq seq-start pos ref alt]
  (let [pos* (inc (- pos seq-start))]
    (str (subs ref-seq 0 (dec pos*))
         alt
         (subs ref-seq (+ (dec pos*) (count ref))))))

(defn- alt-exon-ranges
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

(defn- exon-sequence
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
  [fa-rdr chr start end exon-ranges]
  (exon-sequence (fa/read-sequence fa-rdr {:chr chr, :start start, :end end})
                 start end exon-ranges))

(defn- read-sequence-info
  [fa-rdr rg pos ref alt]
  (let [{:keys [chr tx-start tx-end cds-start cds-end exon-ranges strand]} rg
        ref-seq (fa/read-sequence fa-rdr {:chr chr, :start cds-start, :end cds-end})
        alt-seq (alt-sequence ref-seq cds-start pos ref alt)
        alt-exon-ranges* (alt-exon-ranges exon-ranges pos ref alt)
        ref-exon-seq1 (exon-sequence ref-seq cds-start exon-ranges)
        ref-up-exon-seq1 (->> (read-exon-sequence fa-rdr chr tx-start (dec cds-start) exon-ranges)
                              reverse
                              (partition 3)
                              flatten
                              reverse
                              (apply str))
        ref-down-exon-seq1 (->> (read-exon-sequence fa-rdr chr (inc cds-end) tx-end exon-ranges)
                                reverse
                                (partition 3)
                                flatten
                                reverse
                                (apply str))
        alt-exon-seq1 (exon-sequence alt-seq cds-start alt-exon-ranges*)]
    {:ref-prot-seq (codon/amino-acid-sequence (cond-> ref-exon-seq1
                                             (= strand "-") revcomp-bases))
     :alt-prot-seq (codon/amino-acid-sequence (cond-> alt-exon-seq1
                                             (= strand "-") revcomp-bases))
     :alt-tx-prot-seq (codon/amino-acid-sequence
                       (cond-> (str ref-up-exon-seq1 alt-exon-seq1 ref-down-exon-seq1)
                         (= strand "-") revcomp-bases))
     :ini-offset (quot (:position (rg/cds-coord (case strand
                                                  "+" tx-start
                                                  "-" tx-end) rg))
                       3)}))

(defn- protein-position
  ([pos rg] (protein-position pos 0 rg))
  ([pos offset rg]
   (let [cds-coord (rg/cds-coord pos rg)]
     (if (coord/in-exon? cds-coord)
       (inc (quot (dec (+ (:position cds-coord) offset)) 3))))))

(defn- ->protein-variant
  [rg pos ref alt seq-info]
  (if-let [ppos (protein-position pos rg)]
    (let [{:keys [strand]} rg
          {:keys [ref-prot-seq alt-prot-seq]} seq-info
          base-ppos (case strand
                      "+" ppos
                      "-" (protein-position pos (- (count ref)) rg))
          pref (subs ref-prot-seq
                     (dec base-ppos)
                     (case strand
                       "+" (protein-position pos (dec (count ref)) rg)
                       "-" ppos))
          palt (subs alt-prot-seq
                     (dec base-ppos)
                     (case strand
                       "+" (protein-position pos (dec (count alt)) rg)
                       "-" (protein-position pos (- (count alt) (count ref)) rg)))
          [pref-only palt-only offset _] (diff-bases pref palt)
          nprefo (count pref-only)
          npalto (count palt-only)
          [unit ref-repeat ins-repeat] (common/repeat-info ref-prot-seq (+ base-ppos offset) palt-only)]
      {:type (cond
               (or (= (+ base-ppos offset) 1)
                   (= (+ base-ppos offset) (count ref-prot-seq))) :extension
               (not= (subs ref-prot-seq (+ base-ppos (count pref) -1))
                     (subs alt-prot-seq (+ base-ppos (count palt) -1))) :frame-shift
               (or (and (zero? nprefo) (zero? npalto))
                   (and (= nprefo 1) (= npalto 1))) :substitution
               (and (pos? nprefo) (zero? npalto)) :deletion
               (and (pos? nprefo) (pos? npalto)) (if (= base-ppos 1)
                                                   :extension
                                                   :indel)
               (some? unit) (cond
                              (and (= ref-repeat 1) (= ins-repeat 1)) :duplication
                              (or (> ref-repeat 1) (> ins-repeat 1)) :repeated-seqs)
               (and (zero? nprefo) (pos? npalto)) :insertion
               :else (throw (IllegalArgumentException. "Unsupported variant")))
       :pos base-ppos
       :ref pref
       :alt palt})))

(defn- protein-substitution
  [ppos pref palt]
  (let [[_ _ offset _] (diff-bases pref palt)]
    (mut/protein-substitution (mut/->long-amino-acid (last pref))
                              (coord/protein-coordinate (+ ppos offset))
                              (mut/->long-amino-acid (last palt)))))

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
  (let [[_ ins offset _] (diff-bases pref palt)
        [unit ref-repeat ins-repeat] (common/repeat-info (:ref-prot-seq seq-info) (+ ppos offset) ins)
        nunit (count unit)
        start (+ ppos offset (- (* nunit ref-repeat)))
        end (dec (+ start nunit))]
    (mut/protein-repeated-seqs (mut/->long-amino-acid (first ins))
                               (coord/protein-coordinate start)
                               (if (< start end) (mut/->long-amino-acid (last ins)))
                               (if (< start end) (coord/protein-coordinate end))
                               (+ ref-repeat ins-repeat))))

(defn- protein-frame-shift
  [ppos pref palt seq-info]
  (let [[ppos pref palt] (if (= pref palt)
                           (->> (map vector (:ref-prot-seq seq-info) (:alt-prot-seq seq-info))
                                (drop (dec ppos))
                                (map-indexed #(vector %1 %2))
                                (drop-while (fn [[_ [r a]]] (= r a)))
                                first
                                ((fn [[i [r a]]]
                                   [(+ ppos i) (str r) (str a)])))
                           [ppos pref palt])
        [del ins offset _] (diff-bases pref palt)
        ref (nth (:ref-prot-seq seq-info) (dec (+ ppos offset)))
        alt (nth (:alt-prot-seq seq-info) (dec (+ ppos offset)))
        ter-site (-> (:alt-prot-seq seq-info)
                     (subs (dec (+ ppos offset)))
                     (string/index-of "*"))]
    (mut/protein-frame-shift (mut/->long-amino-acid ref)
                             (coord/protein-coordinate (+ ppos offset))
                             (mut/->long-amino-acid alt)
                             (if ter-site
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
  [fa-rdr rg pos ref alt]
  (let [seq-info (read-sequence-info fa-rdr rg pos ref alt)]
    (if-let [{mut-type :type, ppos :pos, pref :ref, palt :alt}
             (->protein-variant rg pos ref alt seq-info)]
      (case mut-type
        :substitution (protein-substitution ppos pref palt)
        :deletion (protein-deletion ppos pref palt)
        :duplication (protein-duplication ppos pref palt)
        :insertion (protein-insertion ppos pref palt seq-info)
        :indel (protein-indel ppos pref palt)
        :repeated-seqs (protein-repeated-seqs ppos pref palt seq-info)
        :frame-shift (protein-frame-shift ppos pref palt seq-info)
        :extension (protein-extension ppos pref palt seq-info)))))

(defn ->hgvs
  [{:keys [pos ref alt]} fa-rdr rg]
  (if-let [mutation (mutation fa-rdr rg pos ref alt)]
    (hgvs/hgvs nil :protein mutation)))
