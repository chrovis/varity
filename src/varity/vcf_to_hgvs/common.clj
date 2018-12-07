(ns varity.vcf-to-hgvs.common
  (:require [clojure.string :as string]
            [cljam.io.sequence :as cseq]
            [varity.ref-gene :as rg]))

(defn- left-common-len
  "Returns the number of common characters on the left side of s1 and s2."
  [s1 s2]
  (->> (map #(if (= %1 %2) 1 0) s1 s2)
       (take-while #(= % 1))
       (reduce +)))

(defn diff-bases
  "Compares VCF-style base(s) s1 and s2, returning a vector of
  [bases-only-in-s1 bases-only-in-s2 left-common-len right-common-len]."
  [s1 s2]
  (let [l (left-common-len s1 s2)
        ls1 (subs s1 l)
        ls2 (subs s2 l)
        r (left-common-len (reverse ls1) (reverse ls2))]
    [(subs ls1 0 (- (count ls1) r))
     (subs ls2 0 (- (count ls2) r))
     l
     r]))

(def max-repeat-unit-size 100)

;; AGT => ["AGT"]
;; AGTAGT => ["AGT" "AGTAGT"]
(defn- repeat-units
  [s]
  (case (count s)
    0 []
    1 [s]
    (let [limit (min (/ (count s) 2) max-repeat-unit-size)
          xf (comp (map #(partition-all % s))
                   (filter #(apply = %))
                   (map first)
                   (map #(apply str %)))]
      (conj (into [] xf (range 1 (inc limit))) s))))

;; seq*        pos bases
;; MVSTSTHQ... 3   ST    => 2
(defn forward-shift
  [seq* pos bases]
  (if (= (subs seq* (dec pos) (+ (dec pos) (count bases))) bases)
    (if (string/blank? bases)
      0
      (let [step (count (first (repeat-units bases)))]
        (->> seq*
             (drop (dec pos))
             (partition (count bases) step)
             (take-while #(= (apply str %) bases))
             count
             dec
             (* step))))
    (throw (ex-info "The bases is not found on the position."
                    {:type ::invalid-bases
                     :sequence seq*
                     :position pos
                     :bases bases}))))

;; seq*        pos bases
;; MVSTSTHQ... 5   ST    => 2
(defn backward-shift
  [seq* pos bases]
  (if (= (subs seq* (dec pos) (+ (dec pos) (count bases))) bases)
    (if (string/blank? bases)
      0
      (let [rbases (string/reverse bases)
            n (count bases)
            step (count (first (repeat-units rbases)))
            tweak (if (= step 1) #(+ (dec n) %) #(* step %))]
        (->> seq*
             (take (dec pos))
             reverse
             (partition (count rbases) step)
             (take-while #(= (apply str %) rbases))
             count
             tweak)))
    (throw (ex-info "The bases is not found on the position."
                    {:type ::invalid-bases
                     :sequence seq*
                     :position pos
                     :bases bases}))))

;; ...CAGTC... 8 AGT => ["AGT" 1 1]
;; ...CAGTC... 8 AGTAGT => ["AGT" 1 2]
;; ...CAGTAGTC... 11 AGTAGT => ["AGT" 2 2]
(defn repeat-info
  [seq* pos ins]
  (->> (repeat-units ins)
       (map (fn [unit]
              [unit
               (->> (subs seq* 0 (dec pos))
                    (reverse)
                    (partition (count unit))
                    (take-while #(= % (reverse unit)))
                    (count))
               (/ (count ins) (count unit))]))
       (remove (comp zero? second))
       (first)))

;; ("AGT" "AGTCTG" "AGTCTGAAA" ...)
(defn read-sequence-stepwise
  [seq-rdr {:keys [chr start end]} step]
  (->> (range)
       (map (fn [i]
              (let [start* (+ start (* i step))
                    end* (min (dec (+ start* step)) end)]
                (if (<= start* end)
                  (cseq/read-sequence seq-rdr {:chr chr, :start start, :end end*})))))
       (take-while some?)))

;; ("AGT" "CTGAGT" "AAACTGAGT" ...)
(defn read-sequence-stepwise-backward
  [seq-rdr {:keys [chr start end]} step]
  (->> (range)
       (map (fn [i]
              (let [end* (- end (* i step))
                    start* (max (inc (- end* step)) start)]
                (if (>= end* start)
                  (cseq/read-sequence seq-rdr {:chr chr, :start start*, :end end})))))
       (take-while some?)))

;; + ...CAGTAGTAGTC... 7 T TAGT => 13 T TAGT
;; + ...CAGTAGTAGTC... 7 TAGT T => 10 TAGT T
;; - ...CAGTAGTAGTC... 7 T TAGT => 4 C CAGT
;; - ...CAGTAGTAGTC... 7 TAGT T => 4 CAGT C
(defn- normalize-variant*
  [{:keys [pos ref alt] :as var} seq* strand]
  (let [[ref-only alt-only offset _] (diff-bases ref alt)
        type* (cond
                (and (zero? (count ref-only)) (pos? (count alt-only))) :ins
                (and (pos? (count ref-only)) (zero? (count alt-only))) :del
                :else :other)]
    (if (#{:ins :del} type*)
      (let [bases (if (string/blank? ref-only) alt-only ref-only)
            seq' (if (= type* :ins)
                   (str (subs seq* 0 (+ pos offset -1)) bases (subs seq* (+ pos offset -1)))
                   seq*)
            move (case strand
                   :forward (forward-shift seq' (+ pos offset) bases)
                   :reverse (- (backward-shift seq' (+ pos offset) bases)))
            comm (subs seq* (+ pos move -1) (+ pos move offset -1))]
        {:pos (+ pos move)
         :ref (str comm ref-only)
         :alt (str comm alt-only)})
      var)))

(def normalization-read-sequence-step 1000)

(defn- normalize-variant-forward
  [{:keys [chr pos ref alt]} seq-rdr rg]
  (->> (read-sequence-stepwise seq-rdr
                               {:chr chr
                                :start pos
                                :end (+ (:tx-end rg) rg/max-tx-margin)}
                               normalization-read-sequence-step)
       (keep (fn [seq*]
               (when (>= (count seq*) (count ref))
                 (let [nvar (normalize-variant* {:pos 1, :ref ref, :alt alt} seq* :forward)]
                   (if (<= (dec (:pos nvar)) (- (count seq*) (max (count ref) (count alt))))
                     (-> nvar
                         (assoc :chr chr)
                         (update :pos + pos -1)))))))
       (first)))

(defn- normalize-variant-backward
  [{:keys [chr pos ref alt]} seq-rdr rg]
  (->> (read-sequence-stepwise-backward seq-rdr
                                        {:chr chr
                                         :start (- (:tx-start rg) rg/max-tx-margin)
                                         :end (+ pos (count ref) -1)}
                                        normalization-read-sequence-step)
       (keep (fn [seq*]
               (let [offset (- (count seq*) (count ref))]
                 (when (>= offset 0)
                   (let [nvar (normalize-variant* {:pos (inc offset), :ref ref, :alt alt} seq* :reverse)]
                     (if (> (:pos nvar) (max (count ref) (count alt)))
                       (-> nvar
                           (assoc :chr chr)
                           (update :pos + (- pos offset) -1))))))))
       (first)))

(defn normalize-variant
  "Normalizes a VCF-style variant based on surrounding sequence. For example,
  {:pos 7, :ref T, :alt TAGT} on ...CAGTAGTAGTC... is equivalent to
  {:pos 13, :ref T, :alt TAGT}. The latter is normalized.
  e.g.  ...CAGTAGTAGTC... 7 T TAGT => 13 T TAGT"
  [variant seq-rdr rg]
  (case (:strand rg)
    :forward (normalize-variant-forward variant seq-rdr rg)
    :reverse (normalize-variant-backward variant seq-rdr rg)))

(defn alt-sequence
  "Returns sequence a variant applied."
  [ref-seq seq-start pos ref alt]
  (let [pos* (inc (- pos seq-start))]
    (str (subs ref-seq 0 (max 0 (dec pos*)))
         alt
         (subs ref-seq (min (count ref-seq)
                            (+ (dec pos*) (count ref)))))))
