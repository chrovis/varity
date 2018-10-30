(ns varity.vcf-to-hgvs.genome
  (:require [clojure.pprint :as pp]
            [clojure.string :as string]
            [cljam.io.sequence :as cseq]
            [varity.vcf-to-hgvs.common :as common]))

(defn- sequence-pstring
  [seq* start end {:keys [pos ref alt]}]
  (let [[ref-up _ ref-down] (common/split-string-at seq*
                                                    [(- pos start)
                                                     (+ (- pos start) (count ref))])
        nmut (max (count ref) (count alt))
        ticks (->> (iterate #(+ % 10) (inc (- start (mod start 10))))
                   (take-while #(<= % end))
                   (remove #(< % start)))
        tick-intervals (conj
                        (->> ticks
                             (map #(+ % (if (< pos %) (max 0 (- (count alt) (count ref))) 0)))
                             (map #(- % start))
                             (partition 2 1)
                             (mapv #(- (second %) (first %))))
                        0)]
    (string/join
     \newline
     [(apply pp/cl-format nil
             (apply str "~V@T" (repeat (count ticks) "~VA"))
             (- (first ticks) start) (interleave tick-intervals ticks))
      (pp/cl-format nil "~A~VA~A" ref-up nmut ref ref-down)
      (pp/cl-format nil "~A~VA~A" ref-up nmut alt ref-down)
      (pp/cl-format nil "~V@T~V@{~A~:*~}" (- pos start) nmut "^")])))

(defn debug-string
  [{:keys [pos ref] :as var} seq-rdr rg]
  (let [start (max 1 (- pos 20))
        end (+ pos (count ref) 19)
        seq* (cseq/read-sequence seq-rdr {:chr (:chr rg), :start start, :end end})
        [ticks refs alts muts] (string/split-lines
                                (sequence-pstring seq* start end var))]
    (string/join \newline [(str "         " ticks)
                           (str "    ref: " refs)
                           (str "    alt: " alts)
                           (str "         " muts)])))
