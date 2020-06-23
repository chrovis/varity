(ns varity.chain
  "Functions to read, index, and search the chain format. See
  https://genome.ucsc.edu/goldenPath/help/chain.html for the detail chain format
  specifications."
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [cljam.util.chromosome :refer [normalize-chromosome-key]]
            [proton.core :refer [as-long]]
            [varity.util :as util])
  (:import [clojure.lang Sorted]))

(defn- update-multi
  [m ks f]
  (reduce #(update %1 %2 f) m ks))

;;; Chain file loading

(defn- header-line?
  [line]
  (if (some? line)
    (string/starts-with? line "chain")
    false))

;; chain score tName tSize    tStrand tStart   tEnd     qName qSize     qStrand qStart   qEnd     id
;; chain 4900  chrY  58368225 +       25985403 25985638 chr5  151006098 -       43257292 43257528 1
(defn- parse-header-line
  [line]
  (letfn [(parse-strand [^String s] (case (.charAt s 0) \+ :forward \- :reverse))]
    (-> (zipmap [:score
                 :t-name :t-size :t-strand :t-start :t-end
                 :q-name :q-size :q-strand :q-start :q-end
                 :id]
                (drop 1 (string/split line #"\s")))
        (update-multi [:t-strand :q-strand] parse-strand)
        (update-multi [:score
                       :t-size :t-start :t-end
                       :q-size :q-start :q-end
                       :id]
                      as-long))))

;; size dt dq
;; 16   0  2
;; 60   4  0
;; 10   0  4
;; 70
(defn- parse-data-line
  [line]
  (-> (zipmap [:size :dt :dq]
              (string/split line #"\s"))
      (update-multi [:size :dt :dq] as-long)))

(defn- split-lines
  [lines]
  (if (seq lines)
    (if (header-line? (first lines))
      (let [[t d] (split-with (complement header-line?) (rest lines))]
        (cons {:header (parse-header-line (first lines))
               :data (->> t (remove empty?) (mapv parse-data-line))}
              (lazy-seq (split-lines d))))
      (split-lines (rest lines)))))

(defn load-chain
  "Loads f (e.g. hg19ToHg38.over.chain(.gz)), returning the all contents as a
  sequence."
  [f]
  (with-open [rdr (io/reader (util/compressor-input-stream f))]
    (->> (line-seq rdr)
         split-lines
         doall)))

;;; Indexing

(defn- cumsum-chain [{:keys [header data]}]
  (loop [results (transient [])
         d (first data)
         r (next data)
         curr-t-start (inc (:t-start header))
         curr-q-start (inc (:q-start header))]
    (if-not d
      (persistent! results)
      (recur (->> (assoc d :t-start curr-t-start
                         :q-start curr-q-start
                         :t-end (+ curr-t-start (:size d))
                         :q-end (+ curr-q-start (:size d)))
                  (conj! results))
             (first r)
             (next r)
             (if (:dt d)
               (+ curr-t-start (:size d) (:dt d))
               (:t-end header))
             (if (:dq d)
               (+ curr-q-start (:size d) (:dq d))
               (:q-end header))))))

(defn- index-chain [chain]
  (->> chain
       cumsum-chain
       (map (juxt :t-start identity))
       (into (sorted-map))
       (assoc chain :data)))

(defn index
  "Creates chain index for search."
  [chains]
  (->> (group-by (comp :t-name :header) chains)
       (into {} (map
                 (fn [[chr xs]]
                   [(normalize-chromosome-key chr)
                    (mapv index-chain xs)])))))

;;; Search

(defn- in-block?
  [pos {:keys [^Sorted data header] :as chain}]
  (when (<= (inc (:t-start header)) pos (:t-end header))
    (when-let [[start m] (first (. data seqFrom pos false))]
      (when (<= start pos (dec (+ start (:size m))))
        (assoc chain :in-block m)))))

(def ^:private normalize-chr (memoize normalize-chromosome-key))

(defn search-chains
  "Searches chain entries with chr and pos using the index, returning the
  results as a sequence. See also varity.chain/index."
  [chr pos chain-idx]
  (->> (chain-idx (normalize-chr chr))
       (keep (partial in-block? pos))))

(defn search-containing-chains
  "Calculates a list of chains that contain the given interval."
  [chr start end chain-idx]
  (->> (chain-idx (normalize-chr chr))
       (filter #(and (<= (get-in % [:header :t-start]) start)
                     (<= end (get-in % [:header :t-end]))))
       (sort-by (comp :score :header) >)))

(defn search-overlap-blocks
  "Calculates a list of blocks that overlap the given interval."
  [start end blocks-idx]
  (->> (subseq blocks-idx <= end)
       (map second)
       (filter #(>= (:t-end %) start))))
