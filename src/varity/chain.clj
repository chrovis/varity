(ns varity.chain
  "Functions to read, index, and search the chain format. See
  https://genome.ucsc.edu/goldenPath/help/chain.html for the detail chain format
  specifications."
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [cljam.util.chromosome :refer [normalize-chromosome-key]]
            [proton.core :refer [as-long]])
  (:import [org.apache.commons.compress.compressors
            CompressorStreamFactory CompressorException]))

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
  (-> (zipmap [:score
               :t-name :t-size :t-strand :t-start :t-end
               :q-name :q-size :q-strand :q-start :q-end
               :id]
              (drop 1 (string/split line #"\s")))
      (update-multi [:score
                     :t-size :t-start :t-end
                     :q-size :q-start :q-end
                     :id]
                    as-long)))

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
  (with-open [^java.io.Reader rdr (try (-> (CompressorStreamFactory.)
                                           (.createCompressorInputStream (io/input-stream f))
                                           io/reader)
                                       (catch CompressorException _ (io/reader f)))]
    (->> (line-seq rdr)
         split-lines
         doall)))

;;; Indexing

(defn index
  "Creates chain index for search."
  [chains]
  (->> (group-by (comp :t-name :header) chains)
       doall))

;;; Search

(defn- in-block?
  [pos {:keys [header data]}]
  (let [rpos (- pos (:t-start header))]
    (and (>= rpos 0)
         (loop [[{:keys [^long size ^long dt]} & r] data, s 0]
           (let [e (+ s size)]
             (if (< rpos s)
               false
               (if (< rpos e)
                 true
                 (if (nil? dt)
                   false
                   (recur r (+ e dt))))))))))

(defn search-chains
  "Searches chain entries with chr and pos using the index, returning the
  results as a sequence. See also varity.chain/index."
  [chr pos chain-idx]
  (->> (get chain-idx (normalize-chromosome-key chr))
       (filter (partial in-block? pos))))
