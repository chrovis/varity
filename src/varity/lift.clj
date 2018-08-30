(ns varity.lift
  "Functions to convert a genome coordinate between assemblies."
  (:require [varity.chain :as ch]))

(defn- convert-coord*
  [pos {{:keys [q-name q-strand q-size]} :header
        {:keys [q-start t-start]} :in-block}]
  (let [new-pos (+ q-start (- pos t-start))]
    {:chr q-name
     :pos (if (= q-strand :reverse)
            (inc (- q-size new-pos))
            new-pos)}))

(defn convert-coord
  "Converts {:chr :pos} between different assemblies, returning a new
  {:chr :pos}. The last argument, chain, must be a path to
  srcToDest.over.chain(.gz) or chain index."
  [{:keys [chr pos]} chain]
  (let [chidx (if (map? chain)
                chain
                (ch/index (ch/load-chain chain)))]
    (some->> (ch/search-chains chr pos chidx)
             (sort-by (comp :score :header) >)
             first
             (convert-coord* pos))))
