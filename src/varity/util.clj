(ns varity.util
  "Utilities."
  (:require [clojure.java.io :as io])
  (:import [java.io InputStream]
           [org.apache.commons.compress.compressors
            CompressorStreamFactory CompressorException]))

(def ^:private comp-base-map
  {\A \T
   \C \G
   \G \C
   \T \A
   \N \N
   \a \t
   \c \g
   \g \c
   \t \a
   \n \n})

(defn revcomp-bases
  "Returns reverse complementary bases string."
  [s]
  {:pre [(string? s)]}
  (->> (reverse s)
       (map (partial get comp-base-map))
       (apply str)))

(defn ^InputStream compressor-input-stream
  "Returns an compressor input stream from f, autodetecting the compressor type
  from the first few bytes of f. Returns java.io.BufferedInputStream if the
  compressor type is not known. Should be used inside with-open to ensure the
  InputStream is properly closed."
  [f]
  (let [is (io/input-stream f)]
    (try
      (-> (CompressorStreamFactory. true)
          (.createCompressorInputStream is))
      (catch CompressorException _
        is))))
