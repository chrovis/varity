(ns varity.util)

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
