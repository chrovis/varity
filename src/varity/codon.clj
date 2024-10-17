(ns varity.codon
  "Handles codon."
  (:require [clojure.string :as string]))

;; See https://en.wikipedia.org/wiki/DNA_codon_table
(def ^:private codon-amino-acid-map*
  {"TTT" "F"
   "TTC" "F"
   "TTA" "L"
   "TTG" "L"
   "CTT" "L"
   "CTC" "L"
   "CTA" "L"
   "CTG" "L"
   "ATT" "I"
   "ATC" "I"
   "ATA" "I"
   "ATG" "M"
   "GTT" "V"
   "GTC" "V"
   "GTA" "V"
   "GTG" "V"
   "TCT" "S"
   "TCC" "S"
   "TCA" "S"
   "TCG" "S"
   "CCT" "P"
   "CCC" "P"
   "CCA" "P"
   "CCG" "P"
   "ACT" "T"
   "ACC" "T"
   "ACA" "T"
   "ACG" "T"
   "GCT" "A"
   "GCC" "A"
   "GCA" "A"
   "GCG" "A"
   "TAT" "Y"
   "TAC" "Y"
   "TAA" "*"
   "TAG" "*"
   "CAT" "H"
   "CAC" "H"
   "CAA" "Q"
   "CAG" "Q"
   "AAT" "N"
   "AAC" "N"
   "AAA" "K"
   "AAG" "K"
   "GAT" "D"
   "GAC" "D"
   "GAA" "E"
   "GAG" "E"
   "TGT" "C"
   "TGC" "C"
   "TGA" "*"
   "TGG" "W"
   "CGT" "R"
   "CGC" "R"
   "CGA" "R"
   "CGG" "R"
   "AGT" "S"
   "AGC" "S"
   "AGA" "R"
   "AGG" "R"
   "GGT" "G"
   "GGC" "G"
   "GGA" "G"
   "GGG" "G"})

(def ^:private uncertain-codon-map
  (let [bases ["A" "T" "G" "C" "N"]
        uncertain-base "N"
        other-base-combinations (for [a bases b bases] (vector a b))
        first-uncertain-base-codons (map #(str uncertain-base (first %) (second %)) other-base-combinations)
        second-uncertain-base-codons (map #(str (first %) uncertain-base (second %)) other-base-combinations)
        third-uncertain-base-codons (map #(str (first %) (second %) uncertain-base) other-base-combinations)
        uncertain-codons (set (concat first-uncertain-base-codons second-uncertain-base-codons third-uncertain-base-codons))]
    (into {} (map #(vector % "X") uncertain-codons))))

(def codon-amino-acid-map (merge codon-amino-acid-map* uncertain-codon-map))

(def ^:private amino-acid-codon-map
  (->> (group-by second codon-amino-acid-map)
       (map #(update % 1 (partial mapv first)))
       (into {})))

(defn codon->amino-acid
  "Converts three-character nucleotide string s into a single-letter amino acid.

  e.g.
    (genes->codon \"TCA\") => \"S\"
    (genes->codon \"tca\") => \"S\"
    (genes->codon \"TAA\") => \"*\"
    (genes->codon \"TCAG\") => AssertionError"
  [s]
  {:pre [(= (count s) 3)]}
  (get codon-amino-acid-map (string/upper-case s)))

(defn amino-acid->codons
  "Returns possible nucleotide strings corresponding to amino acid s. s must be
  one-character String or Character.

  e.g.
    (codon->genes \"L\") => [\"CTG\" \"CTC\" \"CTT\" \"CTA\" \"TTG\" \"TTA\"]
    (codon->genes \\*) => [\"TAG\" \"TAA\" \"TGA\"]"
  [s]
  {:pre [(or (and (string? s) (= (count s) 1)) (char? s))]}
  (let [s (cond-> s (char? s) str)]
    (get amino-acid-codon-map (string/upper-case s))))

(defn amino-acid-sequence
  "Converts nucleotide sequence into amino acid sequence."
  [s]
  (->> (string/upper-case s)
       (re-seq #".{3}")
       (map codon-amino-acid-map)
       (apply str)))
