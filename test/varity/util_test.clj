(ns varity.util-test
  (:require [clojure.test :refer :all]
            [varity.util :refer :all]))

(deftest revcomp-bases-test
  (are [s e] (= (revcomp-bases s) e)
    "ATGNC" "GNCAT"
    "atgnc" "gncat"
    "" "")
  (are [x] (thrown? Error (revcomp-bases x))
    nil
    2))
