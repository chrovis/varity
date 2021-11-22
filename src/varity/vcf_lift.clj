(ns varity.vcf-lift
  (:require [cljam.io.vcf.util.normalize :as norm]
            [cljam.io.vcf.util :as vcf-util]
            [cljam.io.vcf.util.check :as check]
            [cljam.util.sequence :as util-seq]
            [cljam.util.chromosome :as chr]
            [varity.chain :as ch]))

(defn- calc-new-interval [chr beg end indexed-chains]
  (some #(when-let [overlap-blocks
                    (seq (ch/search-overlap-blocks beg end (:data %)))]
           (let [start-offset (max 0 (- beg (:t-start (first overlap-blocks))))
                 end-offset (max 0 (- (:t-end (last overlap-blocks)) end))
                 to-start (+ (:q-start (first overlap-blocks)) start-offset)
                 to-end (- (:q-end (last overlap-blocks)) end-offset)
                 strand (get-in % [:header :q-strand])]
             (if (= strand :reverse)
               {:chr (get-in % [:header :q-name])
                :start (inc (- (get-in % [:header :q-size]) to-end))
                :end (inc (- (get-in % [:header :q-size]) to-start))
                :strand strand}
               {:chr (get-in % [:header :q-name])
                :start to-start :end to-end :strand strand})))
        (ch/search-containing-chains chr beg end indexed-chains)))

(defn liftover-variant*
  "`liftover-variant` without checking reference sequence and realignment."
  [indexed-chain {:keys [chr pos ref info] :as variant}]
  (when-let [{:keys [chr start strand]} (calc-new-interval
                                         chr pos (dec (+ pos (count ref)))
                                         indexed-chain)]
    (cond-> (assoc variant :chr chr :pos start)
      (:END info) (assoc-in [:info :END] (dec (+ start (count ref))))
      (= strand :reverse) (-> (update :ref util-seq/revcomp)
                              (update :alt #(map util-seq/revcomp %))))))

(defn liftover-variant
  "Lifts over a `variant` to another coordinate. Returns a converted and
   realigned variant if succeeded to lift over, otherwise nil.
   `indexed-chain` supplies conversions between sequences and must be indexed by
   `varity.chain/index`. `target-seq-reader` is a reader instance of the target
   sequence that the variant coordinate is converted to."
  [target-seq-reader indexed-chain variant]
  (when-let [x (liftover-variant* indexed-chain variant)]
    (when (check/same-ref? target-seq-reader x)
      (norm/realign target-seq-reader x))))

(defn liftover-variants
  "Lifts over variants data from chains applied chain/index and seq-reader.
  Variants that succeed in liftover will be put in :success,
  and those that failed will be put in :failure."
  [seq-reader indexed-chains variants]
  (let [separated-variants
        (group-by
         (fn [variant]
           (->> (:alt variant)
                (map #(:type (vcf-util/inspect-allele (:ref variant) %)))
                (every? #{:mnv :ref :complex :insertion :deletion :snv})))
         variants)]
    (->> (get separated-variants true)
         (reduce (fn [m x]
                   (if-let [new-variant (liftover-variant seq-reader
                                                          indexed-chains x)]
                     (update m :success conj new-variant)
                     (update m :failure conj x)))
                 {:success [], :failure (vec (get separated-variants false))})
         (reduce-kv #(->> %3
                          (sort-by (juxt (comp chr/chromosome-order-key :chr)
                                         :pos :ref (comp vec :alt)))
                          seq
                          (assoc %1 %2))
                    {}))))
