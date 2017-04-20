# varity

[![Clojars Project](https://img.shields.io/clojars/v/varity.svg)](https://clojars.org/varity)
[![Build Status](https://travis-ci.org/chrovis/varity.svg?branch=master)](https://travis-ci.org/chrovis/varity)

Variant translation library for Clojure.

## Features

* VCF variant -> HGVS
* HGVS -> VCF variant
* Conversion between assemblies

## Installation

With Leiningen/Boot:

```clojure
[varity "0.2.0"]
```

## Usage

### Notice

All positions are represented as one-based number, and all ranges are
represented as one-based closed intervals. For example,

```clojure
{:pos 3}
```

represents the third position from the start, and

```clojure
{:chr "chr1", :start 1, :end 3}
```

represents the first three bases of chromosome 1.

### VCF variant to HGVS

`varity.vcf-to-hgvs` provides functions to convert a VCF-style variant into HGVS.
The returned HGVS is data structure of [clj-hgvs](https://github.com/chrovis/clj-hgvs).

```clojure
(require '[varity.vcf-to-hgvs :as v2h])

(v2h/vcf-variant->hgvs "path/to/hg38.fa" "path/to/refGene.txt.gz"
                       "chr7" 55191822 "T" "G")
;;=> ({:cdna {:kind :cdna,
;;            :mutations (#clj_hgvs.mutation.DNASubstitution
;;                        {:alt "G",
;;                         :coord #clj_hgvs.coordinate.CDNACoordinate
;;                          {:offset 0, :position 2573, :region nil},
;;                         :ref "T",
;;                         :type ">"}),
;;            :transcript "NM_005228"},
;;     :protein {:kind :protein,
;;               :mutations (#clj_hgvs.mutation.ProteinSubstitution
;;                           {:alt "Arg",
;;                            :coord #clj_hgvs.coordinate.ProteinCoordinate
;;                             {:position 858},
;;                            :ref "Leu"}),
;;               :transcript nil}})
```

Use `clj-hgvs.core/format` to obtain HGVS text.

```clojure
(require '[clj-hgvs.core :as hgvs])

(def l858r (-> (v2h/vcf-variant->protein-hgvs "path/to/hg38.fa"
                                              "path/to/refGene.txt.gz"
                                              "chr7" 55191822 "T" "G")
               first))

(hgvs/format l858r {:amino-acid-format :short})
;;=> p.L858R
```

### HGVS to VCF variants

`varity.hgvs-to-vcf` provides functions to convert HGVS into VCF-style variants.

```clojure
(require '[varity.hgvs-to-vcf :as h2v]
         '[clj-hgvs.core :as hgvs])

(h2v/hgvs->vcf-variants "path/to/hg38.fa" "path/to/refGene.txt.gz" (hgvs/parse "c.2573T>G") "EGFR")
;;=> ({:chr "chr7", :pos 55191822, :ref "T", :alt "G"})

(h2v/hgvs->vcf-variants "path/to/hg38.fa" "path/to/refGene.txt.gz" (hgvs/parse "p.A222V") "MTHFR")
;;=> ({:chr "chr1", :pos 11796321, :ref "G", :alt "A"})
```

### Conversion between assemblies

To convert a genome coordinate between assemblies,

```clojure
(require '[varity.lift :as lift])

(lift/convert-coord "chr1" 743267 "path/to/hg19ToHg38.over.chain.gz")
;;=> {:chr "chr1", :pos 807887}
```

## License

Copyright 2017 [Xcoo, Inc.](https://xcoo.jp/)

Licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).
