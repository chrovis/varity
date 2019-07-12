# varity

[![Clojars Project](https://img.shields.io/clojars/v/varity.svg)](https://clojars.org/varity)
[![Build Status](https://travis-ci.org/chrovis/varity.svg?branch=master)](https://travis-ci.org/chrovis/varity)
[![codecov](https://codecov.io/gh/chrovis/varity/branch/master/graph/badge.svg)](https://codecov.io/gh/chrovis/varity)

Variant translation library for Clojure.

## Features

* VCF variant ⇄ HGVS
* Conversion between assemblies

## Installation

Clojure CLI/deps.edn:

```clojure
varity {:mvn/version "0.5.1"}
```

Leiningen/Boot:

```clojure
[varity "0.5.1"]
```

To use varity with Clojure 1.8, you must include a dependency on
[clojure-future-spec](https://github.com/tonsky/clojure-future-spec).

## Usage

### Documentation

[API Reference](https://chrovis.github.io/varity/)

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

(v2h/vcf-variant->hgvs {:chr "chr7", :pos 55191822, :ref "T", :alt "G"}
                       "path/to/hg38.fa" "path/to/refGene.txt.gz")
;;=> ({:cdna {:kind :cdna,
;;            :mutation #clj_hgvs.mutation.DNASubstitution
;;            {:alt "G",
;;             :coord #clj_hgvs.coordinate.CDNACoordinate
;;             {:offset 0, :position 2573, :region nil},
;;             :ref "T",
;;             :type ">"},
;;            :transcript "NM_005228"},
;;     :protein {:kind :protein,
;;               :mutation #clj_hgvs.mutation.ProteinSubstitution
;;               {:alt "Arg",
;;                :coord #clj_hgvs.coordinate.ProteinCoordinate
;;                {:position 858},
;;                :ref "Leu"},
;;               :transcript nil}})
```

Use `clj-hgvs.core/format` to obtain HGVS text.

```clojure
(require '[clj-hgvs.core :as hgvs])

(def l858r (-> (v2h/vcf-variant->protein-hgvs {:chr "chr7", :pos 55191822, :ref "T", :alt "G"}
                                              "path/to/hg38.fa" "path/to/refGene.txt.gz")
               first))

(hgvs/format l858r {:amino-acid-format :short})
;;=> p.L858R
```

### HGVS to VCF variants

`varity.hgvs-to-vcf` provides functions to convert HGVS into VCF-style variants.

```clojure
(require '[varity.hgvs-to-vcf :as h2v]
         '[clj-hgvs.core :as hgvs])

(h2v/hgvs->vcf-variants (hgvs/parse "NM_005228:c.2573T>G") "path/to/hg38.fa" "path/to/refGene.txt.gz")
;;=> ({:chr "chr7", :pos 55191822, :ref "T", :alt "G"})

(h2v/hgvs->vcf-variants (hgvs/parse "c.2573T>G") "EGFR" "path/to/hg38.fa" "path/to/refGene.txt.gz")
;;=> ({:chr "chr7", :pos 55191822, :ref "T", :alt "G"})

(h2v/hgvs->vcf-variants (hgvs/parse "p.A222V") "MTHFR" "path/to/hg38.fa" "path/to/refGene.txt.gz")
;;=> ({:chr "chr1", :pos 11796321, :ref "G", :alt "A"})
```

### Conversion between assemblies

To convert a genome coordinate between assemblies,

```clojure
(require '[varity.lift :as lift])

(lift/convert-coord {:chr "chr1", :pos 743267} "path/to/hg19ToHg38.over.chain.gz")
;;=> {:chr "chr1", :pos 807887}
```

## License

Copyright 2017-2019 [Xcoo, Inc.](https://xcoo.jp/)

Licensed under the [Apache License, Version 2.0](LICENSE).
