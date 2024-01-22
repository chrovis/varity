# varity

[![Clojars Project](https://img.shields.io/clojars/v/varity.svg)](https://clojars.org/varity)
[![build](https://github.com/chrovis/varity/actions/workflows/build.yml/badge.svg)](https://github.com/chrovis/varity/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/chrovis/varity/branch/master/graph/badge.svg)](https://codecov.io/gh/chrovis/varity)

Variant translation library for Clojure.

## Features

* VCF variant ⇄ HGVS
* Finding HGVS aliases
* Conversion of a genomic coordinates between assemblies

## Installation

Clojure CLI/deps.edn:

```clojure
varity/varity {:mvn/version "0.11.0"}
```

Leiningen/Boot:

```clojure
[varity "0.11.0"]
```

To use varity with Clojure 1.8, you must include a dependency on
[clojure-future-spec](https://github.com/tonsky/clojure-future-spec).

## Breaking changes in 0.11.0

We fixed following `varity.vcf-to-hgvs` implementation.
It is inconvenient to throw an exception in `vcf-variant->protein-hgvs` when variants overlapping a boundary of exon/intron even if `vcf-variant->coding-dna-hgvs` succeed in `vcf-to-hgvs.vcf-variant->hgvs`. So we changed to return nil.

## Breaking changes in 0.10.0

We introduced enhancements to the description of protein changes by `varity.vcf-to-hgvs`, specifically making deletions more clinically meaningful:

1. exon-intron boundary deletions:

The deletion that overlaps the exon-intron boundary will trigger an Exception because alterations affecting the splice sites are predicted to be splicing abnormalities.

2. stop codon deletions:

In cases where deletions contain a stop codon, `varity.vcf-to-hgvs` generates the following outputs based on the alteration sequence:

- If the alteration sequence contains a stop codon, varity outputs as deletion-insertion.
- Otherwise, this outputs `p.?`.

## Breaking changes in 0.9.0

The default value of `:prefer-deletion?` option is changed to `false`.

```clojure
(require '[varity.vcf-to-hgvs :as v2h])

(v2h/vcf-variant->coding-dna-hgvs {:chr "chr7", :pos 140924774, :ref "GGGAGGC", :alt "G"}
                                  "path/to/hg38.fa" "path/to/refGene.txt.gz")
;;=> (#clj-hgvs/hgvs "NM_004333:c.-95GCCTCC[3]")
```

If you hope the previous behavior, specify `:prefer-deletion? true`.

```clojure
(v2h/vcf-variant->coding-dna-hgvs {:chr "chr7", :pos 140924774, :ref "GGGAGGC", :alt "G"}
                                  "path/to/hg38.fa" "path/to/refGene.txt.gz"
                                  {:prefer-deletion? true})
;;=> (#clj-hgvs/hgvs "NM_004333:c.-77_-72delGCCTCC")
```

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
;;=> ({:coding-dna #clj-hgvs/hgvs "NM_005228:c.2573T>G",
;;     :protein #clj-hgvs/hgvs "p.L858R"})
```

Use `clj-hgvs.core/format` to obtain HGVS text.

```clojure
(require '[clj-hgvs.core :as hgvs])

(def l858r (-> (v2h/vcf-variant->protein-hgvs {:chr "chr7", :pos 55191822, :ref "T", :alt "G"}
                                              "path/to/hg38.fa" "path/to/refGene.txt.gz")
               first))

(hgvs/format l858r {:amino-acid-format :long})
;;=> "p.Leu858Arg"
```

### HGVS to VCF variants

`varity.hgvs-to-vcf` provides functions to convert HGVS into VCF-style variants.

```clojure
(require '[varity.hgvs-to-vcf :as h2v]
         '[clj-hgvs.core :as hgvs])

(h2v/hgvs->vcf-variants #clj-hgvs/hgvs "NM_005228:c.2573T>G" "path/to/hg38.fa" "path/to/refGene.txt.gz")
;;=> ({:chr "chr7", :pos 55191822, :ref "T", :alt "G"})

(h2v/hgvs->vcf-variants #clj-hgvs/hgvs "c.2573T>G" "EGFR" "path/to/hg38.fa" "path/to/refGene.txt.gz")
;;=> ({:chr "chr7", :pos 55191822, :ref "T", :alt "G"})

(h2v/hgvs->vcf-variants #clj-hgvs/hgvs "p.A222V" "MTHFR" "path/to/hg38.fa" "path/to/refGene.txt.gz")
;;=> ({:chr "chr1", :pos 11796320, :ref "GG", :alt "CA"}
;;    {:chr "chr1", :pos 11796320, :ref "GG", :alt "AA"}
;;    {:chr "chr1", :pos 11796320, :ref "GG", :alt "TA"}
;;    {:chr "chr1", :pos 11796321, :ref "G", :alt "A"})
```

### Finding HGVS aliases

`varity.hgvs/find-aliases` finds alternative HGVS expressions for the same
variant.

```clojure
(require '[varity.hgvs :as vhgvs])

(vhgvs/find-aliases #clj-hgvs/hgvs "NM_000059:c.162CAA[1]"
                    "path/to/hg38.fa" "path/to/refGene.txt.gz")
;;=> (#clj-hgvs/hgvs "NM_000059:c.162CAA[1]"
;;    #clj-hgvs/hgvs "NM_000059:c.165_167delCAA")
```

### Conversion of a genomic coordinate between assemblies

To convert a genomic coordinate between assemblies,

```clojure
(require '[varity.lift :as lift])

(lift/convert-coord {:chr "chr1", :pos 743267} "path/to/hg19ToHg38.over.chain.gz")
;;=> {:chr "chr1", :pos 807887}
```

## License

Copyright 2017-2022 [Xcoo, Inc.](https://xcoo.jp/)

Licensed under the [Apache License, Version 2.0](LICENSE).

## Acknowledgements

The algorithm of `varity.fusion` was initially developed by Norio Tanaka at Cancer Precision Medicine Center, Japanese Foundation for Cancer Research. We thank him for his scientific insight and technical help.
