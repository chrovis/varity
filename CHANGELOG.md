# Changelog

## [0.12.0] - 2025-07-29

### BREAKING

Dropped support for clojure 1.8 and added support for clojure 1.12 in [#123](https://github.com/chrovis/varity/pull/123).

### Changed

- Refactor and generalize overlap-exon-intron-boundary?. [#93](https://github.com/chrovis/varity/pull/93)
- Add no-effect branch to protein type determine process. [#98](https://github.com/chrovis/varity/pull/98)
- Add prefer-extension option for variant affects initial codon. [#101](https://github.com/chrovis/varity/pull/101)
- Introduce GeneAnnotationIndex as an abstraction of RefGeneIndex. [#102](https://github.com/chrovis/varity/pull/102)
- Add UTR variant determine process. [#104](https://github.com/chrovis/varity/pull/104)
- Support uncertain base and amino acid. [#109](https://github.com/chrovis/varity/pull/109)
- Return unknown when initiation site is affected. [#111](https://github.com/chrovis/varity/pull/111)
- Update clj-hgvs to 0.5.1. [#112](https://github.com/chrovis/varity/pull/112)
- Return no-effect when ter site is inserted around ter site and fix extension determine process. [#115](https://github.com/chrovis/varity/pull/115)
- Add in-frame check process to indel for variant around ter site. [#117](https://github.com/chrovis/varity/pull/117)
- Add three-prime-rule option to vcf->hgvs. [#118](https://github.com/chrovis/varity/pull/118)
- Add extension conditional branch to insertion. [#119](https://github.com/chrovis/varity/pull/119)
- Add prefer-deletion-insertion option to coding-dna. [#120](https://github.com/chrovis/varity/pull/120)
- Add stop codon substitution conditional branch to protein-variant. [#122](https://github.com/chrovis/varity/pull/122)
- Add variant-type option and protein aliases to find-aliases. [#127](https://github.com/chrovis/varity/pull/127)

### Fixed

- Fix exon/intron boundary determine process and frameshift that affects initiation site. [#92](https://github.com/chrovis/varity/pull/92)
- Fix long deletion including stop codon. [#96](https://github.com/chrovis/varity/pull/96)
- Add extension conditional branch and fix protein-extension. [#99](https://github.com/chrovis/varity/pull/99)
- Fix initiation codon in protein-frameshift and protein-extension. [#100](https://github.com/chrovis/varity/pull/100)
- Fix apply-offset for first exon variant. [#103](https://github.com/chrovis/varity/pull/103)
- Fix read-seq-info and termination site variant condition. [#105](https://github.com/chrovis/varity/pull/105)
- Fix the dispatcher to check if an arg extends the varity.ref_gene.GeneAnnotationIndex protocol. [#106](https://github.com/chrovis/varity/pull/106)
- Fix NullPointerException case in vcf-variant->hgvs. [#108](https://github.com/chrovis/varity/pull/108)
- Fix frameshift processing order. [#110](https://github.com/chrovis/varity/pull/110)
- Fix substitution for nonsense variant. [#113](https://github.com/chrovis/varity/pull/113)
- Explicitly install leiningen on ubuntu-latest. [#114](https://github.com/chrovis/varity/pull/114)
- Fix around ter codon variant. [#116](https://github.com/chrovis/varity/pull/116)
- Fix ini site affected variant determine process. [#125](https://github.com/chrovis/varity/pull/125)
- Fix select-variant coord position. [#126](https://github.com/chrovis/varity/pull/126)
- Revert and fix substitution. [#128](https://github.com/chrovis/varity/pull/128)

## [0.11.0] - 2024-01-22

### BREAKING

We fixed the `varity.vcf-to-hgvs` implementation. It is confusing to throw the
exception in `vcf-variant->protein-hgvs` when a variant overlaps the exon-intron
boundaries, even if coding DNA HGVS is available. So we changed the behavior to
return protein HGVS as `nil`.

### Fixed

- Return protein HGVS as `nil` if a variant overlaps exon-intron boundaries. [#89](https://github.com/chrovis/varity/pull/89)


## [0.10.1] - 2024-01-15

### Fixed

- Add conditional branch of unaffected stop codon to protein-extension. [#86](https://github.com/chrovis/varity/pull/86)

## [0.10.0] - 2023-12-27

### BREAKING

We introduced enhancements to the description of protein changes by
`varity.vcf-to-hgvs`, specifically making deletions more clinically meaningful:

1. exon-intron boundary deletions:

The deletion that overlaps the exon-intron boundary will trigger an Exception
because alterations affecting the splice sites are predicted to be splicing
abnormalities.

2. stop codon deletions:

In cases where deletions contain a stop codon, `varity.vcf-to-hgvs` generates
the following outputs based on the alteration sequence:

- If the alteration sequence contains a stop codon, varity outputs as deletion-insertion.
- Otherwise, this outputs `p.?`.

### Fixed

- Fix upstream and downstream sequence of sequence-info and delins process. [#81](https://github.com/chrovis/varity/pull/81)
- Fix boundary of exon/intron determining process. [#82](https://github.com/chrovis/varity/pull/82)

## [0.9.3] - 2023-07-03

### Fixed

- Fix 3'-rule errors for certain sequence patterns. [#75](https://github.com/chrovis/varity/pull/75)
- Fix overlapping condition of boundary of exon/intron. [#77](https://github.com/chrovis/varity/pull/77)

## [0.9.2] - 2023-02-08

### Fixed

- Fix the 3’ rule on a trailing sub-sequence. [#63](https://github.com/chrovis/varity/pull/63)
- Fix conditional branch of frameshift caused by insertion. [#65](https://github.com/chrovis/varity/pull/65)
- Fix the position of first amino acid changed by the frame shift. [#66](https://github.com/chrovis/varity/pull/66)
- Fix VCF to protein HGVS conversions of insertions near splice sites. [#68](https://github.com/chrovis/varity/pull/68)

## [0.9.1] - 2022-10-11

### Added

- Feature to filter out refgene index source. [#61](https://github.com/chrovis/varity/pull/61)

### Fixed

- Fix frameshift conversion failure. [#59](https://github.com/chrovis/varity/pull/59)

## [0.9.0] - 2022-07-31

### BREAKING

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

### Added

- Annotate fusion genes. [#52](https://github.com/chrovis/varity/pull/52)

### Changed

- Logging liftover failure caused by different refs. [#48](https://github.com/chrovis/varity/pull/48)
- Change the default `:prefer-deletion?` to `false`. [#56](https://github.com/chrovis/varity/pull/56)

### Fixed

- Fix the performance of liftover-variants. [#47](https://github.com/chrovis/varity/pull/47)
- Fix alt exon calculation for deletion. [#50](https://github.com/chrovis/varity/pull/50)
- Tweak GENCODE loading performance. [#55](https://github.com/chrovis/varity/pull/55)

## [0.8.0] - 2021-11-15

### Added

- Compatible with GENCODE's GTF and GFF3 format. [#42](https://github.com/chrovis/varity/pull/42)

### Changed

- Compat with acession number with version. [#43](https://github.com/chrovis/varity/pull/43)
- Conversions from/to HGVS with GENCODE ID. [#45](https://github.com/chrovis/varity/pull/45)

## [0.7.1] - 2021-06-22

### Fixed

- Fix reference protein seq in fs-ter-substitution. [#40](https://github.com/chrovis/varity/pull/40)

## [0.7.0] - 2021-01-12

### Added

- Add support for vcf lift-over. [#35](https://github.com/chrovis/varity/pull/35)

### Changed

- Seek intron, 3'-UTR, and 5'-UTR regions from genomic position. [#36](https://github.com/chrovis/varity/pull/36)
- Update clj-hgvs to 0.4.4.

### Fixed

- Error handling in conversion of an invalid coordinate. [#37](https://github.com/chrovis/varity/pull/37)

## [0.6.2] - 2020-07-06

### Changed

- Update dependencies:
    - clj-hgvs 0.4.3
    - cljam 0.8.0
    - commons-compress 1.20
    - tools.logging 1.1.0

### Fixed

- Pass nil to cds-end if cds-start equals cds-end in vcf-to-hgvs. [#34](https://github.com/chrovis/varity/pull/34)

## [0.6.1] - 2020-02-19

### Added

- Support repeated sequences caused by deletion on VCF->HGVS. [#30](https://github.com/chrovis/varity/pull/30)
- Support conversion from repeated sequences to VCF deletion. [#31](https://github.com/chrovis/varity/pull/31)
- Add prefer-insertion? option to VCF->HGVS. [#32](https://github.com/chrovis/varity/pull/32)
- Add function for finding alternative HGVS expressions. [#33](https://github.com/chrovis/varity/pull/33)

### Fixed

- Fix conversion from reverse-strand repeated sequences to VCF variant. [#29](https://github.com/chrovis/varity/pull/29)

## [0.6.0] - 2019-12-06

### BREAKING

- cdna is renamed to coding-dna to avoid misunderstanding. See [#28](https://github.com/chrovis/varity/pull/28) for more information.
    - namespace:
        - `varity.hgvs-to-vcf.cdna` → `varity.hgvs-to-vcf.coding-dna`
        - `varity.vcf-to-hgvs.cdna` → `varity.vcf-to-hgvs.coding-dna`
    - function:
        - e.g. `cdna-hgvs->vcf-variants` → `coding-dna-hgvs->vcf-variants`
    - return:
        - `{:cdna ...}` → `{:coding-dna ...}`

### Changed

- Rename cdna to coding-dna. [#28](https://github.com/chrovis/varity/pull/28)
- Update dependencies:
    - clj-hgvs 0.4.0
    - cljam 0.7.4
    - commons-compress 1.19
    - tools.logging 0.5.0

## [0.5.1] - 2019-07-12

### Added

- Include MNVs in protein HGVS -> VCF conversion results. [#26](https://github.com/chrovis/varity/pull/26)

### Changed

- Throw ex-info on hgvs->vcf for a variant containing ambiguous coordinates. [#27](https://github.com/chrovis/varity/pull/27)
- Update dependencies:
    - clj-hgvs 0.3.1
    - cljam 0.7.2
    - proton 0.1.8

### Fixed

- Remove incorrect VCF variants from protein HGVS to VCF variants conversion results. [#25](https://github.com/chrovis/varity/pull/25)

## [0.5.0] - 2019-01-29

### BREAKING

clj-hgvs 0.3.0 requires clojure 1.9+ because it uses clojure.spec for HGVS
validation. To use varity with clojure 1.8, you must include a dependency on
[clojure-future-spec](https://github.com/tonsky/clojure-future-spec).

### Changed

- Improve performance of vcf-variant->protein-hgvs. [#22](https://github.com/chrovis/varity/pull/22)
- Update clj-hgvs and proton. [#23](https://github.com/chrovis/varity/pull/23)

### Fixed

- Fix pos of hgvs->vcf-variants for genes on reverse strand. [#24](https://github.com/chrovis/varity/pull/24)

## [0.4.2] - 2018-11-20

### Changed

- Use ex-info for varity specific exception. [#20](https://github.com/chrovis/varity/pull/20)

### Fixed

- Fix normalization of a variant containing long ref (strand: +). [#21](https://github.com/chrovis/varity/pull/21)

## [0.4.1] - 2018-11-04

### Added

- Add debug printing option. [#19](https://github.com/chrovis/varity/pull/19)

### Fixed

- Fix normalization of a variant containing long ref. [#17](https://github.com/chrovis/varity/pull/17)
- Support frame shift with initiation codon change. [#18](https://github.com/chrovis/varity/pull/18)

## [0.4.0] - 2018-10-03

### BREAKING

Strand representation is changed from string (`+`, `-`) to keyword (`:forward`,
`:reverse`).

### Added

- Support promoter on variant conversion. [#10](https://github.com/chrovis/varity/pull/10)
- Profile for Clojure 1.10. [#11](https://github.com/chrovis/varity/pull/11)
- Add several functions to work with refGene exon sequences. [#13](https://github.com/chrovis/varity/pull/13)

### Changed

- Return p.? w/ warning when CDS is indivisible by 3. [#12](https://github.com/chrovis/varity/pull/12)
- Use :forward/:reverse as values for the key ':strand'. [#14](https://github.com/chrovis/varity/pull/14)
- Update dependencies:
    - clj-hgvs 0.2.4
    - cljam 0.6.0
    - commons-compress 1.18

### Fixed

- Fix nonsense substitution in del case. [#9](https://github.com/chrovis/varity/pull/9)
- Fix the translation condition of termination substitution with frameshift. [#15](https://github.com/chrovis/varity/pull/15)
- Fix conversion of no change substitution. [#16](https://github.com/chrovis/varity/pull/16)

## [0.3.7] - 2018-04-27

### Added

- Add benchmark code for liftover. [#8](https://github.com/chrovis/varity/pull/8)

### Changed

- Improve the performance of liftover. [#7](https://github.com/chrovis/varity/pull/7)
- Update dependencies:
    - cavia 0.5.0
    - commons-compress 1.16.1
    - proton 0.1.6

## [0.3.6] - 2018-02-06

### Added

- Add function to retrieve exon regions. [#4](https://github.com/chrovis/varity/pull/4)

### Changed

- Update dependencies:
    - commons-compress 1.16
    - proton 0.1.4

### Fixed

- Fix mistranslation from HGVS to vcf variant. [#5](https://github.com/chrovis/varity/pull/5)

## [0.3.5] - 2017-12-26

### Changed

- Parse scores and exon-frames in refGene.txt. [#3](https://github.com/chrovis/varity/pull/3)
- Update dependencies:
    - cavia 0.4.3
    - clj-hgvs 0.2.3
    - cljam 0.5.1
    - proton 0.1.3

### Fixed

- Fix incorrect conversion of nonsense substitution.

## [0.3.4] - 2017-08-21

### Changed

- Update cljam (0.4.1) and proton (0.1.2) version.

### Fixed

- Fix errors for some uncommon variants. [#2](https://github.com/chrovis/varity/pull/2)

## [0.3.3] - 2017-06-30

### Added

- Add tests.
- Support TwoBit format as reference sequence file.

### Fixed

- Fix the process of calculating alt protein sequence.

## [0.3.2] - 2017-06-15

### Changed

- Bump clj-hgvs version up to 0.2.0.

## [0.3.1] - 2017-05-15

### Changed

- Improve cDNA HGVS -> VCF variant support.
- ref-seq info can be used for HGVS -> VCF.
- Bump dependencies version up.

### Fixed

- Fix missed arguments.
- Fix arguments of lift in README.md. [#1](https://github.com/chrovis/varity/pull/1)
- Fix reflection/boxing warnings.

## [0.3.0] - 2017-05-08

### Added

- Normalize chromosome of VCF->HGVS inputs.

### Changed

- Use :chr instead of :chrom as chromosome of ref-gene.
- Use record as ref-gene index instead of map.
- Modify arguments of VCF<->HGVS.
- Modify arguments of lift.
- Bump clj-hgvs version up to 0.1.2.

## 0.2.0 - 2017-04-20

First release.

[Unreleased]: https://github.com/chrovis/varity/compare/0.12.0...HEAD
[0.12.0]: https://github.com/chrovis/varity/compare/0.11.0...0.12.0
[0.11.0]: https://github.com/chrovis/varity/compare/0.10.1...0.11.0
[0.10.1]: https://github.com/chrovis/varity/compare/0.10.0...0.10.1
[0.10.0]: https://github.com/chrovis/varity/compare/0.9.3...0.10.0
[0.9.3]: https://github.com/chrovis/varity/compare/0.9.2...0.9.3
[0.9.2]: https://github.com/chrovis/varity/compare/0.9.1...0.9.2
[0.9.1]: https://github.com/chrovis/varity/compare/0.9.0...0.9.1
[0.9.0]: https://github.com/chrovis/varity/compare/0.8.0...0.9.0
[0.8.0]: https://github.com/chrovis/varity/compare/0.7.1...0.8.0
[0.7.1]: https://github.com/chrovis/varity/compare/0.7.0...0.7.1
[0.7.0]: https://github.com/chrovis/varity/compare/0.6.2...0.7.0
[0.6.2]: https://github.com/chrovis/varity/compare/0.6.1...0.6.2
[0.6.1]: https://github.com/chrovis/varity/compare/0.6.0...0.6.1
[0.6.0]: https://github.com/chrovis/varity/compare/0.5.1...0.6.0
[0.5.1]: https://github.com/chrovis/varity/compare/0.5.0...0.5.1
[0.5.0]: https://github.com/chrovis/varity/compare/0.4.2...0.5.0
[0.4.2]: https://github.com/chrovis/varity/compare/0.4.1...0.4.2
[0.4.1]: https://github.com/chrovis/varity/compare/0.4.0...0.4.1
[0.4.0]: https://github.com/chrovis/varity/compare/0.3.7...0.4.0
[0.3.7]: https://github.com/chrovis/varity/compare/0.3.6...0.3.7
[0.3.6]: https://github.com/chrovis/varity/compare/0.3.5...0.3.6
[0.3.5]: https://github.com/chrovis/varity/compare/0.3.4...0.3.5
[0.3.4]: https://github.com/chrovis/varity/compare/0.3.3...0.3.4
[0.3.3]: https://github.com/chrovis/varity/compare/0.3.2...0.3.3
[0.3.2]: https://github.com/chrovis/varity/compare/0.3.1...0.3.2
[0.3.1]: https://github.com/chrovis/varity/compare/0.3.0...0.3.1
[0.3.0]: https://github.com/chrovis/varity/compare/0.2.0...0.3.0
