# Changelog

## [Unreleased]

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

[Unreleased]: https://github.com/chrovis/varity/compare/0.5.1...HEAD
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
