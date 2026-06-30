# nf-core/mcmicro: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - [06/30/2026]

Initial release of nf-core/mcmicro, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Prelude option to check and validate image metadata and report inconsistencies/errors.
- Cellpose is now default segmentation tool.

### `Fixed`

- No longer uses custom directory structure as input.
- Input is no longer restricted to be in same directory as output.
- Updated illumination to basicpy.
- No longer need custom params.yml. Can pass arguments to CLI or nextflow params-file.

### `Dependencies`

### `Deprecated`

- Unmicst
- Illumination tool.
