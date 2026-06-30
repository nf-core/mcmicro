# nf-core/mcmicro: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/mcmicro/usage](https://nf-co.re/mcmicro/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. We currently accept 2 formats for the input samplesheets. One format is one row per sample and the other is one row per sample per cycle. Use the parameter `input_sample` for one row per sample or the parameter `input_cycle` for one row per sample per cycle, to specify its location. It has to be a comma-separated file with a header row and either two (input_sample) or four (input_cycle) columns as shown in the examples below.

```bash
--input_cycle '[path to one row per sample per cycle samplesheet file]'
```

**OR**

```bash
--input_sample '[path to one row per sample samplesheet file]'
```

### Samplesheet with one row per sample per cycle

The `sample` identifier must be the same for multiple cycles of the same sample. All the files from the same sample will be run in a single run of ashlar in the cycle order that they appear in the samplesheet. If illumination correction is requested using basicpy, each cycle will be corrected separately.

```csv title="samplesheet_cycle.csv"
sample,cycle_number,image_tiles
TEST1,1,/path/to/image/cycif-tonsil-cycle1.ome.tif
TEST1,2,/path/to/image/cycif-tonsil-cycle2.ome.tif
TEST1,3,/path/to/image/cycif-tonsil-cycle3.ome.tif
```

| Column         | Description                                                 |
| -------------- | ----------------------------------------------------------- |
| `sample`       | Custom sample name.                                         |
| `cycle_number` | Integer value of the cycle for the file in the current row. |
| `image_tiles`  | Full path or URL to the input image file.                   |

An [example one row per sample per cycle samplesheet](../assets/samplesheet_1_row_sample_cycle.csv) has been provided with the pipeline.

### Samplesheet with one row per sample

This is similar to the above case except each row just contains a column for each `sample` name and a columnn containing a directory where all the files for a given sample are located. All per-cycle image files in the `image_directory` for a given sample will be run in a single run of ashlar. If illumination correction is requested using basicpy, each cycle will be corrected separately.

```csv title="samplesheet_sample.csv"
sample,image_directory
TEST1,/path/to/image/directory
```

| Column            | Description                                          |
| ----------------- | ---------------------------------------------------- |
| `sample`          | Custom sample name.                                  |
| `image_directory` | Full path to directory containing input image files. |

An [example one row per sample samplesheet](../assets/samplesheet_1_row_sample.csv) has been provided with the pipeline.

## Markersheet input

Each row of the markersheet represents a single channel in the associated sample image. The columns `channel_number`, `cycle_number` and `marker_name` are required.

```csv
channel_number,cycle_number,marker_name
1,1,DNA 1
2,1,Na/K ATPase
3,1,CD3
4,1,CD45RO
```

| Column           | Description                                         |
| ---------------- | --------------------------------------------------- |
| `channel_number` | Integer identifier for the respective channel.      |
| `cycle_number`   | Integer identifier for the image cycle.             |
| `marker_name`    | Name of the marker for the given channel and cycle. |

:::note
`cycle_number` must match the `cycle_number` in the supplied samplesheet.
:::

### optional markersheet columns

| Column                  | Description                                    |
| ----------------------- | ---------------------------------------------- |
| `filter`                | Microscope filter common name.                 |
| `excitation_wavelength` | Excitation wavelength for this channel, in nm. |
| `emission_wavelength`   | Emission wavelength for this channel, in nm.   |

## Running the pipeline

### One row per sample per cycle

The typical command for running the one row per sample per cycle pipeline is as follows:

```bash
nextflow run nf-core/mcmicro --input_cycle ./samplesheet_cycle.csv --outdir ./results --marker_sheet markers.csv -profile docker
```

### One row per sample

The typical command for running the one row per sample pipeline is as follows:

```bash
nextflow run nf-core/mcmicro --input_sample ./samplesheet_sample.csv --outdir ./results --marker_sheet markers.csv -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/mcmicro -profile docker -params-file params.yaml
```

with:

```yaml
input_cycle: "samplesheet_cycle.csv"
outdir: "./output"
marker_sheet: "markers.csv"
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Pipeline stages and associated input parameters

#### Summary

We generate a MultiQC formatted report based on the extracted OME-xml metadata and the marker/sample sheet information so the user can debug their runs.
This step is always run, but we include the option to stop immediately after the data extraction/validation with the option `--prelude`.

Running with `--prelude` is recommended as a first run to ensure your sample/marker sheets are complete and the required metadata is present in your image files.

#### Illumination Correction

Illumination correction can optionally be performed before registration. It is triggered by the `--illumination` flag which can currently only be followed by the option `basicpy`. We plan on supporting other modules for illumination correction in the future.
When `basicpy` is selected the nf-core module basicpy is run on the input image(s). Basicpy is a python package for background and shading correction of optical microscopy images. More information about it can be found on the [basicpy nf-core module website](https://nf-co.re/modules/basicpy/).

#### Registration

Registration is a required step of the pipeline and the only module currently supported is ashlar. Ashlar is a software package that combines multi-tile microscopy images into a high-dimensional mosaic image. More information about ashlar can be found on the [ashlar website](https://labsyspharm.github.io/ashlar/). We plan to support other modules for registration in the future.

#### Background Subtraction

This is an optional step that occurs immediately following registration. It is triggered by the `--backsub` flag. When this flag is selected, the module backsub is run on the output from the registration step. The backsub module performs pixel-by-pixel channel subtraction scaled by exposure times of pre-stitched tif images. More information about it can be found on the [backsub nf-core module website](https://nf-co.re/modules/backsub/).

#### TMA Core Separation

This is an optional step that occurs immediately following background subtration if that optional step was run or after registration if is was not. It is triggered by the `--tma_dearray` flag. When this flag is selected, the coreograph module is run on the output from either the background subtraction step or the registration step if background subtration was not performed. Coreograph separates the input image into a set of images for each of the cores. It uses UNet, a deep learning model, to identify complete/incomplete tissue cores on a tissue microarray. It has been trained on 9 TMA slides of different sizes and tissue types. More information about it can be found on the [coreograph nf-core module website](https://nf-co.re/modules/coreograph/)

#### Segmentation

This is a required step that follows Registration (or the optional Background Subtraction and TMA Core Separation steps if enabled). The workflow will run the mccellpose module by default, but other options are available by using the `--segmentation` flag. The flag should be followed by a single segmentation module name or a comma separated list of names to run multiple segmentation modules in parallel. The available options currently supported are `mccellpose`, `cellpose`, and `mesmer`. More information about each of these modules can be found on their respective nf-core module websites: [cellpose](https://nf-co.re/modules/cellpose/) [deepcell_mesmer](https://nf-co.re/modules/deepcell_mesmer/). (mccellpose is an alternative interface to cellpose that uses much less RAM and should be preferred in most cases)

When `cellpose` is selected as a segmentation method you may also provide a pretrained model to the cellpose module by using the `--cellpose_model` flag followed by a full path or URL to the model file.

#### Quantification

This is a required step that follows segmentation. The workflow currently runs the mcquant module by default. Other quantification modules will be added as options in the future. Mcquant extracts single-cell data given a multi-channel image and a segmentation mask. More information about mcquant can be found on the [mcquant nf-core module website](https://nf-co.re/modules/mcquant/).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/mcmicro
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/mcmicro releases page](https://github.com/nf-core/mcmicro/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
