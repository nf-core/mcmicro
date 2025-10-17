# nf-core/mcmicro: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Directory Structure](#directory-structure)
- [Summary](#summary)
- [Illumination Correction](#illumination-correction)
  - [BaSiCPy](#basicpy)
- [Registration](#registration)
  - [ASHLAR](#ashlar)
- [Background Subtraction](#background-subtraction)
  - [Backsub](#backsub)
- [TMA Core Separation](#tma-core-separation)
  - [Coreograph](#coreograph)
- [Segmentation](#segmentation)
  - [Mccellpose](#mccellpose)
  - [Cellpose](#cellpose)
  - [Mesmer](#mesmer)
- [Quantification](#quantification)
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Directory Structure

```
{outdir}
├── summary
├── backsub
├── illumination_correction
│   └── basicpy
├── multiqc
│   ├── multiqc_data
│   ├── multiqc_plots
│   └── multiqc_report.html
├── pipeline_info
├── quantification
│   └── mcquant
│       └── {segmentation module}
├── registration
│   └── ashlar
├── segmentation
│   └── {segmentation module}
└── tma_dearray
    └── masks

```

### Summary

We generate a MultiQC formatted report with the extracted OME-xml data and the marker sheet information so the user can debug their parameters.

<details>
<summary>Output files</summary>

- summary\_{sample_name}\_samplesheet_mqc.tsv : MultiQC formatted TSV with sample's parameter validation
- summary\_{sample_name}\_xml_mqc.tsv : MultiQC formatted TSV with sample's extracted OME-xml data validation
- summary_markersheet_mqc.tsv : Formatted Markersheet in MultiQC format

</details>

### Illumination Correction

#### BaSiCPy

[BaSiCPy](https://nf-co.re/modules/basicpy/) is a python package for background and shading correction of optical microscopy images. It is developed based on the Matlab version of BaSiC tool with major improvements in the algorithm.

<details>
<summary>Output files</summary>

- {sample_name}-dfp.tif : Tiff fields for dark field illumination correction
- {sample_name}-ffp.tif : Tiff fields for flat field illumination correction

</details>

### Registration

#### ASHLAR

[ASHLAR](https://nf-co.re/modules/ashlar/) combines multi-tile microscopy images into a high-dimensional mosaic image.

<details>
<summary>Output files</summary>

- {sample_name}.ome.tif : A pyramidal, tiled OME-TIFF file created from input images.

</details>

### Background Subtraction

#### Backsub

[Backsub](https://nf-co.re/modules/backsub/) performs a pixel-by-pixel channel subtraction scaled by exposure times of pre-stitched tif images.

<details>
<summary>Output files</summary>

- markers_bs.csv : Marker file adjusted to match the background corrected image
- .backsub.ome.tif : Background corrected pyramidal ome.tif

</details>

### TMA Core Separation

#### Coreograph

[Coreograph](https://nf-co.re/modules/coreograph/) uses UNet, a deep learning model, to identify complete/incomplete tissue cores on a tissue microarray. It has been trained on 9 TMA slides of different sizes and tissue types.

<details>
<summary>Output files</summary>

- {core_number}.tif : Individual cropped tissue core images
- centroidsY-X.txt : A text file listing centroids of each core in format Y, X
- masks/{core_number}\_mask.tif : Binary mask image for each tissue core
- TMA_MAP.tif : A TMA map showing core number labels and mask outlines

</details>

### Segmentation

#### Mccellpose

[Mccellpose] A RAM-efficient wrapper around Cellpose (see below).

<details>
<summary>Output files</summary>

- {sample_name}\_mask.ome.tif : labelled mask output from cellpose in OME-TIFF format

</details>

#### Cellpose

[Cellpose](https://nf-co.re/modules/cellpose/) segments cells in images

<details>
<summary>Output files</summary>

- {sample_name}\_mask.tif : labelled mask output from cellpose in tif format

</details>

#### Mesmer

[Mesmer](https://nf-co.re/modules/deepcell_mesmer/) segmentation for whole-cell

<details>
<summary>Output files</summary>

- {sample_name}\_mask.tif : File containing the mask.

</details>

### Quantification

#### Mcquant

[Mcquant](https://nf-co.re/modules/mcquant/) extracts single-cell data given a multi-channel image and a segmentation mask.

<details>
<summary>Output files</summary>

- {segmenter}/{sample_name}.csv : Single-cell feature table for all selected segmenters.

</details>

### Quality Control

### MultiQC

Aggregate report describing results and QC from the whole pipeline

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

Report metrics generated during the workflow execution

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
