# cbGemma

**cbGemma** is a Nextflow pipeline for generating [UCSC Cell Browser](https://cellbrowser.readthedocs.io/) directories from Pavlidis Lab [GEMMA](https://gemma.msl.ubc.ca/) single-cell datasets. It downloads expression matrices and cell-level metadata for one or more Gemma studies, processes them with `cbScanpy`, and builds interactive Cell Browser web pages using `cbBuild`.

---

## Usage


```bash
git clone https://github.com/<your-username>/cbGemma.git
cd cbGemma
nextflow main.nf --study_names <experiment name 1>  <experiment name 2>
```
Default profile is set to `conda` with a pre-existing environment path. Use `-resume` to resume.

## Outputs
Cellbrowser builds are be published in `/space/gemmaData/cellBrowser` by default per the `params.cb_outdir` parameter.
The resulting web page will be viewable at:

https://dev.gemma.msl.ubc.ca/cellbrowser/experiment/?ds=cellBrowser
