# cbGemma

**cbGemma** is a Nextflow pipeline for generating [UCSC Cell Browser](https://cellbrowser.readthedocs.io/) directories from Pavlidis Lab [GEMMA](https://gemma.msl.ubc.ca/) single-cell datasets. It downloads expression matrices and cell-level metadata for one or more Gemma studies, processes them with `cbScanpy`, and builds interactive Cell Browser web pages using `cbBuild`.

---

## Usage


```bash
git clone https://github.com/<your-username>/cbGemma.git
cd cbGemma
nextflow main.nf --study_names <experiment1> <experiment2>
# optionally add the -resume flag to resume from nextflow cache
```
Default profile is set to `conda` with a pre-existing environment path. Use `-resume` to resume.

## Outputs
Cellbrowser builds are be published in `/space/gemmaData/cellBrowser` by default per the `params.cb_outdir` parameter.
The resulting web page will be viewable on the development server. For example:

https://dev.gemma.msl.ubc.ca/cellbrowser/Velmeshev_et_al.2/?ds=Velmeshev_et_al_2
