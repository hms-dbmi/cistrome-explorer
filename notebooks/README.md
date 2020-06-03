# notebooks

One-off jupyter notebooks used for small data processing tasks can be placed here.

Consider extracting tasks that grow more intense / frequently-used, converting to more robust  [`higlass-server`](https://github.com/higlass/higlass-server/) or [`clodius`](https://github.com/higlass/clodius) processing steps.

## Setup

```sh
conda env create -f environment.yml
```

## Start

```sh
conda activate cistrome-explorer-notebooks
jupyter notebook
```