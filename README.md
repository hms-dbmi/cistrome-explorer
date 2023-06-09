# Cistrome Explorer

Interactive visual analytic tool for exploring epigenomics data w/ associated metadata, powered by [HiGlass](http://higlass.io/) and [Cistrome Data Browser Toolkit](http://dbtoolkit.cistrome.org/)

-   [Demo](http://cisvis.gehlenborglab.org/)
-   [Documentation](http://cisvis.gehlenborglab.org/docs/)

### Development

Install dependencies with yarn ([v1](http://classic.yarnpkg.com)):

```sh
yarn
```

Run development server for demo app:

```sh
yarn start
```

Run tests:

```sh
yarn test
# or
yarn test --watch
```

### Production

Build demo app for production:

```sh
yarn build
```

Build for NPM package:

```sh
yarn build-pkg
```

### Docs

```sh
yarn docs
```

### Data processing resources

-   Notebooks for short analyses or file conversions: [./notebooks/](./notebooks/)
-   Pipeline for combining CistromeDB bigWig files into HiGlass multivec (HDF5-based) files: [./pipelines/cistrome-to-multivec/](./pipelines/cistrome-to-multivec/)

### Presentations and other resources

-   [Video](https://drive.google.com/file/d/1SrtFHrEuJY5zHuPjPkBmPTxgZPRQ0qRR/view) & [Slides](https://drive.google.com/file/d/1Z4tO-lrClZY3P7_n2N3kar5YoQoMNVCh/view?usp=sharing) for [NCI ITCR 2020 virtual poster](https://ncihub.org/groups/itcr/2020_virtual_posters)
-   [Poster](https://drive.google.com/file/d/1r0jPwyTlEYGotsrfD2KbJU5r-OEYU5Q5/view?usp=sharing) for [BioVis@ISMB 2020](http://biovis.net/2020/program_ismb/)

### Related repositories

-   [HiGlass](https://github.com/higlass/higlass)
-   [HiGlass Server](https://github.com/higlass/higlass-server) and our fork [Cistrome Explorer HiGlass Server](https://github.com/hms-dbmi/cistrome-explorer-higlass-server)
    -   To deploy `cistrome-explorer-higlass-server` to AWS ECS please refer to the modified `docker-context/` directory and associated README [here](https://github.com/hms-dbmi/cistrome-explorer-higlass-server/blob/develop/docker-context/README.md).
-   [clodius](https://github.com/higlass/clodius)
-   [pybbi](https://github.com/nvictus/pybbi)
