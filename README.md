# cistrome-explorer

_wrapper component for [HiGlass](http://higlass.io/), with features for viewing [CistromeDB](http://cistrome.org/db/) data and querying [CistromeDB Toolkit](http://dbtoolkit.cistrome.org/)_

- [Demo](http://cisvis.gehlenborglab.org/)
- [Documentation](http://cisvis.gehlenborglab.org/docs/)

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
yarn build-demo
```

Build for NPM package:

```sh
yarn build-pkg
```

The `higlass-meta` package is published to the NPM registry by Travis when the version in `package.json` has been updated and pushed to the `master` branch. To perform this update:
- Check out a new branch for the release,
    - Update the CHANGELOG.md to remove the "in progress" text from the current version heading.
    - Update the version by running `npm version [major | minor | patch]` (note: this will add a git commit and a git tag).
- Make a pull request to merge from the release branch into `master`.

Travis uses the `NPM_EMAIL` and `NPM_TOKEN` variables that can be set using the [web interface](https://travis-ci.org/github/hms-dbmi/cistrome-explorer/settings) (Settings -> Environment Variables).

### Docs

```sh
yarn docs
```

### Creating and loading multivec files with JSON metadata

- [Example](https://github.com/keller-mark/clodius-cistrome-example)

### Other resources
- [Video](https://drive.google.com/file/d/1SrtFHrEuJY5zHuPjPkBmPTxgZPRQ0qRR/view) & [Slides](https://drive.google.com/file/d/1Z4tO-lrClZY3P7_n2N3kar5YoQoMNVCh/view?usp=sharing) for [2020 ICTR virtual poster](https://ncihub.org/groups/itcr/2020_virtual_posters)

### Related repositories

- [HiGlass](https://github.com/higlass/higlass)
- [HiGlass Server](https://github.com/higlass/higlass-server) and our fork [Cistrome Explorer HiGlass Server](https://github.com/hms-dbmi/cistrome-explorer-higlass-server)
