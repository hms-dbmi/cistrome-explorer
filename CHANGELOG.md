

## 0.1.0 - in progress

### Added
- Added radio buttons for selecting between two higlass viewconfig options.
- Added higlass viewconfig corresponding to local `higlass-server` demo tileset containing JSON metadata.
- Configured [documentation](https://github.com/documentationjs/documentation).
- Added rollup configuration for packaging the component for NPM.
- Display metadata values in tooltip upon mouse hover interactions for `TrackRowInfo` component.
- Added `TrackRowLink` component files, and implemented dynamic dimensions based on the `rowLinkPosition` option.

### Changed
- Moved to rendering the row info categorical metadata on a single `<canvas/>` element, rather than colored HTML elements for each row.
- Updated all `d3` imports so that they refer to a specific sub-package (such as `d3-selection` or `d3-scale`).
- Changed rollup config `output.sourcemap` values to `true` rather than `"inline"`, which means that sourcemaps are stored in separate files, reducing the size of the `.js` files.