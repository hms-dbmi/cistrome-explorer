

## 0.1.0 - in progress

### Added
- Added radio buttons for selecting between two higlass viewconfig options.
- Added higlass viewconfig corresponding to local `higlass-server` demo tileset containing JSON metadata.
- Configured [documentation](https://github.com/documentationjs/documentation).
- Added rollup configuration for packaging the component for NPM.
- Display metadata values in tooltip upon mouse hover interactions for `TrackRowInfo` component.
- Added `TrackRowLink` component, and implemented dynamic dimensions based on the `rowLinkPosition` option and dynamic show/hide based on row height. URLs are obtained based on the `rowLinkAttribute` option which defaults to `"url"`.
- Added `TrackColTools` component, shows a button for creating a genome interval selection, and can be positioned above or below `horizontal-multivec` tracks using the `colToolsPosition` option.
- Added a dendrogram visualization to the `TrackRowInfo` component, which is displayed when `type === "tree"` for a row info attribute in the `options` prop.
- Added a quantitative attribute visualization (vertical bar plot) to the `TrackRowInfo` component, which is displayed when `type === "quantitative"` for a row info attribute in the `options` prop.
- Added buttons for sorting according to a particular attribute, which appear upon hovering over the visualization for the attribute in the `TrackRowInfo` component.
- Added a new demo with a HiGlass viewconfig containing two `horizontal-multivec` tracks, where two sets of options for the wrapper component are specified based on the `viewId` and `trackId` of each track.
- Added a mechanism for components to "register" their `draw()` functions, to prepare for drawing visualizations to SVG when the user would like to export the current visual state.

### Changed
- Moved to rendering the row info categorical metadata on a single `<canvas/>` element, rather than colored HTML elements for each row.
- Updated all `d3` imports so that they refer to a specific sub-package (such as `d3-selection` or `d3-scale`).
- Changed rollup config `output.sourcemap` values to `true` rather than `"inline"`, which means that sourcemaps are stored in separate files, reducing the size of the `.js` files.
- Updated the positioning of the genomic interval start and end coordinate elements, so that they are dynamically positioned based on the current genomic region in the track.
- Added a small implementation of renderer-agnostic drawing code, based on the `two.js` JavaScript library, which enables the same code to draw to SVG and canvas elements. For example, when rendering in the browser, canvas is preferred for speed, but when rendering for export to a file, SVG is preferred for resolution.