## 0.2.0 - in progress

### Added
- Added a feature to interactively resize the width of each vertical track by dragging resizers.
- Added a feature to interactively filter rows by keywords and reset all filters.
- Added an option to specify multiple data fields in a single vertical track and implemented stacked bar charts for this case.
- Added a feature to interactively filter rows by selecting a node in the dendrogram visualization.
- Added x-axis ticks for quantitative bar plot visualizations for row info metadata.
- Added transparent white background for visualization title text so that it does not visually interfere with the text for plot data labels.
- Added `measureText()` two function for computing width and height of `TwoText` objects without rendering.
- Added a search suggestion autocomplete feature.
- Added a feature to request for bind TFs in a certain interval using Cistrome DB API and render the result in a simple data table.

### Changed
- Fixed bug which prevented wrapper options for sorting, filtering, highlighting from being used upon initial component render.


## 0.1.0 - 02/19/20

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
- Added a global variable that stores `rowInfo` as well as sorting & highlighting info, using `React.Context`.
- Added a mechanism to sort rows of a HiGrass `horizontal-multivec` track and vertical tracks on the left and right sides which is activated when users click on sorting buttons.
- Added a feature to interactively highlight rows by typing keywords on a text field which appears when a search button is clicked.

### Changed
- Moved to rendering the row info categorical metadata on a single `<canvas/>` element, rather than colored HTML elements for each row.
- Updated all `d3` imports so that they refer to a specific sub-package (such as `d3-selection` or `d3-scale`).
- Changed rollup config `output.sourcemap` values to `true` rather than `"inline"`, which means that sourcemaps are stored in separate files, reducing the size of the `.js` files.
- Updated the positioning of the genomic interval start and end coordinate elements, so that they are dynamically positioned based on the current genomic region in the track.
- Added a small implementation of renderer-agnostic drawing code, based on the `two.js` JavaScript library, which enables the same code to draw to SVG and canvas elements. For example, when rendering in the browser, canvas is preferred for speed, but when rendering for export to a file, SVG is preferred for resolution.