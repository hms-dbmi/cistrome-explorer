## 0.4.0 - in progress

### Added

### Changed
- Changed the table view to be shown as a modal view.
- Updated format of API results (e.g., use consistent letter cases and use scientific notations for long numbers).

## 0.3.0 - 05/27/20

### Added
- Added real data (quality score for each sample) for quantitative bar charts.
- Added a mover component that allows one to interactively reposition the keyword search box.
- Added filter reset buttons on the wrapper tracks to make it more accessible.
- Added `range` condition in filter and highlighting options to enable filtering rows by value ranges.
- Added a range slider to interactively filter rows by determining min and max cutoff values in quantitative bar charts, including stacked bar charts.
- Added a vertical mouse wheel zooming feature for horizontal multivec tracks, activated when the `Y` key is pressed on the keyboard.
- Added a minimum similarity bar in dendrogram to enable filtering rows by similarity scores.
- Added a simple nav bar with links to GitHub and documentations.

### Changed
- Make keyword search box semi-transparent when blurred.
- Changed the color of filter and sort icons when they are applied on a track as a visual indicator.
- Changed the filtering interface for nominal values to use checkboxes.
- Render range sliders and checkboxes for filtering interactions based on the original data.
- Replace `contains` filtering condition to `notOneOf`.
- Rename the filtering condition for `tree` type data from `contains` to `subtree`.
- Reset buttons remove all filters related only to a certain field.
- Each field has only one `filterInfo` for the simplicity.
- Update dendrogram to visually encode similarity distance between nodes.

## 0.2.0 - 03/18/20

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
- Added a context menu for filtering or highlighting rows based on a nomival value in categorical bar charts.
- Added hover styles for the `TrackRowInfoVisLink` visualization component, so that when a link area is hovered, a gray background rect and a text underline are also rendered.
- Added a context menu to add a new track on top or bottom with selected rows in categorical bar charts.
- Added `viewWrapper` react component that support view-specific interactions, such as brusing columns of higlas heatmaps.
- Added an option to render a checkbox for each row of the data table, with a callback function that takes an array of checked row objects.
- Added wrapper visualizations to SVGs exported by HiGlass using the new `.on('createSVG')` callback feature.

### Changed
- Fixed bug which prevented wrapper options for sorting, filtering, highlighting from being used upon initial component render.
- Added resizer elements to their containers and used `d3-drag` for the mouse events.
- Compute color scale for the `TrackRowInfoVisNominalBar` categorical value visualization based on the full `rowInfo` rather than `transformedRowInfo` so that the original color mapping is kept after row filtering events.
- Fixed the genomic interval selection behavior to reflect the support for `fromViewUid: null` for `viewport-projection-horizontal` tracks (new in HiGlass 1.9.0).
- Set the `horizontal-multivec` track `zeroValueColor` option to demonstrate the feature (new in HiGlass 1.9.0).
- Enabled multiple genomic interval selections for a view.


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