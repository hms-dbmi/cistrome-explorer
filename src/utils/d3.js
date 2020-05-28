/*
 * Construct our own d3 object, with only the functions that we are using.
 * This should reduce the size of the final javascript bundle file.
 * See https://github.com/d3/d3/issues/3076
 */ 

import { select, create } from "d3-selection";
import { format } from "d3-format";
import { schemeSet3, schemeTableau10, interpolateViridis } from "d3-scale-chromatic";
import { mouse, event as d3_event } from "d3-selection";
import { scaleLinear, scaleOrdinal, scaleThreshold } from "d3-scale";
import { axisTop, axisBottom } from "d3-axis";
import { scale as vega_scale } from "vega-scale";
import { extent, sum, mean, max, min } from "d3-array";
import { hsl } from "d3-color";
import { hierarchy, cluster } from "d3-hierarchy";
import { Delaunay as delaunay } from "d3-delaunay";
import { drag } from "d3-drag";

/*
 * Same as `d3.scaleBand` but also supports `invert()`.
 * References: 
 * - https://github.com/d3/d3-scale/pull/64#issuecomment-319516140
 * - https://github.com/vega/vega/blob/5d76c88eeead37e4589f9a4278ea1296683f5f68/packages/vega-scale/src/scales/scaleBand.js
 */
const scaleBand = vega_scale("band");

export default {
    get event() { return d3_event; }, // https://stackoverflow.com/a/40048292
    select,
    create,
    format,
    schemeSet3,
    schemeTableau10,
    interpolateViridis,
    mouse,
    scaleLinear,
    scaleOrdinal,
    scaleThreshold,
    axisTop,
    axisBottom,
    scaleBand,
    extent,
    sum,
    mean,
    max, 
    min,
    hsl,
    hierarchy,
    cluster,
    delaunay,
    drag
};