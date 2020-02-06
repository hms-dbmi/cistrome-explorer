/**
 * Construct our own d3 object, with only the functions that we are using.
 * This should reduce the size of the final javascript bundle file.
 * See https://github.com/d3/d3/issues/3076
 */ 

import { select } from "d3-selection";
import { format } from "d3-format";
import { schemeSet3 } from "d3-scale-chromatic";
import { mouse, event as d3_event } from "d3-selection";
import { scaleOrdinal, scaleThreshold } from "d3-scale";
import { hsl } from "d3-color";

import { scale as vega_scale } from "vega-scale";

/*
 * Same as `d3.scaleBand` but also supports `invert()`.
 * References: 
 * - https://github.com/d3/d3-scale/pull/64#issuecomment-319516140
 * - https://github.com/vega/vega/blob/5d76c88eeead37e4589f9a4278ea1296683f5f68/packages/vega-scale/src/scales/scaleBand.js
 */
const scaleBand = vega_scale("band");


export default {
    select,
    format,
    schemeSet3,
    mouse,
    get event() { return d3_event; }, // https://stackoverflow.com/a/40048292
    scaleOrdinal,
    scaleThreshold,
    scaleBand,
    hsl
};