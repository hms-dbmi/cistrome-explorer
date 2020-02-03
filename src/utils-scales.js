import { scale as vega_scale } from 'vega-scale';

/*
 * Same as `d3.scaleBand` but also supports `invert()`.
 * References: 
 * - https://github.com/d3/d3-scale/pull/64#issuecomment-319516140
 * - https://github.com/vega/vega/blob/5d76c88eeead37e4589f9a4278ea1296683f5f68/packages/vega-scale/src/scales/scaleBand.js
 */
export const vega_scaleBand = vega_scale('band');