// Functions here have been copied directly from HiGlass.
// Reference: https://github.com/higlass/higlass/blob/develop/app/scripts/utils/track-utils.js


/**
 * Calculate the tiles that should be visible get a data domain
 * and a tileset info
 *
 * All the parameters except the first should be present in the
 * tileset_info returned by the server.
 * 
 * Reference: https://github.com/higlass/higlass/blob/cdca9990a598c99529cfb6d722c292368f426923/app/scripts/utils/track-utils.js#L148
 *
 * @param zoomLevel: The zoom level at which to find the tiles (can be
 *   calculated using this.calculateZoomLevel, but needs to synchronized across
 *   both x and y scales so should be calculated externally)
 * @param scale: A d3 scale mapping data domain to visible values
 * @param minX: The minimum possible value in the dataset
 * @param maxX: The maximum possible value in the dataset
 * @param maxZoom: The maximum zoom value in this dataset
 * @param maxDim: The largest dimension of the tileset (e.g., width or height)
 *   (roughlty equal to 2 ** maxZoom * tileSize * tileResolution)
 */
const calculateTiles = (zoomLevel, scale, minX, maxX, maxZoom, maxDim) => {
    const zoomLevelFinal = Math.min(zoomLevel, maxZoom);
  
    // the ski areas are positioned according to their
    // cumulative widths, which means the tiles need to also
    // be calculated according to cumulative width
  
    const tileWidth = maxDim / 2 ** zoomLevelFinal;
    const epsilon = 0.0000001;
  
    return range(
      Math.max(0, Math.floor((scale.domain()[0] - minX) / tileWidth)),
      Math.min(
        2 ** zoomLevelFinal,
        Math.ceil((scale.domain()[1] - minX - epsilon) / tileWidth)
      )
    );
};

/**
 * Get the tile's position in its coordinate system.
 *
 * See Tiled1DPIXITrack.js
 * 
 * Reference: https://github.com/higlass/higlass/blob/cdca9990a598c99529cfb6d722c292368f426923/app/scripts/utils/track-utils.js#L405
 */
export const getTilePosAndDimensions = (tilesetInfo, tileId) => {
    const zoomLevel = +tileId.split('.')[0];
    const xTilePos = +tileId.split('.')[1];
  
    // max_width should be substitutable with 2 ** tilesetInfo.max_zoom
    const totalWidth = tilesetInfo.max_width;
  
    const minX = tilesetInfo.min_pos[0];
  
    const tileWidth = totalWidth / 2 ** zoomLevel;
  
    const tileX = minX + xTilePos * tileWidth;
  
    return {
      tileX,
      tileWidth
    };
};