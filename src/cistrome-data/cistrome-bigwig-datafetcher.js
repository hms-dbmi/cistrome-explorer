import { BigWig } from "@gmod/bbi";
import { RemoteFile } from "generic-filehandle";

export const CHROM_SIZE_HG38 = {
    chr1: 248956422,
    chr2: 242193529,
    chr3: 198295559,
    chr4: 190214555,
    chr5: 181538259,
    chr6: 170805979,
    chr7: 159345973,
    chr8: 145138636,
    chr9: 138394717,
    chr10: 133797422,
    chr11: 135086622,
    chr12: 133275309,
    chr13: 114364328,
    chr14: 107043718,
    chr15: 101991189,
    chr16: 90338345,
    chr17: 83257441,
    chr18: 80373285,
    chr19: 58617616,
    chr20: 64444167,
    chr21: 46709983,
    chr22: 50818468,
    chrX: 156040895,
    chrY: 57227415
};

function CistromeBigWigDataFetcher(HGC, ...args) {
    console.log('here');
    if(!new.target) {
        throw new Error("Uncaught TypeError: Class constructor cannot be invoked without \"new\"");
    }
    const cls = class CistromeBigWigDataFetcherClass {
        constructor(dataConfig) {
            dataConfig.chromSizesUrl = "https://aveit.s3.amazonaws.com/higlass/data/sequence/hg38.chrom.sizes"; //https://s3.amazonaws.com/gosling-lang.org/data/hg38.chrom.sizes';
            this.dataConfig = dataConfig;
            this.bwFileHeader = null;
            this.bwFile = null;
            this.TILE_SIZE = 1024;

            this.errorTxt = "";
            this.dataPromises = [];

            // Prepare chromosome interval information
            const chromosomeSizes = CHROM_SIZE_HG38;
            const chromosomeCumPositions = [];
            const chromosomePositions = {};
            let prevEndPosition = 0;

            Object.keys(chromosomeSizes).forEach((chrStr, i) => {
                const positionInfo = {
                    id: i,
                    chr: chrStr,
                    pos: prevEndPosition
                };

                chromosomeCumPositions.push(positionInfo);
                chromosomePositions[chrStr] = positionInfo;

                prevEndPosition += chromosomeSizes[chrStr];
            });
            this.chromSizes = {
                chrToAbs: (chrom, chromPos) => this.chromSizes.chrPositions[chrom].pos + chromPos,
                cumPositions: chromosomeCumPositions,
                chrPositions: chromosomePositions,
                totalLength: prevEndPosition,
                chromLengths: chromosomeSizes
            };

            this.dataPromises.push(this.loadBBI(dataConfig));
        }

        loadBBI(dataConfig) {
            if (dataConfig.url) {
                this.bwFile = new BigWig({
                    filehandle: new RemoteFile(dataConfig.url)
                });
                return this.bwFile.getHeader().then(h => {
                    this.bwFileHeader = h;
                });
            } else {
                console.error("Please enter a \"url\" field to the data config");
                return null;
            }
        }

        tilesetInfo(callback) {
            this.tilesetInfoLoading = true;

            return Promise.all(this.dataPromises)
                .then(() => {
                    this.tilesetInfoLoading = false;

                    let retVal = {};

                    const totalLength = this.chromSizes.totalLength;

                    retVal = {
                        tile_size: this.TILE_SIZE,
                        max_zoom: Math.ceil(Math.log(totalLength / this.TILE_SIZE) / Math.log(2)),
                        max_width: 2 ** Math.ceil(Math.log(totalLength) / Math.log(2)),
                        min_pos: [0],
                        max_pos: [totalLength]
                    };

                    if (callback) {
                        callback(retVal);
                    }

                    return retVal;
                })
                .catch(err => {
                    this.tilesetInfoLoading = false;

                    console.error(err);

                    if (callback) {
                        callback({
                            error: `Error parsing bigwig: ${err}`
                        });
                    }
                });
        }

        fetchTilesDebounced(receivedTiles, tileIds) {
            const tiles = {};

            const validTileIds = [];
            const tilePromises = [];

            for (const tileId of tileIds) {
                const parts = tileId.split(".");
                const z = parseInt(parts[0], 10);
                const x = parseInt(parts[1], 10);

                if (Number.isNaN(x) || Number.isNaN(z)) {
                    console.warn("Invalid tile zoom or position:", z, x);
                    continue;
                }

                validTileIds.push(tileId);
                tilePromises.push(this.tile(z, x));
            }

            Promise.all(tilePromises).then(values => {
                for (let i = 0; i < values.length; i++) {
                    const validTileId = validTileIds[i];
                    tiles[validTileId] = values[i];
                    tiles[validTileId].tilePositionId = validTileId;
                }

                receivedTiles(tiles);
            });
            // tiles = tileResponseToData(tiles, null, tileIds);
            return tiles;
        }

        tile(z, x) {
            return this.tilesetInfo().then(tsInfo => {
                const tileWidth = +tsInfo.max_width / 2 ** +z;

                const recordPromises = [];

                const tile = {
                    tilePos: [x],
                    tileId: `bigwig.${z}.${x}`,
                    zoomLevel: z
                };

                // get the bounds of the tile
                const minXOriginal = tsInfo.min_pos[0] + x * tileWidth;
                let minX = minXOriginal;
                const maxX = tsInfo.min_pos[0] + (x + 1) * tileWidth;

                const basesPerPixel = this.determineScale(minX, maxX);
                const basesPerBin = (maxX - minX) / this.TILE_SIZE;

                const binStarts = [];
                for (let i = 0; i < this.TILE_SIZE; i++) {
                    binStarts.push(minX + i * basesPerBin);
                }

                const { chromLengths, cumPositions } = this.chromSizes;

                cumPositions.forEach(cumPos => {
                    const chromName = cumPos.chr;
                    const chromStart = cumPos.pos;
                    const chromEnd = cumPos.pos + chromLengths[chromName];

                    let startPos, endPos;

                    if (chromStart <= minX && minX < chromEnd) {
                        // start of the visible region is within this chromosome

                        if (maxX > chromEnd) {
                            // the visible region extends beyond the end of this chromosome
                            // fetch from the start until the end of the chromosome
                            startPos = minX - chromStart;
                            endPos = chromEnd - chromStart;
                            recordPromises.push(
                                this.bwFile
                                    .getFeatures(chromName, startPos, endPos, {
                                        scale: 1 / basesPerPixel
                                    })
                                    .then(values => {
                                        values.forEach(v => {
                                            v["startAbs"] = HGC.utils.chrToAbs(chromName, v.start, this.chromSizes);
                                            v["endAbs"] = HGC.utils.chrToAbs(chromName, v.end, this.chromSizes);
                                        });
                                        return values;
                                    })
                            );

                            minX = chromEnd;
                        } else {
                            startPos = Math.floor(minX - chromStart);
                            endPos = Math.ceil(maxX - chromStart);
                            recordPromises.push(
                                this.bwFile
                                    .getFeatures(chromName, startPos, endPos, {
                                        scale: 1 / basesPerPixel
                                    })
                                    .then(values => {
                                        values.forEach(v => {
                                            v["startAbs"] = HGC.utils.chrToAbs(chromName, v.start, this.chromSizes);
                                            v["endAbs"] = HGC.utils.chrToAbs(chromName, v.end, this.chromSizes);
                                        });
                                        return values;
                                    })
                            );
                            return;
                        }
                    }
                });

                return Promise.all(recordPromises).then(v => {
                    const values = v.flat();

                    var dense = [];
                    for (var i = 0; i < this.TILE_SIZE; i++) {
                        dense.push(null);
                    }

                    // Currently we use the same binning strategy in all cases (basesPerBin =>< basesPerBinInFile)
                    binStarts.forEach((curStart, index) => {
                        if (curStart < minXOriginal || curStart > maxX) {
                            return;
                        }
                        const filtered = values
                            .filter(v => {
                                return curStart >= v.startAbs && curStart < v.endAbs;
                            })
                            .map(v => v.score);
                        dense[index] = filtered.length > 0 ? filtered[0] : null;
                    });

                    tile.min_value = Math.min(...dense);
                    tile.max_value = Math.max(...dense);

                    const dde = new HGC.utils.DenseDataExtrema1D(dense);
                    tile.dense = dense;
                    tile.denseDataExtrema = dde;
                    tile.minNonZero = dde.minNonZeroInTile;
                    tile.maxNonZero = dde.maxNonZeroInTile;
                    return tile;
                });
            });
        }

        // We never want to request more than 1024 * 20 elements from the file.
        determineScale(minX, maxX) {
            const reductionLevels = [1];
            const numRequestedElements = maxX - minX;

            this.bwFileHeader.zoomLevels.forEach(z => {
                reductionLevels.push(z.reductionLevel);
            });

            let level;
            reductionLevels.forEach(rl => {
                if (level) return; // we found one

                const numElementsFromFile = numRequestedElements / rl;
                if (numElementsFromFile <= this.TILE_SIZE * 20) {
                    level = rl;
                }
            });

            // return the highest reductionLevel, if we could not find anything better
            return level || reductionLevels.slice(-1)[0];
        }
    };

    return new cls(...args);
}

CistromeBigWigDataFetcher.config = {
    type: "cistrome-bigwig"
};

export default CistromeBigWigDataFetcher;
