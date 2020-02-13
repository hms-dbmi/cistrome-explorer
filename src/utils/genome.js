import { ChromosomeInfo } from 'higlass';

const CHROM_INFO_SERVER = "https://higlass.io/api/v1";
const CHROM_SIZES_URL = `${CHROM_INFO_SERVER}/chrom-sizes/`;
const AVAILABLE_CHROM_SIZES_URL = `${CHROM_INFO_SERVER}/available-chrom-sizes/`;

/**
 * Convert absolute genomic positions to relative genomic positions.
 * See https://docs.higlass.io/javascript_api.html#obtaining-ordered-chromosome-info.
 * @param {string} coordSystem The genome assembly, corresponds to the value of `coordSystem`.
 * @param {number} absStart The genome absolute position of the start of the interval.
 * @param {number} absEnd The genome absolute position of the end of the interval.
 * @returns {promise<array>} A promise that resolves with `[[chrStartName, chrStartPos], [chrEndName, chrEndPos]]`.
 */
export function resolveIntervalCoordinates(coordSystem, absStart, absEnd) {
    return new Promise((resolve, reject) => {
        fetch(AVAILABLE_CHROM_SIZES_URL, {
            credentials: 'same-origin',
            cache: 'force-cache'
        })
        .then(response => {
            return response.json();
        })
        .then(data => {
            const coordSystemInfo = data.results.find(d => d.coordSystem === coordSystem);
            ChromosomeInfo(
                `${CHROM_SIZES_URL}?id=${coordSystemInfo.uuid}`,
                (chromInfo) => {
                    const chrStart = chromInfo.absToChr(absStart);
                    const chrEnd = chromInfo.absToChr(absEnd);
                    resolve([chrStart, chrEnd]);
                }
            );
        })
        .catch(() => {
            reject();
        });
    });
}