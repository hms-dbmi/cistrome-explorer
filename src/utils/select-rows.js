import { matrixToTree } from './tree.js';
import d3 from './d3.js';

/**
 * Generate an array of selected row indices based on sort and filter options.
 * @param {object[]} rowInfo The original/full default-ordered rowInfo array.
 * @param {object} options The track options object, containing sort/filter options.
 * @returns {(number[]|null)} The array of selected indices.
 */
export function selectRows(rowInfo, options) {
    if(options) {
        // Filter
        // ...

        // Sort
        let transformedRowInfo = Array.from(rowInfo.entries());
        if(options.rowSort && options.rowSort.length > 0) {
            let sortOptions = options.rowSort.slice().reverse();
            sortOptions.forEach((d) => {
                const { field, type, order } = d;
                if(type === "tree") {
                    const hierarchyData = matrixToTree(rowInfo.map(d => d[field]));
                    const root = d3.hierarchy(hierarchyData);
                    const leaves = root.leaves().map(l => l.data.i);
                    transformedRowInfo = leaves.map((i) => transformedRowInfo[i]);
                } else if(type === "quantitative") {
                    transformedRowInfo.sort((a, b) => (a[1][field] - b[1][field]) * (order === "ascending" ? 1 : -1));
                } else if(type === "nominal") {
                    transformedRowInfo.sort(function(a, b) {
                        let compared = 0, categoryA = a[1][field].toUpperCase(), categoryB = b[1][field].toUpperCase();
                        if(categoryA > categoryB) {
                            compared = 1;
                        } else {
                            compared = -1;
                        }
                        return compared * (order === "ascending" ? 1 : -1);
                    });
                }
            });
        }
        return transformedRowInfo.map(d => d[0]);
    }

    // Null means select all rows.
    return null;
}

/**
 * Generate an array of highlighted row indices based on search keyword information.
 * @param {object[]} rowInfo The original/full default-ordered rowInfo array.
 * @param {string} field The attribute on which to search.
 * @param {string} type The data type contained in the field value.
 * @param {string} contains The search term.
 * @returns {number[]} The array of highlit indices.
 */
export function highlightRowsFromSearch(rowInfo, field, type, contains) {
    let newHighlitRows = [];
    if(contains === "") {
        newHighlitRows = [];
    } else if(type === "nominal") {
        const rowsWithIndex = Array.from(rowInfo.entries());
        const filteredRows = rowsWithIndex.filter(d => d[1][field].toUpperCase().includes(contains.toUpperCase()));
        newHighlitRows = filteredRows.map(d => d[0]);
    } else if(type === "quantitative") {
        // TODO: Better deal with quantitative data. Need to update Wrapper options for this.
        // refer vega filter, such as lt: https://vega.github.io/vega-lite/docs/filter.html
        const rowsWithIndex = Array.from(rowInfo.entries());
        const filteredRows = rowsWithIndex.filter(d => d[1][field].toString().includes(contains));
        newHighlitRows = filteredRows.map(d => d[0]);
    }
    return newHighlitRows;
}