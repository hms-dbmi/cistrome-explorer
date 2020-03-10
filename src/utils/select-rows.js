import { matrixToTree } from './tree.js';
import d3 from './d3.js';

/**
 * Generate an array of selected row indices based on filter and sort options.
 * @param {object[]} rowInfo The original/full default-ordered rowInfo array.
 * @param {object} options The track options object, containing sort/filter options.
 * @returns {(number[]|null)} The array of selected indices.
 */
export function selectRows(rowInfo, options) {
    if(options) {
        // Filter
        let filteredRowInfo = Array.from(rowInfo.entries());
        if(options.rowFilter && options.rowFilter.length > 0) {
            const filterInfos = options.rowFilter;
            filterInfos.forEach(info => {
                const { field, type, contains } = info;
                const isMultipleFields = Array.isArray(field);
                if(type === "nominal") {
                    filteredRowInfo = filteredRowInfo.filter(d => d[1][field].toString().toUpperCase().includes(contains.toUpperCase()));
                } else if(type === "quantitative") {
                    // TODO: Better deal with quantitative data. Need to update Wrapper options for this.
                    // refer vega filter, such as lt: https://vega.github.io/vega-lite/docs/filter.html
                    if(isMultipleFields) {
                        //...
                    } else {
                        filteredRowInfo = filteredRowInfo.filter(d => d[1][field].toString().includes(contains));
                    }
                } else if(type === "tree") {
                    filteredRowInfo = filteredRowInfo.filter(d => d[1][field].reduce((a, h, i) => a && (i >= contains.length || h === contains[i]), true));
                }
            });
        }
        // Sort
        let transformedRowInfo = filteredRowInfo;
        if(options.rowSort && options.rowSort.length > 0) {
            let sortOptions = options.rowSort.slice().reverse();
            sortOptions.forEach((d) => {
                const { field, type, order } = d;
                const isMultipleFields = Array.isArray(field);
                if(type === "tree") {
                    const hierarchyData = matrixToTree(filteredRowInfo.map(d => d[1][field]));
                    const root = d3.hierarchy(hierarchyData);
                    const leaves = root.leaves().map(l => l.data.i);
                    transformedRowInfo = leaves.map((i) => transformedRowInfo[i]);
                } else if(type === "quantitative") {
                    transformedRowInfo.sort((a, b) => {
                        if(isMultipleFields) {
                            let sumA = 0, sumB = 0;
                            field.forEach(f => sumA += a[1][f]);
                            field.forEach(f => sumB += b[1][f]);
                            return (sumA - sumB) * (order === "ascending" ? 1 : -1);
                        } else {
                            return (a[1][field] - b[1][field]) * (order === "ascending" ? 1 : -1);
                        }
                    });
                } else if(type === "nominal") {
                    transformedRowInfo.sort(function(a, b) {
                        let compared = 0, categoryA = a[1][field].toString().toUpperCase(), categoryB = b[1][field].toString().toUpperCase();
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
    if(Array.isArray(field)) return newHighlitRows;
    if(contains === "") {
        newHighlitRows = [];
    } else if(type === "nominal") {
        const rowsWithIndex = Array.from(rowInfo.entries());
        const filteredRows = rowsWithIndex.filter(d => !d[1][field].toString().toUpperCase().includes(contains.toUpperCase()));
        newHighlitRows = filteredRows.map(d => d[0]);
    } else if(type === "quantitative") {
        // TODO: Better deal with quantitative data. Need to update Wrapper options for this.
        // refer vega filter, such as lt: https://vega.github.io/vega-lite/docs/filter.html
        const rowsWithIndex = Array.from(rowInfo.entries());
        const filteredRows = rowsWithIndex.filter(d => !d[1][field].toString().includes(contains));
        newHighlitRows = filteredRows.map(d => d[0]);
    }
    return newHighlitRows;
}