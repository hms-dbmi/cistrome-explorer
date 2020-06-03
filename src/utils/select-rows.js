import { matrixToTree } from './tree.js';
import d3 from './d3.js';
import { getAggregatedValue } from './aggregate.js';

/**
 * Generate an array of selected row indices based on filter and sort options.
 * @param {object[]} rowInfo The original/full default-ordered rowInfo array.
 * @param {object} options The track options object, containing sort/filter options.
 * @returns {(number[]|null)} The array of selected indices.
 */
export function selectRows(rowInfo, options) {
    if(options) {   
        // Aggregate
        let aggregatedRowInfo = getAggregatedRowInfo(rowInfo, options.rowAggregate);
        // Hereafter, each datum of `aggregatedRowInfo` can be either [number, Object] or [number[], Object[]].
        
        // Filter
        let filteredRowInfo = aggregatedRowInfo;
        if(options.rowFilter && options.rowFilter.length > 0) {
            let filterInfos = options.rowFilter;

            // Should apply `minSimilarity` filter lastly.
            filterInfos.sort((a, b) => a.minSimilarity ? 1 : b.minSimilarity ? -1 : 0);

            filterInfos.forEach(info => {
                const { field, type, notOneOf, range, subtree, minSimilarity } = info;
                const aggFunction = options.rowInfoAttributes?.find(d => d.field === field)?.aggFunction;
                const isMultipleFields = Array.isArray(field);
                if(type === "nominal") {
                    notOneOf.forEach(one => {
                        filteredRowInfo = filteredRowInfo.filter(
                            d => getAggregatedValue(d[1], field, "nominal", aggFunction).toString().toUpperCase() !== one.toUpperCase()
                        );
                    })
                } else if(type === "quantitative") {
                    const [minCutoff, maxCutoff] = range;
                    if(isMultipleFields) {
                        filteredRowInfo = filteredRowInfo.filter(d => {
                            let sum = 0;
                            field.forEach(f => sum += getAggregatedValue(d[1], f, "quantitative", aggFunction));
                            return sum > minCutoff && sum < maxCutoff;
                        });
                    } else {
                        filteredRowInfo = filteredRowInfo.filter(
                            d => getAggregatedValue(d[1], field, "quantitative", aggFunction) > minCutoff 
                                && getAggregatedValue(d[1], field, "quantitative", aggFunction) < maxCutoff
                        );
                    }
                } else if(type === "tree") {
                    if(options.rowAggregate && options.rowAggregate.length >= 1) {
                        // TODO: How to best support aggregation for `tree`?
                        return;
                    }
                    // `tree` type filter can have both the `subtree` and `minSimilarity` filters in a single `filterInfo`.
                    if(subtree) {
                        filteredRowInfo = filteredRowInfo.filter(d => d[1][field].reduce(
                            // TODO: Remove `h === subtree[i]` when we always encode similarity distance in dendrogram.
                            (a, h, i) => a && (i >= subtree.length || h === subtree[i] || h.name === subtree[i]), true)
                        );
                    } 
                    if (minSimilarity) {
                        filteredRowInfo = filteredRowInfo.filter(
                            // Note that leafs' `dist` values are zero.
                            d => d[1][field].map(d => d.dist).filter(d => d <= minSimilarity).length > 1
                        );
                        filteredRowInfo = filteredRowInfo.filter(
                            // When a row has a branch that is closer than `minSimilarity`, but there is
                            // no other rows to be connected (due to a certain filtering status), the row has no connections.
                            // (Refer to https://github.com/hms-dbmi/cistrome-explorer/pull/241)
                            d => {
                                const treeData = d[1][field];
                                for(let i = 0; i < treeData.length; i++) {
                                    const branch = treeData[i];
                                    if(branch.dist > minSimilarity) continue;
                                    const numOfFound = filteredRowInfo.filter(
                                        info => info[1][field].find(b => b.dist === branch.dist && b.name === branch.name) !== undefined
                                    ).length;
                                    if(numOfFound >= 2) {   // Using `2` since filteredRowInfo contains `treeData` itself.
                                        return true;
                                    }
                                }
                                return false;
                            }
                        );
                    }
                }
            });
        }
        let transformedRowInfo = filteredRowInfo;
        if(!transformedRowInfo || transformedRowInfo.length == 0) {
            return [];
        }

        // Sort
        if(options.rowSort && options.rowSort.length > 0) {
            let sortOptions = options.rowSort.slice().reverse();
            sortOptions.forEach((d) => {
                const { field, type, order } = d;
                const isMultipleFields = Array.isArray(field);
                const aggFunction = options.rowInfoAttributes?.find(d => d.field === field)?.aggFunction;
                if(type === "tree") {
                    if(options.rowAggregate && options.rowAggregate.length >= 1) {
                        // TODO: How to best support aggregation for `tree`?
                        return;
                    }
                    const hierarchyData = matrixToTree(filteredRowInfo.map(d => d[1][field]));
                    const root = d3.hierarchy(hierarchyData);
                    const leaves = root.leaves().map(l => l.data.i);
                    transformedRowInfo = leaves.map((i) => transformedRowInfo[i]);
                } else if(type === "quantitative") {
                    transformedRowInfo.sort((a, b) => {
                        if(isMultipleFields) {
                            let sumA = 0, sumB = 0;
                            field.forEach(f => sumA += getAggregatedValue(a[1], f, "quantitative", aggFunction));
                            field.forEach(f => sumB += getAggregatedValue(b[1], f, "quantitative", aggFunction));
                            return (sumA - sumB) * (order === "ascending" ? 1 : -1);
                        } else {
                            const valueA = getAggregatedValue(a[1], field, "quantitative", aggFunction);
                            const valueB = getAggregatedValue(b[1], field, "quantitative", aggFunction);
                            return (valueA - valueB) * (order === "ascending" ? 1 : -1);
                        }
                    });
                } else if(type === "nominal") {
                    transformedRowInfo.sort(function(a, b) {
                        let compared = 0;
                        const categoryA = getAggregatedValue(a[1], field, "nominal", aggFunction).toString().toUpperCase();
                        const categoryB = getAggregatedValue(b[1], field, "nominal", aggFunction).toString().toUpperCase();
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

        // Zoom
        if(options.rowZoom) {
            const { level: zoomLevel, top: zoomTop, numRows } = options.rowZoom;
            const rowUnit = 1.0 / numRows;
            if(0 < zoomLevel && zoomLevel <= 1.0) {
                const topIndex = Math.floor(zoomTop * numRows);
                const numRowsZoomed = Math.ceil(zoomLevel / rowUnit);
                transformedRowInfo = transformedRowInfo.slice(topIndex, topIndex + numRowsZoomed);
            }
        }

        return transformedRowInfo.map(d => d[0]);
    }

    // Null means select all rows.
    return null;
}

/**
 * Aggregate row information using `rowAggregate` options.
 * Data formats are changed from [number, Object] to [number[], Object[]] for aggregated rows.
 * @param {object[]} rowInfo The original/full default-ordered rowInfo array.
 * @param {object} rowAggregate The track option object for aggregating rows.
 * @returns {object[]} The aggregated `rowInfo`.
 */
export function getAggregatedRowInfo(rowInfo, rowAggregate) {
    let aggregatedRowInfo = Array.from(rowInfo.entries());
    if(rowAggregate && rowAggregate.length > 0) {
        const aggregateOptions = rowAggregate.slice().reverse();
        aggregateOptions.forEach(d => {
            const { field, type, oneOf } = d;
            if(type === "nominal") {
                oneOf.forEach(o => {
                    const matchings = aggregatedRowInfo.filter(t => t[1][field] === o);
                    const notMatchings = aggregatedRowInfo.filter(t => t[1][field] !== o);
                    const matchingIndices = matchings.map(t => t[0]);
                    const matchingRowInfos = matchings.map(t => t[1]);

                    if(matchings?.length == 0) return;

                    aggregatedRowInfo = notMatchings;
                    aggregatedRowInfo.push([matchingIndices, matchingRowInfos]);
                });
            }
        });
    }
    return aggregatedRowInfo;
}

/**
 * Generate an array of highlighted row indices based on conditions.
 * @param {object[]} rowInfo The original/full default-ordered rowInfo array.
 * @param {string} field The attribute on which to search.
 * @param {string} type The data type contained in the field value.
 * @param {string|array} conditions The conditions to check for highlighting.
 * @param {object} options The track options object, containing aggregation options.
 * @returns {number[]} The array of highlit indices.
 */
export function highlightRowsFromSearch(rowInfo, field, type, conditions, options) {
    // Aggregation rows in advance to apply highlighting to the aggregated rows.
    const aggregatedRowInfo = getAggregatedRowInfo(rowInfo, options.rowAggregate);
    const aggFunction = options.rowInfoAttributes?.find(d => d.field === field)?.aggFunction;
    let newHighlitRows = [];
    if(conditions === "") {
        newHighlitRows = [];
    } else if(type === "nominal") {
        const filteredRows = aggregatedRowInfo.filter(
            d => !getAggregatedValue(d[1], field, 'nominal', aggFunction).toString().toUpperCase().includes(conditions.toUpperCase())
        );
        newHighlitRows = filteredRows.map(d => d[0]);
    } else if(type === "quantitative") {
        const [minCutoff, maxCutoff] = conditions;
        let filteredRows;
        if(Array.isArray(field)) {
            filteredRows = aggregatedRowInfo.filter(d => {
                let sum = 0;
                field.forEach(f => sum += getAggregatedValue(d[1], f, 'quantitative', aggFunction));
                return sum < minCutoff || sum > maxCutoff;
            });
        } else {
            filteredRows = aggregatedRowInfo.filter(
                d => getAggregatedValue(d[1], field, 'quantitative', aggFunction) < minCutoff 
                    || getAggregatedValue(d[1], field, 'quantitative', aggFunction) > maxCutoff
            );
        }
        newHighlitRows = filteredRows.map(d => d[0]);
    } else if(type === "tree") {
        if(Array.isArray(conditions)) {
            const subtree = conditions;
            const filteredRows = aggregatedRowInfo.filter(
                d => getAggregatedValue(d[1], field, 'tree', aggFunction).reduce(
                // TODO: Remove `h === subtree[i]` when we always encode similarity distance in dendrogram.
                (a, h, i) => a || (i < subtree.length && h !== subtree[i] && h.name !== subtree[i]), false)
            );
            newHighlitRows = filteredRows.map(d => d[0]);
        } else {
            const minSimilarity = conditions;
            const filteredRows = aggregatedRowInfo.filter(
                // Note that leafs' `dist` values are zero.
                d => getAggregatedValue(d[1], field, 'tree', aggFunction).map(d => d.dist).filter(d => d <= minSimilarity).length <= 1
            );
            newHighlitRows = filteredRows.map(d => d[0]);
        }
    }
    return newHighlitRows;
}