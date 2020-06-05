import d3 from './d3.js';

/**
 * Get aggregated values using specified function, such as a mean of quantitative values.
 * @param {object|array} rowInfo Object or array of objects containing information of rows.
 * @param {string} field The name of data field.
 * @param {string} type The type of data field.
 * @param {string} fn A function to aggregate values.
 */
export function getAggregatedValue(rowInfo, field, type, fn = "sum") {
    const isAggregated = Array.isArray(rowInfo);
    if(!isAggregated) {
        return rowInfo[field];
    }

    // Aggregated values
    if(type === "nominal" || type === "link") {
        if(fn === "mostCommon" || fn === "leastCommon") {
            const counts = {};
            rowInfo.forEach(d => {
                counts[d[field]] = 1 + (counts[d[field]] || 0);
            });
            const oneCategory = Object.keys(counts).reduce((a, b) => {
                if(fn === "mostCommon") {
                    return counts[a] > counts[b] ? a : b;
                } else {
                    return counts[a] > counts[b] ? b : a;
                }
            });
            return oneCategory;
        } else {
            // "concat"
            const categories = Array.from(new Set(rowInfo.map(d => d[field])));
            return categories.join(", ");
        } 
    } else if(type === "quantitative") {
        if(fn === "max") {
            return d3.max(rowInfo.map(d => d[field]));
        } else if(fn === "mean") {
            return d3.mean(rowInfo.map(d => d[field]));
        } else {
            return d3.sum(rowInfo.map(d => d[field]));
        } 
    } else if(type === "tree") {
        // `aggFunction` does not mean anything for `tree`.
        // TODO: How to best aggregate?
        return rowInfo.map(d => d[field]);
    }
    return rowInfo.map(d => d[field]);
}