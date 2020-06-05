import d3 from './d3.js';

/**
 * Get aggregated values using specified function, such as a mean of quantitative values.
 * @param {object|array} rowInfo Object or array of objects containing information of rows.
 * @param {string} field The name of data field.
 * @param {string} type The type of data field.
 * @param {string} func A function to aggregate values.
 */
export function getAggregatedValue(rowInfo, field, type, func) {
    const isAggregated = Array.isArray(rowInfo);
    if(!isAggregated) {
        return rowInfo[field];
    }
    const funcName = func ? func : "default";

    // Mapping to the function that accepts an array of values and returns a single representative value.
    const aggFuncMapping = {
        quantitative: {
            max: d3.max,
            min: d3.min,
            mean: d3.mean,
            sum: d3.sum,
            default: d3.mean
        },
        nominal: {
            mostCommon: values => {
                const counts = {};
                values.forEach(d => { counts[d] = 1 + (counts[d] || 0) });
                return Object.keys(counts).reduce((a, b) => counts[a] > counts[b] ? a : b);
            },
            leastCommon: values => {
                const counts = {};
                values.forEach(d => { counts[d] = 1 + (counts[d] || 0) });
                return Object.keys(counts).reduce((a, b) => counts[a] < counts[b] ? a : b);
            },
            concat: values => Array.from(new Set(values)).join(', '),
            default: values => Array.from(new Set(values)).join(', ')
        },
        tree: {
            // `aggFunction` does not mean anything for `tree`.
            // TODO: How to best aggregate?
            default: d => d
        }
    }

    const aggregatedValue = aggFuncMapping?.[type]?.[funcName];
    if(!aggregatedValue) {
        console.warn(`The aggFunction is ill-defined for ${field}.`);
        return rowInfo.map(d => d[field]);
    }
    return aggregatedValue(rowInfo.map(d => d[field]));
}