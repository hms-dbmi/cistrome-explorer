/**
 * TODO: 
 * @param {object|array} rowInfo 
 * @param {string} field 
 * @param {string} type 
 */
export function getAggregatedValue(rowInfo, field, type) {
    if(type === "nominal") {
        return Array.isArray(rowInfo) ? rowInfo[0][field] : rowInfo[field];
    }
}