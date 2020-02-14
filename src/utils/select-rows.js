
/**
 * Generate an array of selected row indices based on sort and filter options.
 */
export function selectRows(rowInfoOrig, sortOptions, filterOptions) {
    if(!sortOptions && !filterOptions) {
        // Null means select all rows.
        return null;
    } else {
        return [1, 3, 5];
    }
    // Filter
    // ...
    
    // Sort
    /*let transformedRowInfo = rowInfo.slice();
    if(options.rowSort && options.rowSort.length > 0) {
        let sortOptions = options.rowSort.slice().reverse();
        sortOptions.forEach((d, i) => {
            const { field, type, order } = d;
            if(type === "tree") {
                // Do nothing for the "tree" type.
            } else if(type === "quantitative") {
                transformedRowInfo.sort((a, b) => (a[field] - b[field]) * (order === "ascending" ? 1 : -1));
            } else {
                transformedRowInfo.sort(function(a, b) {
                    let compared = 0, categoryA = a[field].toUpperCase(), categoryB = b[field].toUpperCase();
                    if(categoryA > categoryB) {
                        compared = 1;
                    } else {
                        compared = -1;
                    }
                    return compared * (order === "ascending" ? 1 : -1);
                });
            }
        });
    }*/
}