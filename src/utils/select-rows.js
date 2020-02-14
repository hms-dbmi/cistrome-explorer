
/**
 * Generate an array of selected row indices based on sort and filter options.
 * @param {object[]} rowInfo The original/full default-ordered rowInfo array.
 * @param {object} options The track options object, containing sort/filter options.
 */
export function selectRows(rowInfo, options) {
    if(!options) {
        // Null means select all rows.
        return null;
    } else {
        // Filter
        // ...

        console.log(options);

        // Sort
        let transformedRowInfo = Array.from(rowInfo.entries());
        if(options.rowSort && options.rowSort.length > 0) {
            let sortOptions = options.rowSort.slice().reverse();
            sortOptions.forEach((d, i) => {
                const { field, type, order } = d;
                if(type === "tree") {
                    // Do nothing for the "tree" type.
                } else if(type === "quantitative") {
                    transformedRowInfo.sort((a, b) => (a[1][field] - b[1][field]) * (order === "ascending" ? 1 : -1));
                } else {
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
}