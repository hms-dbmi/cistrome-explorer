/**
 * Get tree structured data from matrix data.
 * @param {array} inputMatrix Matrix data.
 * @returns {object} Tree structure.
 */
export function matrixToTree(inputMatrix) {
    const group = (matrix, level) => {
        return matrix.reduce((a, h) => {
                if(a.hasOwnProperty(h.levels[level])) {
                    a[h.levels[level]].push(h);
                } else {
                    a[h.levels[level]] = [ h ];
                }
            return a;
        }, {});
    };
    const convert = (matrix, level) => {
        const curr = group(matrix, level);
        return Object.entries(curr).map(([k, v]) => {
            if(v.length === 1) {
                return { name: k, i: v[0].i }
            } else {
                return { name: k, children: convert(v, level+1) }
            }
        });
    };

    const inputMatrixWithI = inputMatrix.map((d, i) => ({ i, levels: d }));
    return { name: "root", children: convert(inputMatrixWithI, 0) };
}