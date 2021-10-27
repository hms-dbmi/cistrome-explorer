/**
 * Get tree structured data from matrix data.
 * @param {array} inputMatrix Matrix data.
 * @returns {object} Tree structure.
 */
export function matrixToTree(inputMatrix) {
    if(inputMatrix.length !== 0 && inputMatrix[0].length !== 0 && typeof inputMatrix[0][0] === "object") {
        // TODO: Ultimately, use the following function when we always encode similarity in dendrogram.
        return matrixToTreeWithDistance(inputMatrix);
    }
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
                return { name: k, i: v[0].i };
            } else {
                return { name: k, children: convert(v, level+1) };
            }
        });
    };
    const inputMatrixWithI = inputMatrix.map((d, i) => ({ i, levels: d }));
    return { name: "root", children: convert(inputMatrixWithI, 0) };
}

/**
 * Get tree structured data from matrix data.
 * @param {array} inputMatrix Matrix data.
 * @returns {object} Tree structure.
 */
export function matrixToTreeWithDistance(inputMatrix) {
    // Reference: https://stackoverflow.com/a/20060322
    const group = (matrix, level) => {
        return matrix.reduce((a, h) => {
            const branch = h.levels[level].name;
            if(a.hasOwnProperty(branch)) {
                a[branch].push(h);
            } else {
                a[branch] = [ h ];
            }
            return a;
        }, {});
    };
    const convert = (matrix, level) => {
        const curr = group(matrix, level);
        return Object.entries(curr).map(([k, v]) => {
            const dist = v[0].levels.find(d => d.name == k).dist;
            if(v.length === 1) {
                return { name: k, dist, i: v[0].i };
            } else {
                return { name: k, dist, children: convert(v, level+1) };
            }
        });
    };
    const inputMatrixWithI = inputMatrix.map((d, i) => ({ i, levels: d }));
    return { name: "root", children: convert(inputMatrixWithI, 0)};
}