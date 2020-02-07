/**
 * Get tree structured data from matrix data.
 * @param {array} inputMatrix Matrix data.
 * @returns {object} Tree structure.
 */
export function matrixToTree(inputMatrix) {
  
    const group = (matrix, level) => {
        return matrix.reduce((a, h) => {
            if(a.hasOwnProperty(h[level])) {
               a[h[level]].push(h);
            } else {
                a[h[level]] = [ h ];
            }
        return a;
      }, {});
    };
    const convert = (matrix, level) => {
        const curr = group(matrix, level);
        return Object.entries(curr).map(([k, v]) => {
            if(v.length === 1) {
                return { name: k }
            } else {
                return { name: k, children: convert(v, level+1) }
            }
        });
    };
    
    return { name: "root", children: convert(inputMatrix, 0) };
}