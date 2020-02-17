/* eslint-env node */

import { 
    matrixToTree
} from './tree.js';

describe('Functions for manipulating hierarchy data structures', () => {
    it('Should convert a matrix-like structure to a tree-like structure', () => {
        const matrix = [
            ["i-0", "i-0", "i-0", "i-0", "4-1"],
            ["i-0", "i-0", "i-0", "i-0", "4-2"],
            ["i-0", "i-0", "i-0", "3-4"],
            ["i-0", "i-0", "2-4"],
            ["i-0", "i-4", "i-4", "i-4", "4-3"],
            ["i-0", "i-4", "i-4", "i-4", "4-4"],
            ["i-0", "i-4", "i-4", "i-4", "4-5"],
            ["i-0", "i-4", "i-4", "i-11", "4-6"],
            ["i-0", "i-4", "i-4", "i-11", "4-7"],
            ["i-0", "i-4", "i-4", "i-11", "4-8"],
            ["i-0", "i-4", "i-29", "3-5"],
            ["i-0", "i-4", "i-29", "3-6"]
        ];
        const tree = matrixToTree(matrix);
        const expectedTree = {"name":"root","children":[{"name":"i-0","children":[{"name":"i-0","children":[{"name":"i-0","children":[{"name":"i-0","children":[{"name":"4-1"},{"name":"4-2"}]},{"name":"3-4"}]},{"name":"2-4"}]},{"name":"i-4","children":[{"name":"i-4","children":[{"name":"i-4","children":[{"name":"4-3"},{"name":"4-4"},{"name":"4-5"}]},{"name":"i-11","children":[{"name":"4-6"},{"name":"4-7"},{"name":"4-8"}]}]},{"name":"i-29","children":[{"name":"3-5"},{"name":"3-6"}]}]}]}]};
        expect(tree).toEqual(expectedTree);
    });

    it('Should convert a matrix-like structure to a tree-like structure, after filtering', () => {
        const matrix = [
            ["i-0", "i-0", "i-0", "i-0", "4-1"],
            ["i-0", "i-0", "i-0", "i-0", "4-2"],
            ["i-0", "i-4", "i-4", "i-11", "4-8"],
            ["i-0", "i-4", "i-29", "3-5"]
        ];
        const tree = matrixToTree(matrix);
        const expectedTree = {"name":"root","children":[{"name":"i-0","children":[{"name":"i-0","children":[{"name":"i-0","children":[{"name":"i-0","children":[{"name":"4-1"},{"name":"4-2"}]}]}]},{"name":"i-4","children":[{"name":"i-4"},{"name":"i-29"}]}]}]};
        expect(tree).toEqual(expectedTree);
    });
});
