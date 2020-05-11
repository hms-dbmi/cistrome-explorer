/* eslint-env node */

import { 
    matrixToTree,
    matrixToTreeWithDistance
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
        const expectedTree = {"name":"root","children":[{"name":"i-0","children":[{"name":"i-0","children":[{"name":"i-0","children":[{"name":"i-0","children":[{i: 0, "name":"4-1"},{i: 1, "name":"4-2"}]},{i: 2, "name":"3-4"}]},{i: 3, "name":"2-4"}]},{"name":"i-4","children":[{"name":"i-4","children":[{"name":"i-4","children":[{i: 4, "name":"4-3"},{i: 5, "name":"4-4"},{i: 6, "name":"4-5"}]},{"name":"i-11","children":[{i: 7, "name":"4-6"},{i: 8, "name":"4-7"},{i: 9, "name":"4-8"}]}]},{"name":"i-29","children":[{i: 10, "name":"3-5"},{i: 11, "name":"3-6"}]}]}]}]};
        expect(tree).toEqual(expectedTree);
    });

    it('Should convert a matrix-like structure with distance values to a tree-like structure', () => {
        const matrix = [
            [{"name": "i-0", "dist" :2}, {"name": "i-10", "dist" :1.632}, {"name": "i-11", "dist" :1.000}, {"name": "3", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-10", "dist" :1.632}, {"name": "i-16", "dist" :0.900}, {"name": "4", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-10", "dist" :1.632}, {"name": "i-11", "dist" :1.000}, {"name": "i-12", "dist" :0.800}, {"name": "i-13", "dist" :0.650}, {"name": "i-14", "dist" :0.350}, {"name": "5", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-10", "dist" :1.632}, {"name": "i-16", "dist" :0.900}, {"name": "6", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-2", "dist" :0.204}, {"name": "i-3", "dist" :0.194}, {"name": "7", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-2", "dist" :0.204}, {"name": "i-3", "dist" :0.194}, {"name": "i-4", "dist" :0.140}, {"name": "8", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-2", "dist" :0.204}, {"name": "i-6", "dist" :0.134}, {"name": "i-7", "dist" :0.100}, {"name": "9", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-10", "dist" :1.632}, {"name": "i-11", "dist" :1.000}, {"name": "i-12", "dist" :0.800}, {"name": "10", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-2", "dist" :0.204}, {"name": "i-6", "dist" :0.134}, {"name": "i-8", "dist" :0.050}, {"name": "i-9", "dist" :0.010}, {"name": "11", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-2", "dist" :0.204}, {"name": "i-6", "dist" :0.134}, {"name": "i-8", "dist" :0.050}, {"name": "i-9", "dist" :0.010}, {"name": "12", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-2", "dist" :0.204}, {"name": "i-6", "dist" :0.134}, {"name": "i-8", "dist" :0.050}, {"name": "13", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-2", "dist" :0.204}, {"name": "i-6", "dist" :0.134}, {"name": "i-7", "dist" :0.100}, {"name": "14", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-2", "dist" :0.204}, {"name": "i-3", "dist" :0.194}, {"name": "i-4", "dist" :0.140}, {"name": "i-5", "dist" :0.110}, {"name": "15", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-10", "dist" :1.632}, {"name": "i-11", "dist" :1.000}, {"name": "i-12", "dist" :0.800}, {"name": "i-13", "dist" :0.650}, {"name": "i-14", "dist" :0.350}, {"name": "i-15", "dist" :0.250}, {"name": "16", "dist": 0}],
            [{"name": "i-0", "dist" :2}, {"name": "i-10", "dist" :1.632}, {"name": "i-11", "dist" :1.000}, {"name": "i-12", "dist" :0.800}, {"name": "i-13", "dist" :0.650}, {"name": "17", "dist": 0}]
        ];
        const tree = matrixToTree(matrix);
        const expectedTree = {"name":"root","children":[{"name":"i-0","dist":2,"children":[{"name":"i-10","dist":1.632,"children":[{"name":"i-11","dist":1,"children":[{"name":"3","dist":0,"i":0},{"name":"i-12","dist":0.8,"children":[{"name":"10","dist":0,"i":7},{"name":"i-13","dist":0.65,"children":[{"name":"17","dist":0,"i":14},{"name":"i-14","dist":0.35,"children":[{"name":"5","dist":0,"i":2},{"name":"i-15","dist":0.25,"i":13}]}]}]}]},{"name":"i-16","dist":0.9,"children":[{"name":"4","dist":0,"i":1},{"name":"6","dist":0,"i":3}]}]},{"name":"i-2","dist":0.204,"children":[{"name":"i-3","dist":0.194,"children":[{"name":"7","dist":0,"i":4},{"name":"i-4","dist":0.14,"children":[{"name":"8","dist":0,"i":5},{"name":"i-5","dist":0.11,"i":12}]}]},{"name":"i-6","dist":0.134,"children":[{"name":"i-7","dist":0.1,"children":[{"name":"9","dist":0,"i":6},{"name":"14","dist":0,"i":11}]},{"name":"i-8","dist":0.05,"children":[{"name":"13","dist":0,"i":10},{"name":"i-9","dist":0.01,"children":[{"name":"11","dist":0,"i":8},{"name":"12","dist":0,"i":9}]}]}]}]}]}]};
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
        const expectedTree = {"name":"root","children":[{"name":"i-0","children":[{"name":"i-0","children":[{"name":"i-0","children":[{"name":"i-0","children":[{i: 0, "name":"4-1"},{i: 1, "name":"4-2"}]}]}]},{"name":"i-4","children":[{i: 2, "name":"i-4"},{i: 3, "name":"i-29"}]}]}]};
        expect(tree).toEqual(expectedTree);
    });
});
