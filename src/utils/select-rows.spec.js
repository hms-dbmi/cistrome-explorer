/* eslint-env node */

import { 
    selectRows
} from './select-rows.js';

describe('Helper functions for producing arrays of row indices for track.options.selectRows', () => {
    it('Should produce correct selectRows array after sorting on quantitative attribute', () => {
        const rowInfo = [
            {"r1":4},
            {"r1":73},
            {"r1":59},
            {"r1":35},
            {"r1":5},
            {"r1":41},
            {"r1":95},
            {"r1":54},
            {"r1":77},
            {"r1":53},
            {"r1":33},
            {"r1":87},
            {"r1":46},
            {"r1":98},
            {"r1":56}
        ];
        const options = {
            rowSort: [{
                field: "r1",
                order: "descending",
                type: "quantitative"
            }]
        };
        const selectedRows = selectRows(rowInfo, options);
        expect(selectedRows).toEqual([13,6,11,8,1,2,14,7,9,12,5,3,10,4,0]);
    });

    it('Should produce correct selectRows array after sorting on nominal attribute', () => {
        const rowInfo = [
            {"t":"Blood"},
            {"t":"Tonsil"},
            {"t":"Blood"},
            {"t":"Tonsil"},
            {"t":"Liver"},
            {"t":"Fetal Liver"},
            {"t":"Bone Marrow"},
            {"t":"Blood"},
            {"t":"Bone Marrow"},
            {"t":"Bone Marrow"},
            {"t":"Bone Marrow"},
            {"t":"Bone Marrow"},
            {"t":"Liver"},
            {"t":"Blood"},
            {"t":"Blood"}
        ];
        const options = {
            rowSort: [{
                field: "t",
                order: "ascending",
                type: "nominal"
            }]
        };
        const selectedRows = selectRows(rowInfo, options);
        expect(selectedRows.slice(0, 5)).toContain(0);
        expect(selectedRows.slice(0, 5)).toContain(2);
        expect(selectedRows.slice(0, 5)).toContain(7);
        expect(selectedRows.slice(0, 5)).toContain(13);
        expect(selectedRows.slice(0, 5)).toContain(14);
    });

    it('Should produce correct selectRows array after filtering on tree attribute using `substree`', () => {
        const rowInfo = [
            {"id": "a","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"b2","dist":0.1},{"name":"node-a","dist":0}]},
            {"id": "b","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"b2","dist":0.1},{"name":"node-b","dist":0}]},
            {"id": "c","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"node-c","dist":0}]},
            {"id": "d","t":[{"name":"root","dist":1},{"name":"b3","dist":0.3},{"name":"node-d","dist":0}]},
            {"id": "e","t":[{"name":"root","dist":1},{"name":"b3","dist":0.3},{"name":"node-e","dist":0}]},
        ];
        const options = {
            rowFilter: [{
                field: "t",
                type: "tree",
                subtree: ["root", "b1"]
            }]
        };
        const selectedRows = selectRows(rowInfo, options);
        expect(selectedRows).toContain(0);
        expect(selectedRows).toContain(1);
        expect(selectedRows).toContain(2);
        expect(selectedRows).not.toContain(3);
        expect(selectedRows).not.toContain(4);
    });

    it('Should produce correct selectRows array after filtering on tree attribute using `minSimilarity`', () => {
        const rowInfo = [
            {"id": "a","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"b2","dist":0.1},{"name":"node-a","dist":0}]},
            {"id": "b","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"b2","dist":0.1},{"name":"node-b","dist":0}]},
            {"id": "c","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"node-c","dist":0}]},
            {"id": "d","t":[{"name":"root","dist":1},{"name":"b3","dist":0.3},{"name":"node-d","dist":0}]},
            {"id": "e","t":[{"name":"root","dist":1},{"name":"b3","dist":0.3},{"name":"node-e","dist":0}]},
        ];
        const options = {
            rowFilter: [{
                field: "t",
                type: "tree",
                minSimilarity: 0.2
            }]
        };
        const selectedRows = selectRows(rowInfo, options);
        expect(selectedRows).toContain(0);
        expect(selectedRows).toContain(1);
        expect(selectedRows).not.toContain(2);
        expect(selectedRows).not.toContain(3);
        expect(selectedRows).not.toContain(4);
    });

    it('Should produce correct selectRows array after filtering on tree attribute using `subtree` and `minSimilarity`', () => {
        const rowInfo = [
            {"id": "a","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"b2","dist":0.1},{"name":"node-a","dist":0}]},
            {"id": "b","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"b2","dist":0.1},{"name":"node-b","dist":0}]},
            {"id": "c","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"node-c","dist":0}]},
            {"id": "d","t":[{"name":"root","dist":1},{"name":"b3","dist":0.3},{"name":"node-d","dist":0}]},
            {"id": "e","t":[{"name":"root","dist":1},{"name":"b3","dist":0.3},{"name":"node-e","dist":0}]},
        ];
        const options = {
            rowFilter: [{
                field: "t",
                type: "tree",
                minSimilarity: 0.2,
                subtree: ["root", "b1"]
            }]
        };
        const selectedRows = selectRows(rowInfo, options);
        expect(selectedRows).toContain(0);
        expect(selectedRows).toContain(1);
        expect(selectedRows).not.toContain(2);
        expect(selectedRows).not.toContain(3);
        expect(selectedRows).not.toContain(4);
    });

    it('Should produce correct selectRows array after filtering on tree and other attributes together', () => {
        const rowInfo = [
            {"id": "a","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"b2","dist":0.1},{"name":"node-a","dist":0}]},
            {"id": "b","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"b2","dist":0.1},{"name":"node-b","dist":0}]},
            {"id": "c","t":[{"name":"root","dist":1},{"name":"b1","dist":0.5},{"name":"node-c","dist":0}]},
            {"id": "d","t":[{"name":"root","dist":1},{"name":"b3","dist":0.3},{"name":"node-d","dist":0}]},
            {"id": "e","t":[{"name":"root","dist":1},{"name":"b3","dist":0.3},{"name":"node-e","dist":0}]},
        ];
        const options = {
            rowFilter: [
                {
                    field: "t",
                    type: "tree",
                    minSimilarity: 0.2,
                    subtree: ["root", "b1"]
                },
                {
                    field: "id",
                    type: "nominal",
                    notOneOf: ["a"] // This should also remove node "b" since there is no connection.
                }
            ]
        };
        let selectedRows = selectRows(rowInfo, options);
        expect(selectedRows).not.toContain(0);
        expect(selectedRows).not.toContain(1);
        expect(selectedRows).not.toContain(2);
        expect(selectedRows).not.toContain(3);
        expect(selectedRows).not.toContain(4);

        // Should work identically independant to the order of filters suggested.
        const optionsRev = {
            rowFilter: [
                {
                    field: "id",
                    type: "nominal",
                    notOneOf: ["a"] // This should also remove node "b" since there is no connection.
                },
                {
                    field: "t",
                    type: "tree",
                    minSimilarity: 0.2,
                    subtree: ["root", "b1"]
                }
            ]
        };
        selectedRows = selectRows(rowInfo, optionsRev);
        expect(selectedRows).not.toContain(0);
        expect(selectedRows).not.toContain(1);
        expect(selectedRows).not.toContain(2);
        expect(selectedRows).not.toContain(3);
        expect(selectedRows).not.toContain(4);
    });
});