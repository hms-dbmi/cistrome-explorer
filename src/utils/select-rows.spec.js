/* eslint-env node */

import { 
    selectRows, getAggregatedRowInfo
} from './select-rows.js';
import { getAggregatedValue } from './aggregate.js';

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

    it('Should aggregate rows properly after using nominal attribute', () => {
        const rowInfo = [
            {"c1":"Blood", "r1":1, "c2":"A"},
            {"c1":"Blood", "r1":2, "c2":"B"},
            {"c1":"Blood", "r1":3, "c2":"A"},
            {"c1":"Blood", "r1":4, "c2":"B"},
            {"c1":"Blood", "r1":5, "c2":"A"},
            {"c1":"Bone Marrow", "r1":10, "c2":"A"},
            {"c1":"Bone Marrow", "r1":11, "c2":"B"},
            {"c1":"Bone Marrow", "r1":12, "c2":"A"},
            {"c1":"Bone Marrow", "r1":13, "c2":"B"},
            {"c1":"Bone Marrow", "r1":14, "c2":"A"},
            {"c1":"Liver", "r1":100, "c2":"A"},
            {"c1":"Liver", "r1":101, "c2":"B"},
            {"c1":"Tonsil", "r1":1000, "c2":"A"},
            {"c1":"Tonsil", "r1":1001, "c2":"B"},
            {"c1":"Fetal Liver", "r1":10000, "c2":"A"},
        ];
        const options = {
            rowAggregate: [{
                field: "c1",
                type: "nominal",
                oneOf: ["Blood"]
            }]
        };
        const selectedRows = selectRows(rowInfo, options);
        expect(selectedRows).toEqual([5,6,7,8,9,10,11,12,13,14,[0,1,2,3,4]]);

        const aggregatedRowInfos = getAggregatedRowInfo(rowInfo, options.rowAggregate).map(d => d[1]);
        const aggregatedRowInfo = aggregatedRowInfos.find(d => Array.isArray(d));   // "Blood"

        let aggregatedValue = getAggregatedValue(aggregatedRowInfo, "c1", "nominal", "sum");
        expect(aggregatedValue).toEqual("Blood");
        
        aggregatedValue = getAggregatedValue(aggregatedRowInfo, "c2", "nominal", "sum");
        expect(aggregatedValue).toEqual("A, B");
        
        aggregatedValue = getAggregatedValue(aggregatedRowInfo, "c2", "nominal", "max");
        expect(aggregatedValue).toEqual("A");
        
        aggregatedValue = getAggregatedValue(aggregatedRowInfo, "r1", "quantitative", "max");
        expect(aggregatedValue).toEqual(5);

        aggregatedValue = getAggregatedValue(aggregatedRowInfo, "r1", "quantitative", "mean");
        expect(aggregatedValue).toEqual(3);
    });

    it('Should produce correct selectRows array after aggregating, filtering, and sorting', () => {
        const rowInfo = [
            {"c1":"Blood", "r1":1, "c2":"A"},
            {"c1":"Blood", "r1":2, "c2":"B"},
            {"c1":"Blood", "r1":3, "c2":"A"},
            {"c1":"Blood", "r1":4, "c2":"B"},
            {"c1":"Blood", "r1":5, "c2":"A"},
            {"c1":"Bone Marrow", "r1":10, "c2":"A"},
            {"c1":"Bone Marrow", "r1":11, "c2":"B"},
            {"c1":"Bone Marrow", "r1":12, "c2":"A"},
            {"c1":"Bone Marrow", "r1":13, "c2":"B"},
            {"c1":"Bone Marrow", "r1":14, "c2":"A"},
            {"c1":"Liver", "r1":100, "c2":"A"},
            {"c1":"Liver", "r1":101, "c2":"B"},
            {"c1":"Tonsil", "r1":1000, "c2":"A"},
            {"c1":"Tonsil", "r1":1001, "c2":"B"},
            {"c1":"Fetal Liver", "r1":10000, "c2":"A"},
        ];
        const options = {
            rowInfoAttributes: [{
                field: "c1",
                type: "nominal",
                aggFunction: "sum"
            },
            {
                field: "c2",
                type: "nominal",
                aggFunction: "sum"
            },
            {
                field: "r1",
                type: "quantitative",
                aggFunction: "sum"
            }],
            rowAggregate: [{
                field: "c1",
                type: "nominal",
                oneOf: ["Blood"]
            },
            {
                field: "c2",
                type: "nominal",
                oneOf: ["A"]
            }],
            rowFilter: [{
                field: "c2",
                type: "nominal",
                notOneOf: ["B"]
            }],
            rowSort: [{
                field: "r1",
                order: "descending",
                type: "quantitative"
            }]
        };
        const selectedRows = selectRows(rowInfo, options);
        expect(selectedRows).toEqual([[5,7,9,10,12,14],[0,1,2,3,4]]);

        const aggregatedRowInfos = getAggregatedRowInfo(rowInfo, options.rowAggregate).map(d => d[1]);
        const aggregatedRowInfo = aggregatedRowInfos.find(d => Array.isArray(d)); // "Blood"

        const aggregatedValue = getAggregatedValue(aggregatedRowInfo, "r1", "quantitative", "mean");
        expect(aggregatedValue).toEqual(3);
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