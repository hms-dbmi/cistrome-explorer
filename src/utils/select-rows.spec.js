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
});
