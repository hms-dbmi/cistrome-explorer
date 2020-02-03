/* eslint-env node */

import { getTracksIdsFromViewConfig } from './utils-viewconf.js';

import hgDemoViewConfig1 from './configs/horizontal-multivec.json';
import hgDemoViewConfig2 from './configs/horizontal-multivec-local.json';

describe('Utilities for processing higlass view config objects', () => {
    it('Should find all horizontal-multivec track IDs', () => {
        const trackIds1 = getTracksIdsFromViewConfig(hgDemoViewConfig1);
        expect(trackIds1.length).toEqual(1);
        expect(trackIds1[0].length).toEqual(2);
        expect(trackIds1[0][0]).toEqual("UiHlCoxRQ-aITBDi5j8b_w");
        expect(trackIds1[0][1]).toEqual("cistrome-track");

        const trackIds2 = getTracksIdsFromViewConfig(hgDemoViewConfig2);
        expect(trackIds2.length).toEqual(1);
        expect(trackIds2[0].length).toEqual(2);
        expect(trackIds2[0][0]).toEqual("UiHlCoxRQ-aITBDi5j8b_l");
        expect(trackIds2[0][1]).toEqual("cistrome-track-local");
    });

});
