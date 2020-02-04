/* eslint-env node */

import { getTracksIdsFromViewConfig, updateViewConfigOnSelectGenomicInterval } from './utils-viewconf.js';

import hgDemoViewConfig1 from './viewconfigs/horizontal-multivec.json';
import hgDemoViewConfig2 from './viewconfigs/horizontal-multivec-local.json';

describe('Utilities for processing higlass view config objects', () => {
    it('Should find all horizontal-multivec track IDs', () => {
        const trackIds1 = getTracksIdsFromViewConfig(hgDemoViewConfig1);
        expect(trackIds1.length).toEqual(1);
        expect(trackIds1[0].length).toEqual(3);
        expect(trackIds1[0][0]).toEqual("cistrome-view");
        expect(trackIds1[0][1]).toEqual("cistrome-track");
        expect(trackIds1[0][2]).toBeNull();

        const trackIds2 = getTracksIdsFromViewConfig(hgDemoViewConfig2);
        expect(trackIds2.length).toEqual(1);
        expect(trackIds2[0].length).toEqual(3);
        expect(trackIds2[0][0]).toEqual("cistrome-view-local");
        expect(trackIds2[0][1]).toEqual("cistrome-track-local");
        expect(trackIds2[0][2]).toBeNull()
    });

    it('Should update the view config to create a genomic interval selection', () => {
        const newViewConfig = updateViewConfigOnSelectGenomicInterval(hgDemoViewConfig1, "cistrome-view", "cistrome-track");

        expect(newViewConfig.views.length).toEqual(1);
        expect(newViewConfig.views[0].tracks.center.length).toEqual(1);
        expect(newViewConfig.views[0].uid).toEqual("cistrome-view-with-projection");
        expect(newViewConfig.views[0].tracks.center[0].type).toEqual("combined");
        expect(newViewConfig.views[0].tracks.center[0].contents.length).toEqual(2);
        expect(newViewConfig.views[0].tracks.center[0].contents[0].type).toEqual("horizontal-multivec");
        expect(newViewConfig.views[0].tracks.center[0].contents[0].uid).toEqual("cistrome-track");
        expect(newViewConfig.views[0].tracks.center[0].contents[1].type).toEqual("viewport-projection-horizontal");
        expect(newViewConfig.views[0].tracks.center[0].contents[1].uid).toEqual("cistrome-track-projection");
    });

});
