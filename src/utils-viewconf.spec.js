/* eslint-env node */

import { 
    getHMTrackIdsFromViewConfig,
    getSiblingProjectionTracksFromViewConfig,
    updateViewConfigOnSelectGenomicInterval
} from './utils-viewconf.js';

import hgDemoViewConfig1 from './viewconfigs/horizontal-multivec-1.json';
import hgDemoViewConfig2 from './viewconfigs/horizontal-multivec-2.json';
import hgDemoViewConfig3 from './viewconfigs/horizontal-multivec-3.json';

describe('Utilities for processing higlass view config objects', () => {
    it('Should find all horizontal-multivec track IDs', () => {
        const trackIds1 = getHMTrackIdsFromViewConfig(hgDemoViewConfig1);
        expect(trackIds1.length).toEqual(1);
        expect(trackIds1[0].length).toEqual(3);
        expect(trackIds1[0][0]).toEqual("cistrome-view-1");
        expect(trackIds1[0][1]).toEqual("cistrome-track-1");
        expect(trackIds1[0][2]).toBeNull();

        const trackIds2 = getHMTrackIdsFromViewConfig(hgDemoViewConfig2);
        expect(trackIds2.length).toEqual(1);
        expect(trackIds2[0].length).toEqual(3);
        expect(trackIds2[0][0]).toEqual("cistrome-view-2");
        expect(trackIds2[0][1]).toEqual("cistrome-track-2");
        expect(trackIds2[0][2]).toBeNull()
    });

    it('Should update the view config to create a genomic interval selection', () => {
        const newViewConfig = updateViewConfigOnSelectGenomicInterval(hgDemoViewConfig1, "cistrome-view-1", "cistrome-track-1");

        expect(newViewConfig.views.length).toEqual(1);
        expect(newViewConfig.views[0].tracks.center.length).toEqual(1);
        expect(newViewConfig.views[0].uid).toEqual("cistrome-view-1-with-projection");
        expect(newViewConfig.views[0].tracks.center[0].type).toEqual("combined");
        expect(newViewConfig.views[0].tracks.center[0].contents.length).toEqual(2);
        expect(newViewConfig.views[0].tracks.center[0].contents[0].type).toEqual("horizontal-multivec");
        expect(newViewConfig.views[0].tracks.center[0].contents[0].uid).toEqual("cistrome-track-1");
        expect(newViewConfig.views[0].tracks.center[0].contents[1].type).toEqual("viewport-projection-horizontal");
        expect(newViewConfig.views[0].tracks.center[0].contents[1].uid).toEqual("cistrome-track-1-projection");
    });

    it('Should find all viewport-projection-horizontal track siblings of a particular horizontal-multivec track', () => {
        const siblingTrackIds = getSiblingProjectionTracksFromViewConfig(hgDemoViewConfig3, "cistrome-track-1");
        expect(siblingTrackIds.length).toEqual(1);
        expect(siblingTrackIds[0]).toEqual(["cistrome-view-1-with-projection", "cistrome-track-1-projection", "cistrome-track-1-combined"]);

        const siblingTrackIdsEmpty = getSiblingProjectionTracksFromViewConfig(hgDemoViewConfig3, "some-unknown-track");
        expect(siblingTrackIdsEmpty.length).toEqual(0);
    });

});
