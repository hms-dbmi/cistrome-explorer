/* eslint-env node */

import { 
    getHMTrackIdsFromViewConfig,
    getSiblingVPHTrackIdsFromViewConfig,
    updateViewConfigOnSelectGenomicInterval
} from './viewconf.js';

import hgDemoViewConfig1 from '../viewconfigs/horizontal-multivec-1.json';
import hgDemoViewConfig2 from '../viewconfigs/horizontal-multivec-2.json';
import hgDemoViewConfig3 from '../viewconfigs/horizontal-multivec-3.json';

describe('Utilities for processing higlass view config objects', () => {
    it('Should find all horizontal-multivec track IDs', () => {
        const trackIds1 = getHMTrackIdsFromViewConfig(hgDemoViewConfig1);
        expect(trackIds1.length).toEqual(1);
        expect(trackIds1[0].viewId).toEqual("cistrome-view-1");
        expect(trackIds1[0].trackId).toEqual("cistrome-track-1");

        const trackIds2 = getHMTrackIdsFromViewConfig(hgDemoViewConfig2);
        expect(trackIds2.length).toEqual(1);
        expect(trackIds2[0].viewId).toEqual("cistrome-view-2");
        expect(trackIds2[0].trackId).toEqual("cistrome-track-2");
    });

    it('Should update the view config to create a genomic interval selection', () => {
        const newViewConfig = updateViewConfigOnSelectGenomicInterval(hgDemoViewConfig1, "cistrome-view-1", "cistrome-track-1");

        expect(newViewConfig.views.length).toEqual(1);
        expect(newViewConfig.views[0].tracks.center.length).toEqual(1);
        expect(newViewConfig.views[0].uid).toEqual("cistrome-view-1-with-col-projection");
        expect(newViewConfig.views[0].tracks.center[0].type).toEqual("combined");
        expect(newViewConfig.views[0].tracks.center[0].contents.length).toEqual(2);
        expect(newViewConfig.views[0].tracks.center[0].contents[0].type).toEqual("horizontal-multivec");
        expect(newViewConfig.views[0].tracks.center[0].contents[0].uid).toEqual("cistrome-track-1");
        expect(newViewConfig.views[0].tracks.center[0].contents[1].type).toEqual("viewport-projection-horizontal");
        expect(newViewConfig.views[0].tracks.center[0].contents[1].uid).toEqual("cistrome-track-1-col-projection");
    });

    it('Should find all viewport-projection-horizontal track siblings of a particular horizontal-multivec track', () => {
        const siblingTrackIds = getSiblingVPHTrackIdsFromViewConfig(hgDemoViewConfig3, "cistrome-track-1");
        expect(siblingTrackIds.length).toEqual(1);
        expect(siblingTrackIds[0].viewId).toEqual("cistrome-view-1-with-col-projection");
        expect(siblingTrackIds[0].trackId).toEqual("cistrome-track-1-col-projection");
        expect(siblingTrackIds[0].combinedTrackId).toEqual("cistrome-track-1-combined");

        const siblingTrackIdsEmpty = getSiblingVPHTrackIdsFromViewConfig(hgDemoViewConfig3, "some-unknown-track");
        expect(siblingTrackIdsEmpty.length).toEqual(0);
    });

});
