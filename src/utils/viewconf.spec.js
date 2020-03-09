/* eslint-env node */
import cloneDeep from 'lodash/cloneDeep';

import { 
    getHMTrackIdsFromViewConfig,
    getSiblingVPHTrackIdsFromViewConfig,
    updateViewConfigOnSelectGenomicInterval,
    updateViewConfigOnSelectRowsByTrack,
    getHMSelectedRowsFromViewConfig,
    getAllViewAndTrackPairs
} from './viewconf.js';

import hgDemoViewConfig1 from '../viewconfigs/horizontal-multivec-1.json';
import hgDemoViewConfig2 from '../viewconfigs/horizontal-multivec-2.json';
import hgDemoViewConfig3 from '../viewconfigs/horizontal-multivec-3.json';
import hgDemoViewConfig4 from '../viewconfigs/horizontal-multivec-4.json';
import hgDemoViewConfig5 from '../viewconfigs/horizontal-multivec-5.json';

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
        expect(newViewConfig.views[0].tracks.whole.length).toEqual(1);
        expect(newViewConfig.views[0].tracks.whole[0].type).toEqual("viewport-projection-horizontal");
        expect(newViewConfig.views[0].tracks.whole[0].fromViewUid).toEqual(null);
    });

    it('Should find all viewport-projection-horizontal track siblings of a particular horizontal-multivec track', () => {
        const siblingTrackIds = getSiblingVPHTrackIdsFromViewConfig(hgDemoViewConfig3, "cistrome-view-1");
        expect(siblingTrackIds.length).toEqual(1);
        expect(siblingTrackIds[0].viewId).toEqual("cistrome-view-1");
        expect(siblingTrackIds[0].trackId).toEqual("cistrome-view-1-projection-track-1");

        const siblingTrackIdsEmpty = getSiblingVPHTrackIdsFromViewConfig(hgDemoViewConfig3, "some-unknown-view");
        expect(siblingTrackIdsEmpty.length).toEqual(0);
    });

    it('Should set the selectRows option value of a particular horizontal-multivec track', () => {
        const currViewConfig = cloneDeep(hgDemoViewConfig2);
        const newViewConfig = updateViewConfigOnSelectRowsByTrack(currViewConfig, [1, 2, 3], "cistrome-view-2", "cistrome-track-2");
        expect(newViewConfig.views[0].tracks.center[0].options.selectRows).toEqual([1, 2, 3]);
    });

    it('Should get the selectRows option value of a particular horizontal-multivec track', () => {
        const selectedRows = getHMSelectedRowsFromViewConfig(hgDemoViewConfig5, "cistrome-view-5", "cistrome-track-5");
        expect(selectedRows).toEqual([5, 3, 1, 4]);
    });
    
    it('Should get viewId and trackId pairs of horizontal-multivec tracks', () => {
        const searchedPairs = getAllViewAndTrackPairs(hgDemoViewConfig4);
        expect(searchedPairs.length).toEqual(2);
        expect(searchedPairs.find(d => d.viewId === "cistrome-view-4-1" && d.trackId === "cistrome-track-4-1") !== undefined).toEqual(true);
    });
});
