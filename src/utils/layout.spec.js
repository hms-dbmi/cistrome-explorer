/* eslint-env node */

import { 
    getNumOfTracks
} from './layout.js';

describe('Functions for calculating the number of tracks that should be visualized', () => {
    it('Should correctly calculate the number of tracks to visualize', () => {
        let isLeft = false;
        let resolveYScales = [];
        let numTracks = getNumOfTracks(isLeft, resolveYScales);
        expect(numTracks).toEqual(0);

        isLeft = false;
        resolveYScales = [false, false, false, false, false];
        numTracks = getNumOfTracks(isLeft, resolveYScales);
        expect(numTracks).toEqual(5);

        isLeft = true;
        resolveYScales = [true];
        numTracks = getNumOfTracks(isLeft, resolveYScales);
        expect(numTracks).toEqual(2);

        isLeft = false;
        resolveYScales = [true, true, true];
        numTracks = getNumOfTracks(isLeft, resolveYScales);
        expect(numTracks).toEqual(6);

        isLeft = false;
        resolveYScales = [true, false, true];
        numTracks = getNumOfTracks(isLeft, resolveYScales);
        expect(numTracks).toEqual(6);

        isLeft = false;
        resolveYScales = [false, false, true];
        numTracks = getNumOfTracks(isLeft, resolveYScales);
        expect(numTracks).toEqual(4);
    });
});
