/* eslint-env node */

import { validateWrapperOptions, processWrapperOptions } from './utils-options.js';

describe('Utilities for processing wrapper component options', () => {
    it('Should validate options object when incorrect', () => {
        const valid = validateWrapperOptions(
            {
                rowInfoPosition: "top"
            }
        );
        expect(valid).toBe(false);
    });

    it('Should validate options object when correct', () => {
        const valid = validateWrapperOptions(
            {
                rowInfoPosition: "right"
            }
        );
        expect(valid).toBe(true);
    });

    it('Should validate options array when incorrect', () => {
        const valid = validateWrapperOptions([
            {
                rowInfoPosition: "top"
            }
        ]);
        expect(valid).toBe(false);
    });

    it('Should validate options array when correct', () => {
        const valid = validateWrapperOptions([
            {
                viewId: "viewA",
                trackId: "trackA",
                rowInfoPosition: "left"
            }
        ]);
        expect(valid).toBe(true);
    });

    it('Should process options object', () => {
        const processedOptions = processWrapperOptions(
            {
                rowInfoPosition: "left"
            }
        );

        expect(Object.keys(processedOptions)).toEqual(['default']);
        expect(processedOptions.default.rowInfoPosition).toEqual("left");
    });

    it('Should process options array with global default', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: "default",
                trackId: "default",
                rowInfoPosition: "left"
            }
        ]);

        expect(Object.keys(processedOptions)).toEqual(['default']);
        expect(processedOptions.default.rowInfoPosition).toEqual("left");
    });

    it('Should process options array with view default', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: "default",
                trackId: "default",
                rowInfoPosition: "left"
            },
            {
                viewId: "viewA",
                trackId: "default",
                rowInfoPosition: "right"
            }
        ]);

        expect(Object.keys(processedOptions)).toEqual(['default', 'viewA']);
        expect(processedOptions.default.rowInfoPosition).toEqual("left");
        expect(processedOptions.viewA.default.rowInfoPosition).toEqual("right");
    });

    it('Should process options array with view and track defaults', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: "default",
                trackId: "default",
                rowInfoPosition: "left"
            },
            {
                viewId: "viewA",
                trackId: "default",
                rowInfoPosition: "right"
            },
            {
                viewId: "viewA",
                trackId: "trackA",
                rowInfoPosition: "hidden"
            },
            {
                viewId: "viewA",
                trackId: "trackB",
                rowInfoPosition: "left"
            },
            {
                viewId: "viewB",
                trackId: "default",
                rowInfoPosition: "hidden"
            },
        ]);

        expect(Object.keys(processedOptions)).toEqual(['default', 'viewA', 'viewB']);
        expect(processedOptions.default.rowInfoPosition).toEqual("left");
        expect(processedOptions.viewA.default.rowInfoPosition).toEqual("right");
        expect(processedOptions.viewA.trackA.rowInfoPosition).toEqual("hidden");
        expect(processedOptions.viewA.trackB.rowInfoPosition).toEqual("left");
        expect(processedOptions.viewB.default.rowInfoPosition).toEqual("hidden");
    });
});
