/* eslint-env node */

import { validateWrapperOptions, processWrapperOptions } from './options.js';

describe('Utilities for processing wrapper component options', () => {
    it('Should validate options object when incorrect', () => {
        const valid = validateWrapperOptions(
            {
                colToolsPosition: "left"
            }
        );
        expect(valid).toBe(false);
    });

    it('Should validate options object when correct', () => {
        const valid = validateWrapperOptions(
            {
                colToolsPosition: "bottom"
            }
        );
        expect(valid).toBe(true);
    });

    it('Should validate options array when incorrect', () => {
        const valid = validateWrapperOptions([
            {
                colToolsPosition: "left"
            }
        ]);
        expect(valid).toBe(false);
    });

    it('Should validate options array when correct', () => {
        const valid = validateWrapperOptions([
            {
                viewId: "viewA",
                trackId: "trackA",
                colToolsPosition: "bottom"
            }
        ]);
        expect(valid).toBe(true);
    });

    it('Should process options object', () => {
        const processedOptions = processWrapperOptions(
            {
                colToolsPosition: "bottom"
            }
        );

        expect(Object.keys(processedOptions)).toEqual(['default']);
        expect(processedOptions.default.colToolsPosition).toEqual("bottom");
    });

    it('Should process options array with global default', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: "default",
                trackId: "default",
                colToolsPosition: "bottom"
            }
        ]);

        expect(Object.keys(processedOptions)).toEqual(['default']);
        expect(processedOptions.default.colToolsPosition).toEqual("bottom");
    });

    it('Should process options array with view default', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: "default",
                trackId: "default",
                colToolsPosition: "bottom"
            },
            {
                viewId: "viewA",
                trackId: "default",
                colToolsPosition: "top"
            }
        ]);

        expect(Object.keys(processedOptions)).toEqual(['default', 'viewA']);
        expect(processedOptions.default.colToolsPosition).toEqual("bottom");
        expect(processedOptions.viewA.default.colToolsPosition).toEqual("top");
    });

    it('Should process options array with view and track defaults', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: "default",
                trackId: "default",
                colToolsPosition: "top"
            },
            {
                viewId: "viewA",
                trackId: "default",
                colToolsPosition: "bottom"
            },
            {
                viewId: "viewA",
                trackId: "trackA",
                colToolsPosition: "hidden"
            },
            {
                viewId: "viewA",
                trackId: "trackB",
                colToolsPosition: "top"
            },
            {
                viewId: "viewB",
                trackId: "default",
                colToolsPosition: "hidden"
            },
        ]);

        expect(Object.keys(processedOptions)).toEqual(['default', 'viewA', 'viewB']);
        expect(processedOptions.default.colToolsPosition).toEqual("top");
        expect(processedOptions.viewA.default.colToolsPosition).toEqual("bottom");
        expect(processedOptions.viewA.trackA.colToolsPosition).toEqual("hidden");
        expect(processedOptions.viewA.trackB.colToolsPosition).toEqual("top");
        expect(processedOptions.viewB.default.colToolsPosition).toEqual("hidden");
    });
});
