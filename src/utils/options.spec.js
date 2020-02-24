/* eslint-env node */

import { 
    validateWrapperOptions, 
    processWrapperOptions,
    updateWrapperOptions,
    getTrackWrapperOptions,
    DEFAULT_OPTIONS_KEY 
} from './options.js';

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

    it('Should add sorting options to global default', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: DEFAULT_OPTIONS_KEY,
                trackId: DEFAULT_OPTIONS_KEY,
                colToolsPosition: "bottom"
            },
            {
                viewId: "viewA",
                trackId: DEFAULT_OPTIONS_KEY,
                colToolsPosition: "top"
            }
        ]);
        const updatedOptions = updateWrapperOptions(
            processedOptions,
            [{
                field: "groupA",
                type: "nominal",
                order: "ascending"
            }],
            "rowSort",
            DEFAULT_OPTIONS_KEY,
            DEFAULT_OPTIONS_KEY,
            { isReplace: true }
        );
        const globalDefaultOptions = updatedOptions[DEFAULT_OPTIONS_KEY];
        expect(globalDefaultOptions.rowSort.length).toBe(1);
        expect(globalDefaultOptions.rowSort[0].field).toBe("groupA");
        expect(globalDefaultOptions.rowSort[0].type).toBe("nominal");
        expect(globalDefaultOptions.rowSort[0].order).toBe("ascending");
    });

    it('Should update sorting options in global default', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: DEFAULT_OPTIONS_KEY,
                trackId: DEFAULT_OPTIONS_KEY,
                colToolsPosition: "bottom",
                rowSort: [{
                    field: "groupB",
                    type: "quantitative",
                    order: "descending"
                },
                {
                    field: "groupC",
                    type: "quantitative",
                    order: "descending"
                }]
            },
            {
                viewId: "viewA",
                trackId: DEFAULT_OPTIONS_KEY,
                colToolsPosition: "top"
            }
        ]);
        const updatedOptions = updateWrapperOptions(
            processedOptions,
            [{
                field: "groupA",
                type: "nominal",
                order: "ascending"
            }],
            "rowSort",
            DEFAULT_OPTIONS_KEY,
            DEFAULT_OPTIONS_KEY,
            { isReplace: true }
        );
        const globalDefaultOptions = updatedOptions[DEFAULT_OPTIONS_KEY];
        expect(globalDefaultOptions.rowSort.length).toBe(1);
        expect(globalDefaultOptions.rowSort[0].field).toBe("groupA");
        expect(globalDefaultOptions.rowSort[0].type).toBe("nominal");
        expect(globalDefaultOptions.rowSort[0].order).toBe("ascending");
    });

    it('Should add filtering options in the list to track-global default', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: DEFAULT_OPTIONS_KEY,
                trackId: DEFAULT_OPTIONS_KEY,
                colToolsPosition: "bottom",
                rowSort: [{
                    field: "groupA",
                    type: "quantitative",
                    order: "descending"
                }],
                rowFilter: [{
                    field: "groupB",
                    type: "quantitative",
                    contains: "substringB"
                }]
            },
            {
                viewId: "viewA",
                trackId: DEFAULT_OPTIONS_KEY,
                colToolsPosition: "top"
            }
        ]);
        const updatedOptions = updateWrapperOptions(
            processedOptions,
            {
                field: "groupC",
                type: "nominal",
                contains: "substringC"
            },
            "rowFilter",
            "viewA",
            DEFAULT_OPTIONS_KEY,
            { isReplace: false }
        );
        const trackGlobalOptions = updatedOptions["viewA"][DEFAULT_OPTIONS_KEY];
        expect(trackGlobalOptions.rowFilter.length).toBe(2);
        expect(trackGlobalOptions.rowFilter.filter(d => d.field === "groupC").length).toBe(1);
        expect(trackGlobalOptions.rowFilter.filter(d => d.type === "nominal").length).toBe(1);
        expect(trackGlobalOptions.rowFilter.filter(d => d.contains === "substringC").length).toBe(1);
    });

    it('Should add sorting info to non-global options', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: "viewA",
                trackId: "trackA",
                colToolsPosition: "top"
            }
        ]);
        const updatedOptions = updateWrapperOptions(
            processedOptions,
            [{
                field: "groupA",
                type: "nominal",
                order: "ascending"
            }],
            "rowSort",
            "viewA",
            "trackA",
            { isReplace: true }
        );
        const nonGlobalOptions = updatedOptions["viewA"]["trackA"];
        expect(nonGlobalOptions !== undefined).toBe(true);
        expect(nonGlobalOptions.rowSort.length).toBe(1);
        expect(nonGlobalOptions.rowSort[0].field).toBe("groupA");
        expect(nonGlobalOptions.rowSort[0].type).toBe("nominal");
        expect(nonGlobalOptions.rowSort[0].order).toBe("ascending");
    });

    it('Should return the processed options object for a particular track', () => {
        const trackOptionsAA = getTrackWrapperOptions({
            viewA: { trackA: { colToolsPosition: "top" } },
            [DEFAULT_OPTIONS_KEY]: { colToolsPosition: "bottom" }
        }, "viewA", "trackA");
        expect(trackOptionsAA.colToolsPosition).toEqual("top");

        const trackOptionsAB = getTrackWrapperOptions({
            viewA: { trackA: { colToolsPosition: "top" } },
            [DEFAULT_OPTIONS_KEY]: { colToolsPosition: "bottom" }
        }, "viewA", "trackB");
        expect(trackOptionsAB.colToolsPosition).toEqual("bottom");

        const trackOptionsBA = getTrackWrapperOptions({
            viewA: { trackA: { colToolsPosition: "top" } },
            [DEFAULT_OPTIONS_KEY]: { colToolsPosition: "bottom" }
        }, "viewB", "trackA");
        expect(trackOptionsBA.colToolsPosition).toEqual("bottom");

        const trackOptionsCA = getTrackWrapperOptions({
            viewA: { trackA: { colToolsPosition: "top" } },
            viewC: { [DEFAULT_OPTIONS_KEY]: { colToolsPosition: "hidden" } },
            [DEFAULT_OPTIONS_KEY]: { colToolsPosition: "bottom" }
        }, "viewC", "trackA");
        expect(trackOptionsCA.colToolsPosition).toEqual("hidden");
    });
});