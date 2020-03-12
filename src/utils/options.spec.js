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
                viewId: 12
            }
        );
        expect(valid).toBe(false);
    });

    it('Should validate options object when correct', () => {
        const valid = validateWrapperOptions(
            {
                viewId: "uid"
            }
        );
        expect(valid).toBe(true);
    });

    it('Should validate options array when incorrect', () => {
        const valid = validateWrapperOptions([
            {
                viewId: ["uid1", "uid2"]
            }
        ]);
        expect(valid).toBe(false);
    });

    it('Should validate options array when correct', () => {
        const valid = validateWrapperOptions([
            {
                viewId: "viewA",
                trackId: "trackA"
            }
        ]);
        expect(valid).toBe(true);
    });

    it('Should process options object', () => {
        const processedOptions = processWrapperOptions(
            {
                rowHighlight: {field: "fieldA", type: "nominal", contains: "abc"}
            }
        );

        expect(Object.keys(processedOptions)).toEqual(['default']);
        expect(processedOptions.default.rowHighlight.field).toEqual("fieldA");
    });

    it('Should process options array with global default', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: "default",
                trackId: "default",
                rowHighlight: {field: "fieldA", type: "nominal", contains: "abc"}
            }
        ]);

        expect(Object.keys(processedOptions)).toEqual(['default']);
        expect(processedOptions.default.rowHighlight.field).toEqual("fieldA");
    });

    it('Should process options array with view default', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: "default",
                trackId: "default",
                rowHighlight: {field: "fieldA", type: "nominal", contains: "abc"}
            },
            {
                viewId: "viewA",
                trackId: "default",
                rowHighlight: {field: "fieldB", type: "nominal", contains: "abc"}
            }
        ]);

        expect(Object.keys(processedOptions)).toEqual(['default', 'viewA']);
        expect(processedOptions.default.rowHighlight.field).toEqual("fieldA");
        expect(processedOptions.viewA.default.rowHighlight.field).toEqual("fieldB");
    });

    it('Should process options array with view and track defaults', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: "default",
                trackId: "default",
                rowHighlight: {field: "fieldDD", type: "nominal", contains: "abc"}
            },
            {
                viewId: "viewA",
                trackId: "default",
                rowHighlight: {field: "fieldAD", type: "nominal", contains: "abc"}
            },
            {
                viewId: "viewA",
                trackId: "trackA",
                rowHighlight: {field: "fieldAA", type: "nominal", contains: "abc"}
            },
            {
                viewId: "viewA",
                trackId: "trackB",
                rowHighlight: {field: "fieldAB", type: "nominal", contains: "abc"}
            },
            {
                viewId: "viewB",
                trackId: "default",
                rowHighlight: {field: "fieldBD", type: "nominal", contains: "abc"}
            },
        ]);

        expect(Object.keys(processedOptions)).toEqual(['default', 'viewA', 'viewB']);
        expect(processedOptions.default.rowHighlight.field).toEqual("fieldDD");
        expect(processedOptions.viewA.default.rowHighlight.field).toEqual("fieldAD");
        expect(processedOptions.viewA.trackA.rowHighlight.field).toEqual("fieldAA");
        expect(processedOptions.viewA.trackB.rowHighlight.field).toEqual("fieldAB");
        expect(processedOptions.viewB.default.rowHighlight.field).toEqual("fieldBD");
    });

    it('Should add sorting options to global default', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: DEFAULT_OPTIONS_KEY,
                trackId: DEFAULT_OPTIONS_KEY,
                rowHighlight: {field: "fieldDD", type: "nominal", contains: "abc"}
            },
            {
                viewId: "viewA",
                trackId: DEFAULT_OPTIONS_KEY,
                rowHighlight: {field: "fieldAD", type: "nominal", contains: "abc"}
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

    it('Should add sorting options to global default when there is no track-specific options', () => {
        const processedOptions = processWrapperOptions([
            {
                viewId: DEFAULT_OPTIONS_KEY,
                trackId: DEFAULT_OPTIONS_KEY,
                rowHighlight: {field: "fieldDD", type: "nominal", contains: "abc"}
            },
            {
                viewId: "viewA",
                trackId: DEFAULT_OPTIONS_KEY,
                rowHighlight: {field: "fieldAD", type: "nominal", contains: "abc"}
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
                rowHighlight: {field: "fieldAD", type: "nominal", contains: "abc"}
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
                rowHighlight: {field: "fieldAD", type: "nominal", contains: "abc"}
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
                rowHighlight: {field: "fieldAA", type: "nominal", contains: "abc"}
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
            viewA: { trackA: { rowHighlight: {field: "fieldAA", type: "nominal", contains: "abc"} } },
            [DEFAULT_OPTIONS_KEY]: { rowHighlight: {field: "fieldDD", type: "nominal", contains: "abc"} }
        }, "viewA", "trackA");
        expect(trackOptionsAA.rowHighlight.field).toEqual("fieldAA");

        const trackOptionsAB = getTrackWrapperOptions({
            viewA: { trackA: { rowHighlight: {field: "fieldAA", type: "nominal", contains: "abc"} } },
            [DEFAULT_OPTIONS_KEY]: { rowHighlight: {field: "fieldDD", type: "nominal", contains: "abc"} }
        }, "viewA", "trackB");
        expect(trackOptionsAB.rowHighlight.field).toEqual("fieldDD");

        const trackOptionsBA = getTrackWrapperOptions({
            viewA: { trackA: { rowHighlight: {field: "fieldAA", type: "nominal", contains: "abc"} } },
            [DEFAULT_OPTIONS_KEY]: { rowHighlight: {field: "fieldDD", type: "nominal", contains: "abc"} }
        }, "viewB", "trackA");
        expect(trackOptionsBA.rowHighlight.field).toEqual("fieldDD");

        const trackOptionsCA = getTrackWrapperOptions({
            viewA: { trackA: { rowHighlight: {field: "fieldAA", type: "nominal", contains: "abc"} } },
            viewC: { [DEFAULT_OPTIONS_KEY]: { rowHighlight: {field: "fieldCD", type: "nominal", contains: "abc"} } },
            [DEFAULT_OPTIONS_KEY]: { rowHighlight: {field: "fieldDD", type: "nominal", contains: "abc"} }
        }, "viewC", "trackA");
        expect(trackOptionsCA.rowHighlight.field).toEqual("fieldCD");
    });
});