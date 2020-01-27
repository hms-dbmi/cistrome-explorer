/* eslint-env node */

import { validateWrapperOptions, processWrapperOptions } from './options-schema.js';

describe('Utilities for processing wrapper component options', () => {
  it('Should validate options object', () => {
    const valid = validateWrapperOptions({
        rowInfoPosition: "top"
    });
    expect(valid).toBe(false);
  });
});
