import resolve from '@rollup/plugin-node-resolve';
import json from '@rollup/plugin-json';
import commonjs from '@rollup/plugin-commonjs';
import { babel } from '@rollup/plugin-babel';
import scss from 'rollup-plugin-scss';
import { terser } from 'rollup-plugin-terser';

import { join } from 'path';
import React from 'react';

import pkg from './package.json';


// Constants for output files:
const SRC_DIR = 'src';
const BUILD_DIR = 'build-pkg';
const OUTPUT_CSS = 'index.css';

const outputDefaults = {
    name: 'HiGlassMeta',
    // We want sourcemap files to be created for debugging purposes.
    // https://rollupjs.org/guide/en/#outputsourcemap
    sourcemap: true,
    // Since we want HiGlass, React, and ReactDOM to be externals,
    // we need to tell the bundle how these libraries can be found as global variables.
    // Reference: https://rollupjs.org/guide/en/#outputglobals
    globals: {
        higlass: 'hglib',
        react: 'React',
        'react-dom': 'ReactDOM',
    }
};

const pkgConfig = {
    input: join(SRC_DIR, 'index.js'),
    output: [
        {
            file: join(BUILD_DIR, 'index.min.js'),
            // Reference: https://rollupjs.org/guide/en/#outputformat
            format: 'iife',
            plugins: [
                terser()
            ],
            ...outputDefaults
        },
        {
            file: join(BUILD_DIR, 'index.umd.js'),
            // Reference: https://rollupjs.org/guide/en/#outputformat
            format: 'umd',
            ...outputDefaults
        },
        {
            file: join(BUILD_DIR, 'index.cjs.js'),
            // Reference: https://rollupjs.org/guide/en/#outputformat
            format: 'cjs',
            ...outputDefaults
        },
        {
            file: join(BUILD_DIR, 'index.esm.js'),
            // Reference: https://rollupjs.org/guide/en/#outputformat
            format: 'es',
            ...outputDefaults
        }
    ],
    plugins: [
        // Tell Rollup how to resolve packages in node_modules.
        // Reference: https://github.com/rollup/plugins/tree/master/packages/commonjs#using-with-rollupplugin-node-resolve
        resolve({
            browser: true,
        }),
        // Tell Rollup how to handle CSS and SCSS imports.
        scss({
            output: join(BUILD_DIR, OUTPUT_CSS),
        }),
        // Tell Rollup how to handle JSON imports.
        json(),
        // Need to convert CommonJS modules in node_modules to ES6.
        // Reference: https://github.com/rollup/plugins/tree/master/packages/node-resolve#using-with-rollupplugin-commonjs
        commonjs({
            // Using this RegEx rather than 'node_modules/**' is suggested, to enable symlinks.
            // Reference: https://github.com/rollup/plugins/tree/master/packages/commonjs#usage-with-symlinks
            include: /node_modules/,
            namedExports: {
                // Need to explicitly tell Rollup how to handle imports like `React, { useState }`
                // Reference: https://github.com/rollup/rollup-plugin-commonjs/issues/407#issuecomment-527837831
                // Reference: https://github.com/facebook/react/issues/11503
                'node_modules/react/index.js': Object.keys(React)
            }
        }),
        // Tell Rollup to compile our source files with Babel.
        // Note: This plugin respects Babel config files by default.
        // Reference: https://github.com/rollup/plugins/tree/master/packages/babel
        babel({
            // The 'runtime' option is recommended when bundling libraries.
            // Reference: https://github.com/rollup/plugins/tree/master/packages/babel#babelhelpers
            babelHelpers: 'runtime',
            // Reference: https://github.com/rollup/plugins/issues/381#issuecomment-627215009
            skipPreflightCheck: true,
            // Only transpile our source code.
            // Reference: https://github.com/rollup/plugins/tree/master/packages/babel#extensions
            exclude: 'node_modules/**'
        })
    ],
    // We do not to include HiGlass, React, or ReactDOM in the bundle.
    external: Object.keys(pkg.peerDependencies)
};

export default pkgConfig;