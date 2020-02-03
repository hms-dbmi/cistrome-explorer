import pkg from './package.json';

import resolve from '@rollup/plugin-node-resolve';
import json from '@rollup/plugin-json';
import commonjs from '@rollup/plugin-commonjs';
import babel from 'rollup-plugin-babel';
import scss from 'rollup-plugin-scss';
import { terser } from 'rollup-plugin-terser';

import { join } from './rollup.utils.js';

// Constants for output files:
const SRC_DIR = 'src';
const BUILD_DIR = 'build-pkg';
const OUTPUT_CSS = 'index.css';

const outputDefaults = {
    name: 'CistromeHGW',
    sourcemap: 'inline',
    globals: {
        higlass: 'hglib'
    }
};

const pkgConfig = {
    input: join(SRC_DIR, 'index.js'),
    output: [
        {
            file: join(BUILD_DIR, 'index.min.js'),
            format: 'iife',
            plugins: [
                terser()
            ],
            ...outputDefaults
        },
        {
            file: join(BUILD_DIR, 'index.umd.js'),
            format: 'umd',
            ...outputDefaults
        },
        {
            file: join(BUILD_DIR, 'index.cjs.js'),
            format: 'cjs',
            ...outputDefaults
        },
        {
            file: join(BUILD_DIR, 'index.esm.js'),
            format: 'es',
            ...outputDefaults
        }
    ],
    plugins: [
        resolve({
            browser: true,
        }),
        scss({
            output: join(BUILD_DIR, OUTPUT_CSS),
        }),
        json(),
        commonjs({
            include: [
              'node_modules/**',
            ]
        }),
        babel({
            runtimeHelpers: true,
            exclude: 'node_modules/**' // only transpile our source code
        })
    ],
    external: Object.keys(pkg.peerDependencies)
};

export default pkgConfig;