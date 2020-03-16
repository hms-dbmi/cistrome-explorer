import pkg from './package.json';
import merge from 'lodash/merge';

import resolve from '@rollup/plugin-node-resolve';
import json from '@rollup/plugin-json';
import commonjs from '@rollup/plugin-commonjs';
import html from '@rollup/plugin-html';
import babel from 'rollup-plugin-babel';
import scss from 'rollup-plugin-scss';
import visualizer from 'rollup-plugin-visualizer';
import { terser } from 'rollup-plugin-terser';

import { join } from './rollup.utils.js';
import { htmlFromTemplate } from './rollup.demo.utils.js';

// Constants for output files:
const SRC_DIR = 'src';
const BUILD_DIR = 'build-demo';
const OUTPUT_JS = {
    'production': 'index.min.js',
    'development': 'index.js',
};
const OUTPUT_CSS = 'index.css';
const OUTPUT_HTML = 'index.html';

// Only import dev server if necessary:
let serve = (() => {});
let livereload = (() => {});
if(process.env.NODE_ENV === 'development') {
    serve = require('rollup-plugin-serve');
    livereload = require('rollup-plugin-livereload');
}

// The base rollup configuration. To be merged with dev or prod object.
const baseConfig = {
    input: join(SRC_DIR, 'demo', 'index.js'),
    output: {
        file: join(BUILD_DIR, OUTPUT_JS[process.env.NODE_ENV]),
        sourcemap: true,
        globals: {
            higlass: 'hglib'
        }
    },
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
        }),
        html({
            title: pkg.name,
            publicPath: '/',
            fileName: OUTPUT_HTML,
            template: ({ publicPath, title }) => {
                return htmlFromTemplate({
                    publicPath: (process.env.NODE_ENV === 'production' ? publicPath : './'),
                    title: title,
                    nodeEnv: process.env.NODE_ENV,
                    cssFile: OUTPUT_CSS,
                    jsFile: OUTPUT_JS[process.env.NODE_ENV],
                });
            }
        })
    ],
    external: Object.keys(pkg.peerDependencies)
};

const devConfig = {
    output: {
        format: 'umd',
    },
    plugins: [
        ...baseConfig.plugins.map(() => {}),
        visualizer({
            filename: join(BUILD_DIR, 'stats.html')
        }),
        serve({
            port: 8000,
            contentBase: BUILD_DIR
        }),
        livereload(BUILD_DIR)
    ]
};

const prodConfig = {
    output: {
        format: 'iife',
        plugins: [
            terser()
        ],
    }
};

export default merge(
    (process.env.NODE_ENV === 'development' ? devConfig : prodConfig), 
    baseConfig
);