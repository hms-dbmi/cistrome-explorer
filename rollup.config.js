import pkg from './package.json';
import merge from 'lodash/merge';

import resolve from '@rollup/plugin-node-resolve';
import json from '@rollup/plugin-json';
import commonjs from '@rollup/plugin-commonjs';
import babel from 'rollup-plugin-babel';
import scss from 'rollup-plugin-scss';
import visualizer from 'rollup-plugin-visualizer';
import { terser } from 'rollup-plugin-terser';

const BUILD_DIR = 'dist';

let serve = (() => {});
let livereload = (() => {});
if(process.env.NODE_ENV === 'development') {
    serve = require('rollup-plugin-serve');
    livereload = require('rollup-plugin-livereload');
}

const defaultOutputConfig = {
    sourcemap: 'inline',
    globals: {
        higlass: 'hglib'
    }
};

const defaultConfig = {
    input: 'src/index.js',
    output: [],
    plugins: [
        resolve({
            browser: true,
        }),
        scss({
            output: `${BUILD_DIR}/index.css`,
        }),
        json(),
        commonjs({
            include: [
              'node_modules/**',
            ],
            namedExports: {
              'node_modules/react/index.js': ['useState', 'useEffect', 'useRef', 'cloneElement'],
              'node_modules/higlass/dist/hglib.js': ['HiGlassComponent']
            },
        }),
        babel({
            runtimeHelpers: true,
            exclude: 'node_modules/**' // only transpile our source code
        }),
    ],
    external: Object.keys(pkg.peerDependencies)
};

const devConfig = {
    output: [
        {
            ...defaultOutputConfig,
            file: `${BUILD_DIR}/index.js`,
            format: 'umd',
        },
        {
            ...defaultOutputConfig,
            file: `${BUILD_DIR}/index.esm.js`,
            format: 'es',
        }
    ],
    plugins: [
        ...defaultConfig.plugins.map(() => {}),
        ...[
            visualizer({
                filename: `${BUILD_DIR}/stats.html`
            }),
            serve({
                port: 8000,
                contentBase: BUILD_DIR
            }),
            livereload(BUILD_DIR)
        ]
    ]
};

const prodConfig = {
    output: [
        {
            ...defaultOutputConfig,
            file: `${BUILD_DIR}/index.min.js`,
            format: 'iife',
            plugins: [
                terser()
            ],
        }
    ],
    plugins: []
};

export default merge(
    (process.env.NODE_ENV === 'development' ? devConfig : prodConfig), 
    defaultConfig
);