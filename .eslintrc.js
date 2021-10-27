module.exports = {
    "env": {
        "browser": true,
        "es2021": true,
        jest: true
    },
    "extends": [
        "eslint:recommended",
        "plugin:react/recommended"
    ],
    "parserOptions": {
        "ecmaFeatures": {
            "jsx": true
        },
        "ecmaVersion": 13,
        "sourceType": "module"
    },
    "plugins": [
        "react",
        "jest"
    ],
    settings: {
        react: {
            version: "detect"
        }
    },
    "rules": {
        "indent": [
            "error",
            4
        ],
        "react/prop-types": 0,
        "no-unused-vars": "warn",
        "no-prototype-builtins": "warn",
        "no-constant-condition": 0,
        "react/display-name": 0,
        "linebreak-style": [
            "error",
            "unix"
        ],
        "quotes": [
            "error",
            "double"
        ],
        "semi": [
            "error",
            "always"
        ]
    }
};
