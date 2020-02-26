module.exports = {
    verbose: true,
    moduleFileExtensions: [ "js" ],
    transformIgnorePatterns: [
        "/node_modules/(?!higlass-register|higlass-multivec).+\\.js$"
    ],
    transform: {
        "^.+\\.js$": "babel-jest"
    },
    moduleNameMapper: {
        "\\.(css|scss)$": "identity-obj-proxy"
    },
};