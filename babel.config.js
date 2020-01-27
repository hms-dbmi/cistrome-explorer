module.exports = {
    env: {
        test: {
            presets: [
                ["@babel/preset-env", { modules: "commonjs" }],
                "@babel/preset-react"
            ]
        },
        production: {
            presets: [
                ["@babel/preset-env", { modules: false }],
                "@babel/preset-react"
            ]
        },
        development: {
            presets: [
                ["@babel/preset-env", { modules: false }],
                "@babel/preset-react"
            ]
        }
    }
};