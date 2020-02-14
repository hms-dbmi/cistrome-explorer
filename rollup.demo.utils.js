const peerScripts = {
    'development': `
        <script crossorigin type="text/javascript" src="https://unpkg.com/react@16/umd/react.development.js"></script>
        <script crossorigin type="text/javascript" src="https://unpkg.com/react-dom@16/umd/react-dom.development.js"></script>
        <script crossorigin type="text/javascript" src="https://unpkg.com/pixi.js@5/dist/pixi.js"></script>
        <script crossorigin type="text/javascript" src="https://unpkg.com/react-bootstrap@0.32.1/dist/react-bootstrap.js"></script>
        <script crossorigin type="text/javascript" src="https://unpkg.com/higlass@1.8.4/dist/hglib.js"></script>
    `,
    'production': `
        <script crossorigin type="text/javascript" src="https://unpkg.com/react@16/umd/react.production.min.js"></script>
        <script crossorigin type="text/javascript" src="https://unpkg.com/react-dom@16/umd/react-dom.production.min.js"></script>
        <script crossorigin type="text/javascript" src="https://unpkg.com/pixi.js@5/dist/pixi.min.js"></script>
        <script crossorigin type="text/javascript" src="https://unpkg.com/react-bootstrap@0.32.1/dist/react-bootstrap.min.js"></script>
        <script crossorigin type="text/javascript" src="https://unpkg.com/higlass@1.8.4/dist/hglib.min.js"></script>
    `
};

export function htmlFromTemplate({ title, publicPath, nodeEnv, cssFile, jsFile }) {
    return `<!DOCTYPE html>
<html>
    <head lang="en">
        <title>${title}</title>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1" />
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
        <link rel="stylesheet" href="https://unpkg.com/higlass@1.8.4/dist/hglib.css">
        <link rel="stylesheet" href="${publicPath}${cssFile}"/>
    </head>
    <body>
        <noscript>You need to enable JavaScript to run this app.</noscript>
        
        <div id="root"></div>
        
        ${peerScripts[nodeEnv]}
        
        <script type="text/javascript" src="${publicPath}${jsFile}"></script>
    </body>
</html>`;
};