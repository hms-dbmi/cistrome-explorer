import d3 from './d3.js';

/**
 * Wrap the SVG element exported by HiGlass, using the wrapper component draw functions.
 * @param {SVGElement} hgSvgNode The SVG node passed to the HiGlass JS API `.on('createSVG')` callback.
 * @param {object} drawInfo Mapping of unique keys to objects with draw functions and options.
 * `{ draw, options: { top, left, width, height } }`
 * @returns {SVGElement} The manipulated SVG node, containing both the HiGlass SVG contents 
 * and the wrapper component SVG contents.
 */
export function wrapSvg(hgSvgNode, drawInfo) {
    const hgSvg = d3.select(hgSvgNode);
    const svg = d3.create("svg");

    // Compute dimensions.
    const hgWidth = parseInt(hgSvg.attr("width"));
    const hgHeight = parseInt(hgSvg.attr("height"));
    const minX = Math.min(0, ...Object.values(drawInfo).map(d => (d.options && d.options.left ? d.options.left : 0)));
    const hgOffsetX = minX * -1;
    const maxX = Math.max(hgOffsetX + hgWidth, ...Object.values(drawInfo).map(d => (d.options && d.options.left && d.options.width ? (hgOffsetX + d.options.left + d.options.width) : 0)));
    const maxY = hgHeight;

    // Set the height and width of the outer SVG.
    svg
        .attr("width", maxX)
        .attr("height", maxY)
        .attr("xmlns", "http://www.w3.org/2000/svg")
        .attr("version", "1.1");

    // Append the HiGlass exported SVG contents, putting them in their own <g> element.
    const hgG = svg.append("g")
        .attr("class", "higlass-g")
        .attr("width", hgSvg.attr("width"))
        .attr("height", hgSvg.attr("height"))
        .attr("transform", `translate(${hgOffsetX},${0})`);
    hgG.html(hgSvgNode.innerHTML);

    // Append the wrapper component visualizations, with a separate <g> element for each draw function.
    for(const [key, d] of Object.entries(drawInfo)) {
        const { top = 0, left = 0, width = 0, height = 0 } = d.options || {};
        const draw = d.draw;

        const tempSvg = d3.create("svg")
            .attr("width", width)
            .attr("height", height);
        
        draw(tempSvg.node());
        
        const g = svg.append("g")
            .attr("class", "wrapper-g")
            .attr("width", width)
            .attr("height", height)
            .attr("transform", `translate(${hgOffsetX+left},${top})`);
        g.html(tempSvg.node().innerHTML);
    }
    
    // Uncomment the following line for development purposes.
    //throw new Error("break");
    return svg.node();
}