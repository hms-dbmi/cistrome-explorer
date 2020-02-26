import d3 from './d3.js';
import { getRetinaRatio } from './canvas.js';

/**
 * Represents a rectangle to be rendered.
 * @param {number} x
 * @param {number} y
 * @param {number} width
 * @param {number} height
 */
export class TwoRectangle {
    constructor(x, y, width, height) {
        this.x = x;
        this.y = y;
        this.width = width;
        this.height = height;

        /** @member {string} */
        this.stroke = null;
        /** @member {string} */
        this.fill = "#000";
        /** Width of the stroke line if stroke is not null. 
         * @member {number} */
        this.linewidth = 1;
        /** @member {number} */
        this.opacity = 1;
        /** In radians. 
         * @member {number} */
        this.rotation = null;
    }
}

/**
 * Represents a circle to be rendered.
 * @param {number} x
 * @param {number} y
 * @param {number} radius
 */
export class TwoCircle {
    constructor(x, y, radius) {
        this.x = x;
        this.y = y;
        this.radius = radius;

        /** @member {string} */
        this.stroke = null;
        /** @member {string} */
        this.fill = "#000";
        /** Width of the stroke line if stroke is not null. 
         * @member {number} */
        this.linewidth = 1;
        /** @member {number} */
        this.opacity = 1;
    }
}

/**
 * Represents a line to be rendered.
 * @param {number} x1
 * @param {number} y1
 * @param {number} x2
 * @param {number} y2
 */
export class TwoLine {
    constructor(x1, y1, x2, y2) {
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;

        /** @member {string} */
        this.stroke = "#000";
        /** Width of the stroke line if stroke is not null. 
         * @member {number} */
        this.linewidth = 1;
        /** @member {number} */
        this.opacity = 1;
    }
}

/**
 * Represents a path to be rendered.
 * @param {number[]} points
 */
export class TwoPath {
    constructor(points) {
        this.points = points;

        /** @member {string} */
        this.stroke = "#000";
        /** Width of the stroke line if stroke is not null. 
         * @member {number} */
        this.linewidth = 1;
        /** @member {number} */
        this.opacity = 1;
    }
}

/**
 * Represents text to be rendered.
 * @param {number} x
 * @param {number} y
 * @param {number} width
 * @param {number} height
 * @param {string} text
 */
export class TwoText {
    constructor(x, y, width, height, text) {
        this.x = x;
        this.y = y;
        this.width = width;
        this.height = height;
        this.text = text;

        /** @member {string} */
        this.fill = "#000";
        /** @member {number} */
        this.fontsize = 14;
        /** @member {string} */
        this.font = "Arial,sans-serif";
        /** Corresponds to canvas `context.textAlign`.
         * Possible values: "start", "middle", "end".
         * @member {string} */
        this.align = "middle";
        /** Corresponds to canvas `context.textBaseline`.
         * Possible values: "alphabetic", "top", "middle", "bottom".
         * @member {string} */
        this.baseline = "alphabetic";
        /** @member {number} */
        this.opacity = 1;
        /** In radians. 
         * @member {number} */
        this.rotation = null;
        /** How text that overflows the bounding box should be dealt with.
         * Possible values: null, "clip", "ellipsis".
         * @member {string} */
        this.overflow = null;
    }


}

/**
 * This class is based on the two.js library: https://two.js.org/
 * It has functionality for rendering to both SVG and canvas.
 * @param {object} config
 * @param {number} config.width The width of the domElement
 * @param {number} config.height The height of the domElement.
 * @param {object} config.domElement The canvas or SVG element.
 */
export default class Two {

    constructor({ width, height, domElement }) {
        this.width = width;
        this.height = height;
        this.domElement = domElement;

        this.elements = [];

        if(!domElement) {
            return;
        }

        switch(domElement.nodeName.toLowerCase()) {
            case 'canvas':
                this.init = this.initCanvas.bind(this);
                this.update = this.updateCanvas.bind(this);
                this.teardown = this.teardownCanvas.bind(this);
                this.measureText = this.measureTextCanvas.bind(this);
                break;
            case 'svg':
                this.init = this.initSvg.bind(this);
                this.update = this.updateSvg.bind(this);
                this.teardown = this.teardownSvg.bind(this);
                this.measureText = this.measureTextSvg.bind(this);
                break;
            default:
                console.warn("Unknown DOM element type.");
        }

        this.init();
    }

    initSvg() {
        this.svg = d3.select(this.domElement);
        this.svg.selectAll("g").remove();
        this.svg
            .attr("width", this.width)
            .attr("height", this.height);
        
        this.g = this.svg
            .append("g")
                .attr("width", this.width)
                .attr("height", this.height);
    }

    initCanvas() {
        const context = this.domElement.getContext('2d');
        const ratio = getRetinaRatio(context);
        const scaledWidth = this.width * ratio;
        const scaledHeight = this.height * ratio;
        this.domElement.setAttribute("width", scaledWidth);
        this.domElement.setAttribute("height", scaledHeight);
        context.scale(ratio, ratio);

        this.context = context;
    }

    /**
     * Create a new rectangle.
     * @param {number} x The x-coordinate for the top left corner of the rect.
     * @param {number} y The y-coordinate for the top left corner of the rect.
     * @param {number} width The width for the rect.
     * @param {number} height The height for the rect.
     * @returns {TwoRectangle} Instance of new rectangle.
     */
    makeRect(x, y, width, height) {
        const rect = new TwoRectangle(x, y, width, height);
        this.elements.push(rect);
        return rect;
    }

    /**
     * Create a new circle.
     * @param {number} x The x-coordinate for the center of the circle.
     * @param {number} y The y-coordinate for the center of the circle.
     * @param {number} radius The radius for the circle.
     * @returns {TwoCircle} Instance of new circle.
     */
    makeCircle(x, y, radius) {
        const circle = new TwoCircle(x, y, radius);
        this.elements.push(circle);
        return circle;
    }

    /**
     * Create a new line.
     * @param {number} x1 The x-coordinate for the line start point.
     * @param {number} y1 The y-coordinate for the line start point.
     * @param {number} x2 The x-coordinate for the line end point.
     * @param {number} y2 The y-coordinate for the line end point.
     * @returns {TwoLine} Instance of new line.
     */
    makeLine(x1, y1, x2, y2) {
        const line = new TwoLine(x1, y1, x2, y2);
        this.elements.push(line);
        return line;
    }

    /**
     * Create a new path.
     * @param {...number} args Coordinates x1, y1, x2, y2, x3, y3, etc.
     * @returns {TwoPath} Instance of new path.
     */
    makePath(...args) {
        const points = [];
        for(let i = 0; i < args.length; i += 2) {
            points.push([ args[i], args[i+1] ]);
        }

        const path = new TwoPath(points);
        this.elements.push(path);
        return path;
    }

    /**
     * Create a new text.
     * @param {number} x The x-coordinate for the anchor point of the text.
     * @param {number} y The y-coordinate for the anchor point of the text.
     * @param {number} width The width for the text.
     * @param {number} height The height for the text.
     * @returns {TwoText} Instance of new text.
     */
    makeText(x, y, width, height, text, withRect = false) {
        const obj = new TwoText(x, y, width, height, text);
        this.elements.push(obj);
        return obj;
    }

    /**
     * Append a TwoText, TwoCircle, TwoRect, TwoPath, etc. to the current list of shape elements.
     * @param {*} obj 
     */
    append(obj) {
        this.elements.push(obj);
    }

    /**
     * Render to the DOM element.
     */
    update() {

    }

    updateSvg() {
        const g = this.g;

        this.elements.forEach((d) => {
            if(d instanceof TwoRectangle) {
                const rect = g.append("rect")
                    .attr("x", d.x)
                    .attr("y", d.y)
                    .attr("width", d.width)
                    .attr("height", d.height)
                    .attr("opacity", d.opacity);
                
                if(d.fill != null) {
                    rect
                        .attr("fill", d.fill);
                } else {
                    rect
                        .attr("fill", "transparent");
                }
                if(d.stroke != null) {
                    rect
                        .attr("stroke-width", d.linewidth)
                        .attr("stroke", d.stroke);
                }
                if(d.rotation != null) {
                    rect
                        .attr("transform", `rotate(${d.rotation * 180/Math.PI},${d.x + d.width/2},${d.y + d.height/2})`);
                }
            } else if(d instanceof TwoCircle) {
                const circle = g.append("circle")
                    .attr("cx", d.x)
                    .attr("cy", d.y)
                    .attr("r", d.radius)
                    .attr("opacity", d.opacity);
                
                if(d.fill != null) {
                    circle
                        .attr("fill", d.fill);
                } else {
                    circle
                        .attr("fill", "transparent");
                }
                if(d.stroke != null) {
                    circle
                        .attr("stroke-width", d.linewidth)
                        .attr("stroke", d.stroke);
                }
            } else if(d instanceof TwoLine) {
                const line = g.append("line")
                    .attr("x1", d.x1)
                    .attr("y1", d.y1)
                    .attr("x2", d.x2)
                    .attr("y2", d.y2)
                    .attr("opacity", d.opacity);
                if(d.stroke != null) {
                    line
                        .attr("stroke-width", d.linewidth)
                        .attr("stroke", d.stroke);
                }
            } else if(d instanceof TwoPath) {
                const path = g.append("path")
                    .attr("opacity", d.opacity);
                
                let pathD = "";
                if(d.points.length > 1) {
                    d.points.forEach((p, i) => {
                        if(i == 0) {
                            pathD += `M ${p[0]} ${p[1]}`;
                        } else {
                            pathD += `L ${p[0]} ${p[1]}`;
                        }
                    });
                }

                if(d.stroke != null) {
                    path
                        .attr("stroke-width", d.linewidth)
                        .attr("stroke", d.stroke);
                }
                path.attr("d", pathD);
            } else if(d instanceof TwoText) {
                const text = g.append("text")
                    .attr("x", d.x)
                    .attr("y", d.y)
                    .attr("text-anchor", d.align)
                    .attr("dominant-baseline", 
                        (d.baseline === "top" ? "text-before-edge" : (d.baseline === "bottom" ? "text-after-edge" : d.baseline))
                    )
                    .attr("opacity", d.opacity)
                    .attr("fill", d.fill)
                    .attr("font-size", d.fontsize)
                    .attr("font-family", d.font)
                    .text(d.text);
                
                let content = d.text;
                if(d.overflow === "clip") {
                    while(content.length > 0 && text.node().getComputedTextLength() > d.width) {
                        content = content.substring(0, content.length - 1);
                        text.text("content");
                    }
                } else if(d.overflow === "ellipsis") {
                    if(text.node().getComputedTextLength() > d.width) {
                        while(content.length > 0 && text.node().getComputedTextLength() > d.width) {
                            content = content.substring(0, content.length - 1);
                            text.text(content + "...");
                        }
                    }
                }
                
                if(d.rotation != null) {
                    text
                        .attr("transform", `rotate(${d.rotation * 180/Math.PI},${d.x},${d.y})`);
                }
            }
        });
    }

    updateCanvas() {
        const context = this.context;
        context.clearRect(0, 0, this.width, this.height);

        this.elements.forEach((d) => {
            context.lineWidth = d.linewidth;
            context.globalAlpha = d.opacity;

            if(d instanceof TwoRectangle) {
                if(d.rotation !== null) {
                    context.save();
                    context.translate(d.x + d.width/2, d.y + d.height/2);
                    context.rotate(d.rotation);
                    context.translate(-(d.x + d.width/2), -(d.y + d.height/2));
                }
                if(d.fill !== null) {
                    context.fillStyle = d.fill;
                    context.fillRect(d.x, d.y, d.width, d.height);
                }
                if(d.stroke !== null) {
                    context.strokeStyle = d.stroke;
                    context.strokeRect(d.x, d.y, d.width, d.height);
                }
                if(d.rotation !== null) {
                    context.restore();
                }
            } else if(d instanceof TwoCircle) {
                if(d.fill !== null) {
                    context.fillStyle = d.fill;
                }
                if(d.stroke !== null) {
                    context.strokeStyle = d.stroke;
                }
                context.beginPath();
                context.arc(d.x, d.y, d.radius, 0, 2*Math.PI);
                if(d.fill !== null) {
                    context.fill();
                }
                if(d.stroke !== null) {
                    context.stroke();
                }
            } else if(d instanceof TwoLine) {
                if(d.stroke !== null) {
                    context.strokeStyle = d.stroke;
                }
                context.beginPath();
                context.moveTo(d.x1, d.y1);
                context.lineTo(d.x2, d.y2);
                if(d.stroke !== null) {
                    context.stroke();
                }
            } else if(d instanceof TwoPath) {
                if(d.stroke !== null) {
                    context.strokeStyle = d.stroke;
                }
                if(d.points.length > 1) {
                    context.beginPath();
                    d.points.forEach((p, i) => {
                        if(i == 0) {
                            context.moveTo(p[0], p[1]);
                        } else {
                            context.lineTo(p[0], p[1]);
                        }
                    });
                }
                if(d.stroke !== null) {
                    context.stroke();
                }
            } else if(d instanceof TwoText) {
                if(d.rotation !== null) {
                    context.save();
                    context.translate(d.x, d.y);
                    context.rotate(d.rotation);
                    context.translate(-d.x, -d.y);
                }
                context.font = `${d.fontsize}px ${d.font}`;
                context.fillStyle = d.fill;
                context.textAlign = (d.align === "middle" ? "center" : d.align);
                context.textBaseline = d.baseline;

                let content = d.text;
                if(d.overflow === "clip") {
                    while(content.length > 0 && context.measureText(content).width > d.width) {
                        content = content.substring(0, content.length - 1);
                    }
                } else if(d.overflow === "ellipsis") {
                    if(context.measureText(content).width > d.width) {
                        while(content.length > 0 && context.measureText(content + "...").width > d.width) {
                            content = content.substring(0, content.length - 1);
                        }
                        content = content + "...";
                    }
                }

                context.fillText(content, d.x, d.y);
                if(d.rotation !== null) {
                    context.restore();
                }
            }
        });
    }

    /**
     * Compute width and height for a particular text element.
     * @param {TwoText} d The object representing the text to be measured.
     * @returns {object} Object containing the values `width` and `height`.
     */
    measureText(d) {
        return { width: 0, height: 0 };
    }

    measureTextSvg(d) {
        const svg = d3.create("svg");
        const g = svg.append("g");

        const text = g.append("text")
            .attr("x", d.x)
            .attr("y", d.y)
            .attr("text-anchor", d.align)
            .attr("dominant-baseline", 
                (d.baseline === "top" ? "text-before-edge" : (d.baseline === "bottom" ? "text-after-edge" : d.baseline))
            )
            .attr("opacity", d.opacity)
            .attr("fill", d.fill)
            .attr("font-size", d.fontsize)
            .attr("font-family", d.font)
            .text(d.text);
        
        let content = d.text;
        if(d.overflow === "clip") {
            while(content.length > 0 && text.node().getComputedTextLength() > d.width) {
                content = content.substring(0, content.length - 1);
                text.text("content");
            }
        } else if(d.overflow === "ellipsis") {
            if(text.node().getComputedTextLength() > d.width) {
                while(content.length > 0 && text.node().getComputedTextLength() > d.width) {
                    content = content.substring(0, content.length - 1);
                    text.text(content + "...");
                }
            }
        }
        
        if(d.rotation != null) {
            text
                .attr("transform", `rotate(${d.rotation * 180/Math.PI},${d.x},${d.y})`);
        }

        // Measure the dimensions.
        const dims = { width: 0, height: 0 };
        try {
            const { width, height } = text.node().getBBox();
            dims.width = width;
            dims.height = height;
        } catch(e) {
            console.log(e);
        }

        return dims;
    }

    measureTextCanvas(d) {
        const context = this.context;
        if(d.rotation !== null) {
            context.save();
            context.translate(d.x, d.y);
            context.rotate(d.rotation);
            context.translate(-d.x, -d.y);
        }
        context.font = `${d.fontsize}px ${d.font}`;
        context.fillStyle = d.fill;
        context.textAlign = (d.align === "middle" ? "center" : d.align);
        context.textBaseline = d.baseline;

        let content = d.text;
        if(d.overflow === "clip") {
            while(content.length > 0 && context.measureText(content).width > d.width) {
                content = content.substring(0, content.length - 1);
            }
        } else if(d.overflow === "ellipsis") {
            if(context.measureText(content).width > d.width) {
                while(content.length > 0 && context.measureText(content + "...").width > d.width) {
                    content = content.substring(0, content.length - 1);
                }
                content = content + "...";
            }
        }

        // Measure the dimensions.
        const dims = { width: 0, height: 0 };
        try {
            const metrics = context.measureText(content);
            dims.width = metrics.width;
            dims.height = (metrics.actualBoundingBoxAscent + metrics.actualBoundingBoxDescent);
        } catch(e) {
            console.log(e);
        }
        
        if(d.rotation !== null) {
            context.restore();
        }

        return dims;
    }

    

    /**
     * Clean up the DOM element (remove event listeners, etc.).
     */
    teardown() {

    }

    teardownSvg() {

    }

    teardownCanvas() {
        const canvasSelection = d3.select(this.domElement);
        canvasSelection.on("mousemove", null);
        canvasSelection.on("mouseout", null);
        canvasSelection.on("click", null);
    }
}