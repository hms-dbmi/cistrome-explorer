import d3 from './d3.js';
import { getRetinaRatio } from './canvas.js';


class TwoRectangle {
    constructor(x, y, width, height) {
        this.x = x;
        this.y = y;
        this.width = width;
        this.height = height;

        this.stroke = null;
        this.fill = "#000";
        this.linewidth = 1;
        this.opacity = 1;
        this.rotation = null;
    }
}

class TwoCircle {
    constructor(x, y, size) {
        this.x = x;
        this.y = y;
        this.size = size;

        this.stroke = null;
        this.fill = "#000";
        this.linewidth = 1;
        this.opacity = 1;
    }
}

class TwoLine {
    constructor(x1, y1, x2, y2) {
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;

        this.stroke = "#000";
        this.linewidth = 1;
        this.opacity = 1;
    }
}

class TwoPath {
    constructor(points) {
        this.points = points;

        this.stroke = "#000";
        this.linewidth = 1;
        this.opacity = 1;
    }
}

class TwoText {
    constructor(x, y, width, height, text) {
        this.x = x;
        this.y = y;
        this.width = width;
        this.height = height;
        this.text = text

        this.fill = "#000";
        this.fontsize = 14;
        this.font = "Arial,sans-serif";
        this.align = "middle"; // options: "start", "middle", "end"
        this.baseline = "alphabetic"; // options "alphabetic", "top", "middle", "bottom"
        this.linewidth = 1;
        this.opacity = 1;
        this.rotation = null;
    }
}

/**
 * This class is based on the two.js library: https://two.js.org/
 * It has functionality for rendering to both SVG and canvas.
 * @param {object} config
 * @param {number} config.width
 * @param {number} config.height
 * @param {object} config.domElement
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
                break;
            case 'svg':
                this.init = this.initSvg.bind(this);
                this.update = this.updateSvg.bind(this);
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
     * @returns {object} Instance of new `TwoRectangle`.
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
     * @param {number} size The diameter for the circle.
     * @returns {object} Instance of new `TwoCircle`.
     */
    makeCircle(x, y, size) {
        const circle = new TwoCircle(x, y, size);
        this.elements.push(circle);
        return circle;
    }

    /**
     * Create a new line.
     * @param {number} x1 The x-coordinate for the line start point.
     * @param {number} y1 The y-coordinate for the line start point.
     * @param {number} x2 The x-coordinate for the line end point.
     * @param {number} y2 The y-coordinate for the line end point.
     * @returns {object} Instance of new `TwoLine`.
     */
    makeLine(x1, y1, x2, y2) {
        const line = new TwoLine(x1, y1, x2, y2);
        this.elements.push(line);
        return line;
    }

    /**
     * Create a new path.
     * @param {...number} coord Coordinates x1, y1, x2, y2, x3, y3, etc.
     * @returns {object} Instance of new `TwoPath`.
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
     * @returns {object} Instance of new `TwoText`.
     */
    makeText(x, y, width, height, text) {
        const obj = new TwoText(x, y, width, height, text);
        this.elements.push(obj);
        return obj;
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
                    .attr("r", d.size)
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
                context.arc(d.x, d.y, d.size, 0, 2*Math.PI);
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
                context.fillText(d.text, d.x, d.y);
                if(d.rotation !== null) {
                    context.restore();
                }
            }
        });
    }
}