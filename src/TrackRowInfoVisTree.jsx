import React, { useRef, useCallback, useEffect, useState, useMemo, useLayoutEffect } from 'react';
import range from 'lodash/range';
import PubSub from 'pubsub-js';

import d3 from './utils/d3.js';
import Two from './utils/two.js';
import { CONTEXT_MENU_TYPE, EVENT } from './utils/constants.js';
import { TooltipContent, destroyTooltip } from './Tooltip.jsx';
import { drawVisTitle } from './utils/vis.js';

import TrackRowInfoControl from './TrackRowInfoControl.jsx';
import { rgbToHex, generateNextUniqueColor } from './utils/color.js';
import { getRetinaRatio } from './utils/canvas.js';
import { modifyItemInArray } from './utils/array.js';
import { getAggregatedValue } from './utils/aggregate.js';
import { drawRowHighlightRect } from './utils/linking.js';
import { AnnDataSource, ObsFeatureMatrixAnndataLoader } from '@vitessce/zarr';
import { treegl } from 'tree.gl';
import { ARROW_H, COMPRESS, FILTER } from './utils/icons.js';

export const margin = 5;

/**
 * Component for visualization of row info quantitative or nominal attribute values.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {number} titleHeight The height of the track title.
 * @prop {object} fieldInfo The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 * @prop {boolean} isShowControlButtons Determine if control buttons should be shown.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object[]} transformedRowInfo The `rowInfo` array after aggregating, filtering, and sorting rows.
 * @prop {array} selectedRows The array of selected indices.
 * @prop {array} highlitRows The array of highlit indices.
 * @prop {string} titleSuffix The suffix of a title, information about sorting and filtering status.
 * @prop {object} sortInfo The options for sorting rows of the field used in this track.
 * @prop {object} filterInfo The options for filtering rows of the field used in this track.
 * @prop {function} onAddTrack The function to call upon a track insertion.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onHighlightRows The function to call upon a highlight interaction.
 * @prop {function} onFilterRows The function to call upon a filter interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisTree(props) {
	const {
		left,
		top,
		width,
		height,
		titleHeight,
		url,
		field,
		type,
		alt,
		title,
		color: constantColor,
		colors,
		only,
		aggFunction,
		resolveYScale,
		isLeft,
		isShowControlButtons,
		rowInfo,
		transformedRowInfo,
		selectedRows,
		highlitRows,
		titleSuffix,
		sortInfo,
		filterInfo,
		onAddTrack,
		onSortRows,
		onHighlightRows,
		onFilterRows,
		helpActivated,
		drawRegister
	} = props;

	const divRef = useRef();
	const canvasRef = useRef();
	const titleCanvasRef = useRef();
	const mousePositionRef = useRef();
	const treeglRef = useRef();

	// Data, layouts and styles
	const textAreaWidth = width >= 60 ? 20 : 0;
	const barAreaWidth = width - textAreaWidth;
	const minTrackWidth = 40;
	const fontSize = 10;

	let xScale = d3.scaleLinear();
	const yScale = d3.scaleBand().domain(range(transformedRowInfo.length)).range([titleHeight, height]);
	const rowHeight = yScale.bandwidth();

	// Array to store information for mouse events, such as unique color.
	let colorToInfo = [];
	let cnt = 1,
		fields = [field];
    useLayoutEffect(() => {
		const two = new Two({
			width,
			height,
			domElement: titleCanvasRef.current
		});
		drawVisTitle(title, { two, isLeft, width, height, titleSuffix });
		two.update();

        canvasRef.current.width = width - 30;
        canvasRef.current.height = height - 30;
        const attributes = {
            a: { customDistance: [0, 15] },
        
            b: { customDistance: [0, 5] },
            c: { customDistance: [5, 5 + 18] },
            
            d: { customDistance: [0, 112] },
            e: { customDistance: [112, 112 + 350] }
        }
        treeglRef.current = treegl(canvasRef.current, {
            data: {
                tree: '(b:0.1,(d:0.2,e:0.5)c:0.2)a:0.1;',
                attributes
            },
            encoding: {
                layout: 'linear',
                nodes: 'all',
                nodeSize: 10,
                edgeLength: 'customDistance'
            },
            style: {
                padding: 1,
                edge: 'orthogonal',
                edgeToBaseline: false,
				background: "#F6F6F6"
            },
        }).clicked(opt => {
			const { x, y, spec, node } = opt;
			const { fold, remove, attributes: currAttr } = spec.data;
			if(node?.name) {
				// Find the node name
				const mapping = {
					a: 'Astrocyte, Excitatory, Inhibitory',
					b: 'Astrocyte',
					c: 'Excitatory, Inhibitory',
					d: 'Inhibitory',
					e: 'Excitatory'
				}
				const clickedNode = node.name;

				// collapse the node
				const collapseNode = () => {
					if(clickedNode === 'a') {
						const newFold = [...fold, 'a'];
						onFilterRows(
							field,
							'nominal',
							[...newFold.map(node => mapping[node])],
							false
						);
						treeglRef.current.update({
							data: {
								fold: newFold,
								// remove: [...remove],
								// attributes: structuredClone(currAttr)
							},
							transition: {
								duration: 0
							}
						});
						return;
					}
					const editedAttr = structuredClone(currAttr);
					const editMap = {
						b: 'c',
						c: 'b',
						d: 'e',
						e: 'd'
					}
					const editNode = editMap[clickedNode];
					editedAttr[editNode].customDistance = [0, editedAttr[editNode].customDistance[1] - editedAttr[editNode].customDistance[0]];
					if(clickedNode === 'c') {
						editedAttr.d.customDistance = [editedAttr.d.customDistance[0] + 5, editedAttr.d.customDistance[1] + 5];
						editedAttr.e.customDistance = [editedAttr.e.customDistance[0] + 5, editedAttr.e.customDistance[1] + 5];
					}
					const newFold = [...fold, clickedNode];
					onFilterRows(
						field,
						'nominal',
						[...newFold.map(d => mapping[d])],
						false
					);
					treeglRef.current.update({
						data: {
							fold: newFold,
							remove: [...remove],
							attributes: editedAttr
						},
						transition: {
							duration: 0
						}
					});
				};
				const filterNode = () => {
					const removeMapping = {
						a: ['Astrocyte, Excitatory, Inhibitory', 'Astrocyte', 'Excitatory, Inhibitory', 'Inhibitory', 'Excitatory'],
						b: ['Astrocyte'],
						c: ['Excitatory, Inhibitory','Inhibitory', 'Excitatory'],
						d: ['Inhibitory'],
						e: ['Excitatory']
					}

					const editedAttr = structuredClone(currAttr);
					const editMap = {
						b: 'c',
						d: 'e',
					}
					const editNode = editMap[clickedNode];
					if(editNode) {
						editedAttr[editNode].customDistance = [0, editedAttr[editNode].customDistance[1] - editedAttr[editNode].customDistance[0]];
					}
					const newRemove = [...remove, clickedNode];
					console.log([...fold.map(d => mapping[d]), ...newRemove.flatMap(d => removeMapping[d])]);
					onFilterRows(
						field,
						'nominal',
						[...fold.map(d => mapping[d]), ...newRemove.flatMap(d => removeMapping[d])],
						false
					);
					treeglRef.current.update({
						data: {
							fold: [...fold],
							remove: newRemove,
							attributes: editedAttr
						},
						transition: {
							duration: 0
						}
					});
				}
				setTimeout(() => {
					// click event removes context menus, so this is a workaround to show context menu
					PubSub.publish(EVENT.CONTEXT_MENU, {
						x: mousePositionRef.current.x ?? 0,
						y: mousePositionRef.current.y ?? 0,
						menuType: CONTEXT_MENU_TYPE.NOMINAL_BAR,
						title: `Selected category: ${mapping[clickedNode]}`,
						items: [
							{ title: 'Collapse a Node', icon: COMPRESS, action: () => collapseNode() },
							{ title: 'Filter a Node and Child Nodes', icon: FILTER, action: () => filterNode() }
						]
					});
				}, 100);
			} else {
				// remove filtering
				onFilterRows(
					field,
					'nominal',
					rowInfo.map(d => getAggregatedValue(d, field, 'nominal', aggFunction).toString()),
					true
				);
				return {
					data: {
						fold: [],
						remove: [],
						attributes
					}
				}
			}
			

            
            // const clickedNode = !node?.name ? undefined : node.name;
			
        });
    }, []);

	// const draw = useCallback((domElement, isHidden) => {
	// 	const two = new Two({
	// 		width,
	// 		height,
	// 		domElement
	// 	});

	// 	const isTextLabel = width > minTrackWidth;

	// 	// Scales
	// 	const valueExtent = [0, d3.extent(transformedRowInfo.map(d => getExpressionValue(d, field)))[1]]; // Zero baseline
	// 	xScale = d3.scaleLinear().domain(valueExtent).range([0, barAreaWidth]);
	// 	const colorScale = d3.scaleLinear().domain(valueExtent).range([0, 1]);

	// 	// Render visual components for each row (i.e., bars and texts).
	// 	const textAlign = isLeft ? 'end' : 'start';
	// 	transformedRowInfo.forEach((d, i) => {
	// 		const isSkip = only && only.value !== getExpressionValue(d, only.field);
	// 		if (isSkip) return;

	// 		const value = getExpressionValue(d, field);
	// 		const barTop = yScale(i);
	// 		const barWidth = xScale(value);
	// 		const barLeft = isLeft ? width - barWidth : 0;
	// 		const textLeft = isLeft ? width - barWidth - margin : barWidth + margin;
	// 		const infoForMouseEvent = colorToInfo.find(d => d.field === field && d.rowIndex === i);
	// 		const color = isHidden ? infoForMouseEvent.uniqueColor : d3.interpolateViridis(colorScale(value));

	// 		const rect = two.makeRect(barLeft, barTop, barWidth, rowHeight);
	// 		rect.fill = constantColor || color;

	// 		// Render text labels when the space is enough.
	// 		if (rowHeight >= fontSize && isTextLabel) {
	// 			const text = two.makeText(
	// 				textLeft,
	// 				barTop + rowHeight / 2,
	// 				textAreaWidth,
	// 				rowHeight,
	// 				numberFormatShort(value)
	// 			);
	// 			text.fill = d3.hsl(color).darker(3);
	// 			text.fontsize = fontSize;
	// 			text.align = textAlign;
	// 			text.baseline = 'middle';
	// 			text.overflow = 'clip';
	// 		}

	// 		// Add other vis properties to colorToInfo for mouse events.
	// 		if (!isHidden) {
	// 			colorToInfo = modifyItemInArray(colorToInfo, colorToInfo.indexOf(infoForMouseEvent), {
	// 				...infoForMouseEvent,
	// 				color
	// 			});
	// 		}
	// 	});

	// 	drawVisTitle(title, { two, isLeft, width, height, titleSuffix });

	// 	two.update();
	// 	return two.teardown;
	// });

	// const drawAxis = useCallback(domElement => {
	// 	if (width <= 60) return () => {};
	// 	d3.select(domElement).selectAll('*').remove();

	// 	const axisScale = isLeft ? xScale.domain(xScale.domain().reverse()) : xScale;
	// 	const axis = d3
	// 		.axisBottom(axisScale)
	// 		.ticks(Math.ceil(barAreaWidth / 40))
	// 		.tickFormat(d3.format('.2s'));

	// 	d3.select(domElement)
	// 		.attr('width', width)
	// 		.attr('height', axisHeight + titleHeight)
	// 		.append('g')
	// 		.attr('transform', `translate(${isLeft ? textAreaWidth - 1 : 1}, ${titleHeight})`)
	// 		.call(axis);

	// 	d3.select(domElement)
	// 		.selectAll('text')
	// 		.attr('transform', `translate(${isLeft ? -3 : 3}, 0)`);

	// 	return () => {
	// 		/* Teardown */
	// 	};
	// });

	// drawRegister('TrackRowInfoVisQuantitativeBar', draw);
	// drawRegister('TrackRowInfoVisQuantitativeBarAxis', drawAxis);

	// useEffect(() => {
	// 	const canvas = canvasRef.current;
	// 	const hiddenCanvas = hiddenCanvasRef.current;
	// 	const svg = axisRef.current;
	// 	const div = divRef.current;
	// 	const teardown = draw(canvas);
	// 	const teardownHidden = draw(hiddenCanvas, true);
	// 	const teardownSvg = drawAxis(svg);

	// 	d3.select(canvas).on('mousemove', () => {
	// 		const [mouseX, mouseY] = d3.mouse(canvas);

	// 		const hiddenContext = hiddenCanvasRef.current.getContext('2d');
	// 		const ratio = getRetinaRatio(hiddenContext);
	// 		const uniqueColor = rgbToHex(hiddenContext.getImageData(mouseX * ratio, mouseY * ratio, 1, 1).data);
	// 		const hoveredInfo = colorToInfo.find(d => d.uniqueColor === uniqueColor);

	// 		const mouseViewportX = d3.event.clientX;
	// 		const mouseViewportY = d3.event.clientY;

	// 		if (hoveredInfo) {
	// 			PubSub.publish(EVENT.TOOLTIP, {
	// 				x: mouseViewportX,
	// 				y: mouseViewportY,
	// 				content: (
	// 					<TooltipContent
	// 						title={title}
	// 						value={numberFormatLong(hoveredInfo.value)}
	// 						color={hoveredInfo.color}
	// 					/>
	// 				)
	// 			});
	// 			onHighlightRows(field, 'quantitative', [hoveredInfo.value, hoveredInfo.value]);
	// 		} else {
	// 			destroyTooltip();
	// 			onHighlightRows('');
	// 		}
	// 	});

	// 	// Handle mouse enter and leave.
	// 	d3.select(canvas).on('mouseout', destroyTooltip);
	// 	d3.select(div).on('mouseleave', () => {
	// 		onHighlightRows('');
	// 	});

	// 	// Clean up.
	// 	return () => {
	// 		teardown();
	// 		teardownHidden();
	// 		teardownSvg();
	// 		d3.select(div).on('mouseenter', null);
	// 		d3.select(div).on('mouseleave', null);
	// 	};
	// }, [top, left, width, height, transformedRowInfo, isShowControlButtons]);

	return (
		<div
			ref={divRef}
			style={{
				top: `${top}px`,
				position: 'relative',
				width: `${width}px`,
				height: `${height}px`,
			}}
			onMouseMove={e => mousePositionRef.current = { x: e.clientX, y: e.clientY }}
		>
			<canvas
				ref={titleCanvasRef}
				style={{
					top: 0,
					left: 0,
					width: `${width}px`,
					height: `${height}px`,
					position: 'absolute',
				}}
			/>
			<canvas
				ref={canvasRef}
				style={{
					top: 30,
					left: 0,
					width: `${width - 30}px`,
					height: `${height - 30}px`,
					position: 'absolute',
				}}
			/>
			{/* <TrackRowInfoControl
                isLeft={isLeft}
                isVisible={isShowControlButtons}
                field={field}
                type={type}
                title={title}
                aggFunction={aggFunction}
                top={titleHeight}
                searchLeft={left}
                sortAsceButtonHighlit={sortInfo && sortInfo.order === "ascending"}
                sortDescButtonHighlit={sortInfo && sortInfo.order === "descending"}
                filterButtonHighlit={filterInfo !== undefined}
                onSortRows={onSortRows}
                onHighlightRows={onHighlightRows}
                onFilterRows={onFilterRows}
                filterInfo={filterInfo}
                rowInfo={rowInfo}
                transformedRowInfo={transformedRowInfo}
                helpActivated={helpActivated}
            /> */}
		</div>
	);
}
