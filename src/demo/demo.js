import { hgDemoViewConfigAtac } from '../viewconfigs/horizontal-multivec-atac';
import { hgDemoViewConfigAtacRevision } from '../viewconfigs/horizontal-multivec-atac-revision';
import { hgDemoViewConfig3k27 } from '../viewconfigs/horizontal-multivec-3k27';
import { hgDemoViewConfig3k27Revision } from '../viewconfigs/horizontal-multivec-3k27-revision';
import { hgDemoViewConfig3K4 } from '../viewconfigs/horizontal-multivec-3k4';
import {
	hgDemoViewConfigMiraMouse,
	hgDemoViewConfigMiraMouse4000,
	hgDemoViewConfigMiraMouse500Smooth
} from '../viewconfigs/horizontal-multivec-mira-mouse';
import { hgDemoViewConfig3K4Revision } from '../viewconfigs/horizontal-multivec-3k4-revision';

const STATES = [
	'Astrocyte, Excitatory, Inhibitory',
	'Astrocyte',
	'Excitatory, Inhibitory',
	'Inhibitory',
	'Excitatory',
	''
];
const STATE_COLORS = ['#D6641E', '#E6A01B', '#3275B4', '#409F7A', '#CC7DAA', 'transparent'];
const TOPIC_COLORS = [
	'#FF9896',
	'#2077B4',
	'#7F7F7F',
	'#AEC7E8',
	'#BCBD22',
	'#F7B6D2',
	'#DBDB8E',
	'#E377C2',
	'#C5B0D5',
	'#D62828',
	'#2CA02C',
	'grey',
	'#9FDBE5'
];
const streamBase = {
	field: [
		'atac_topic_11',
		'atac_topic_5',
		'atac_topic_6',
		'atac_topic_2',
		'atac_topic_8',
		'atac_topic_7',
		'atac_topic_0',
		'atac_topic_9',
		'atac_topic_10',
		'atac_topic_1',
		'atac_topic_12',
		'atac_topic_4',
		'atac_topic_3'
	],
	type: 'quantitative',
	position: 'right',
	width: 30,
	colors: TOPIC_COLORS
};

export const demos = {
	'ATAC (v0)': {
		viewConfig: hgDemoViewConfigAtac,
		options: [
			{
				viewId: 'cistrome-view-atac',
				trackId: 'cistrome-track-atac',
				rowInfoAttributes: [
					{ field: 'Cell Type', type: 'nominal', position: 'right', width: 200 },
					// {field: "Tissue Type", type: "nominal", position: "right", width: 120},
					{ field: 'Clustering', type: 'tree', position: 'right', width: 200 },
					// { field: "FRACTION_PEAKS_IN_EXON", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_PEAKS_IN_INTRON", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_PEAKS_IN_PROMOTER", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_PEAKS_IN_INTERGENIC", type: "quantitative", position: "right", width: 80 },
					{ field: 'PEAKS_TOTAL', type: 'quantitative', title: 'Total peaks', position: 'right', width: 100 },
					// { field: "PEAKS_10FOLD", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_READS_IN_PEAKS", type: "quantitative", position: "right", width: 80 },
					// { field: "FASTQC_MEDIAN", type: "quantitative", position: "right", width: 80 },
					{
						field: 'SEQUENCE_LENGTH',
						type: 'quantitative',
						title: 'Sequence length',
						position: 'right',
						width: 100
					},
					// { field: "FRACTION_PEAKS_IN_DHS", type: "quantitative", position: "right", width: 80 },
					// { field: "PCR_BOTTLENECK_COEF", type: "quantitative", position: "right", width: 80 },
					// { field: "READS_TOTAL", type: "quantitative", position: "right", width: 80 },
					{
						field: 'READS_MAPPED',
						type: 'quantitative',
						title: 'Reads mapped',
						position: 'right',
						width: 100
					}
					// { field: "PHYLOCONS_0_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
					// { field: "PHYLOCONS_500_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
					// { field: "PHYLOCONS_1000_BP_OFFSET", type: "quantitative", position: "right", width: 80 }
					// {field: "Species", type: "nominal", position: "left", width: 120},
					// {field: "cid", type: "nominal-dynamic",  position: "left", title: "Compare Positive and Negative", domain: ["positive", "negative"], range: ["blue", "red"], width: 30},
					// {field: "ID", type: "url", shortName: "📊", position: "right", title: "➕ Add Detail Track", width: 35, addTrackOnClick: true},
					// {field: "Metadata URL", alt: "cid", type: "url", position: "right", width: 70},
				],
				rowAggregate: [
					// {field: "Cell Type", type: "nominal", notOneOf: []},
					// {field: "Tissue Type", type: "nominal", notOneOf: []}
				],
				rowSort: [{ field: 'Clustering', type: 'tree', order: 'ascending' }],
				rowFilter: [
					// ...
				]
			}
		]
	},
	'ATAC (v1)': {
		viewConfig: hgDemoViewConfigAtacRevision,
		options: [
			{
				viewId: 'cistrome-view-atac-revision',
				trackId: 'cistrome-track-atac-revision',
				rowInfoAttributes: [
					{ field: 'Cell Type', type: 'nominal', position: 'right', width: 200 },
					// {field: "Tissue Type", type: "nominal", position: "right", width: 120},
					{ field: 'Clustering', type: 'tree', position: 'right', width: 200 },
					// { field: "FRACTION_PEAKS_IN_EXON", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_PEAKS_IN_INTRON", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_PEAKS_IN_PROMOTER", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_PEAKS_IN_INTERGENIC", type: "quantitative", position: "right", width: 80 },
					{ field: 'PEAKS_TOTAL', type: 'quantitative', title: 'Total peaks', position: 'right', width: 100 },
					// { field: "PEAKS_10FOLD", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_READS_IN_PEAKS", type: "quantitative", position: "right", width: 80 },
					// { field: "FASTQC_MEDIAN", type: "quantitative", position: "right", width: 80 },
					{
						field: 'SEQUENCE_LENGTH',
						type: 'quantitative',
						title: 'Sequence length',
						position: 'right',
						width: 100
					},
					// { field: "FRACTION_PEAKS_IN_DHS", type: "quantitative", position: "right", width: 80 },
					// { field: "PCR_BOTTLENECK_COEF", type: "quantitative", position: "right", width: 80 },
					// { field: "READS_TOTAL", type: "quantitative", position: "right", width: 80 },
					{
						field: 'READS_MAPPED',
						type: 'quantitative',
						title: 'Reads mapped',
						position: 'right',
						width: 100
					}
					// { field: "PHYLOCONS_0_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
					// { field: "PHYLOCONS_500_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
					// { field: "PHYLOCONS_1000_BP_OFFSET", type: "quantitative", position: "right", width: 80 }
					// {field: "Species", type: "nominal", position: "left", width: 120},
					// {field: "cid", type: "nominal-dynamic",  position: "left", title: "Compare Positive and Negative", domain: ["positive", "negative"], range: ["blue", "red"], width: 30},
					// {field: "ID", type: "url", shortName: "📊", position: "right", title: "➕ Add Detail Track", width: 35, addTrackOnClick: true},
					// {field: "Metadata URL", alt: "cid", type: "url", position: "right", width: 70},
				],
				rowAggregate: [
					// {field: "Cell Type", type: "nominal", notOneOf: []},
					// {field: "Tissue Type", type: "nominal", notOneOf: []}
				],
				rowSort: [{ field: 'Clustering', type: 'tree', order: 'ascending' }],
				rowFilter: [
					// ...
				]
			}
		]
	},
	// "H3K27ac (v0)": {
	//     viewConfig: hgDemoViewConfig3k27,
	//     options: [{
	//         viewId: "cistrome-view-3k27",
	//         trackId: "cistrome-track-3k27",
	//         rowInfoAttributes: [
	//             {field: "Cell Type", type: "nominal", position: "right", width: 200},
	//             {field: "Clustering", type: "tree", position: "right", width: 200},
	//             // { field: "FRACTION_PEAKS_IN_EXON", type: "quantitative", position: "right", width: 80 },
	//             // { field: "FRACTION_PEAKS_IN_INTRON", type: "quantitative", position: "right", width: 80 },
	//             // { field: "FRACTION_PEAKS_IN_PROMOTER", type: "quantitative", position: "right", width: 80 },
	//             // { field: "FRACTION_PEAKS_IN_INTERGENIC", type: "quantitative", position: "right", width: 80 },
	//             { field: "PEAKS_TOTAL", type: "quantitative", title: 'Total peaks', position: "right", width: 100 },
	//             // { field: "PEAKS_10FOLD", type: "quantitative", position: "right", width: 80 },
	//             // { field: "FRACTION_READS_IN_PEAKS", type: "quantitative", position: "right", width: 80 },
	//             // { field: "FASTQC_MEDIAN", type: "quantitative", position: "right", width: 80 },
	//             { field: "SEQUENCE_LENGTH", type: "quantitative", title: 'Sequence length', position: "right", width: 100 },
	//             // { field: "FRACTION_PEAKS_IN_DHS", type: "quantitative", position: "right", width: 80 },
	//             // { field: "PCR_BOTTLENECK_COEF", type: "quantitative", position: "right", width: 80 },
	//             // { field: "READS_TOTAL", type: "quantitative", position: "right", width: 80 },
	//             { field: "READS_MAPPED", type: "quantitative", title: 'Reads mapped', position: "right", width: 100 },
	//             // { field: "PHYLOCONS_0_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
	//             // { field: "PHYLOCONS_500_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
	//             // { field: "PHYLOCONS_1000_BP_OFFSET", type: "quantitative", position: "right", width: 80 }
	//         ],
	//         rowAggregate: [
	//         ],
	//         rowSort: [
	//             {field: "Clustering", type: "tree", order: "ascending"}
	//         ],
	//         rowFilter: [
	//         ]
	//     }]
	// },
	'H3K27ac (v1)': {
		viewConfig: hgDemoViewConfig3k27Revision,
		options: [
			{
				viewId: 'cistrome-view-3k27-revision',
				trackId: 'cistrome-track-3k27-revision',
				rowInfoAttributes: [
					{ field: 'Cell Type', type: 'nominal', position: 'right', width: 200 },
					{ field: 'Clustering', type: 'tree', position: 'right', width: 200 },
					// { field: "FRACTION_PEAKS_IN_EXON", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_PEAKS_IN_INTRON", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_PEAKS_IN_PROMOTER", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_PEAKS_IN_INTERGENIC", type: "quantitative", position: "right", width: 80 },
					{ field: 'PEAKS_TOTAL', type: 'quantitative', title: 'Total peaks', position: 'right', width: 100 },
					// { field: "PEAKS_10FOLD", type: "quantitative", position: "right", width: 80 },
					// { field: "FRACTION_READS_IN_PEAKS", type: "quantitative", position: "right", width: 80 },
					// { field: "FASTQC_MEDIAN", type: "quantitative", position: "right", width: 80 },
					{
						field: 'SEQUENCE_LENGTH',
						type: 'quantitative',
						title: 'Sequence length',
						position: 'right',
						width: 100
					},
					// { field: "FRACTION_PEAKS_IN_DHS", type: "quantitative", position: "right", width: 80 },
					// { field: "PCR_BOTTLENECK_COEF", type: "quantitative", position: "right", width: 80 },
					// { field: "READS_TOTAL", type: "quantitative", position: "right", width: 80 },
					{
						field: 'READS_MAPPED',
						type: 'quantitative',
						title: 'Reads mapped',
						position: 'right',
						width: 100
					}
					// { field: "PHYLOCONS_0_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
					// { field: "PHYLOCONS_500_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
					// { field: "PHYLOCONS_1000_BP_OFFSET", type: "quantitative", position: "right", width: 80 }
				],
				rowAggregate: [],
				rowSort: [{ field: 'Clustering', type: 'tree', order: 'ascending' }],
				rowFilter: []
			}
		]
	},
	// "H3K4me3 (v0)": {
	//     viewConfig: hgDemoViewConfig3K4,
	//     options: [{
	//         viewId: "cistrome-view-3k4",
	//         trackId: "cistrome-track-3k4",
	//         rowInfoAttributes: [
	//             {field: "Cell Type", type: "nominal", position: "right", width: 200},
	//             {field: "Clustering", type: "tree", position: "right", width: 200},
	//             { field: "PEAKS_TOTAL", type: "quantitative", title: 'Total peaks', position: "right", width: 100 },
	//             { field: "SEQUENCE_LENGTH", type: "quantitative", title: 'Sequence length', position: "right", width: 100 },
	//             { field: "READS_MAPPED", type: "quantitative", title: 'Reads mapped', position: "right", width: 100 },
	//         ],
	//         rowAggregate: [],
	//         rowSort: [
	//             {field: "Clustering", type: "tree", order: "ascending"}
	//             // {field: "Cell Type", type: "nominal", order: "ascending"}
	//         ],
	//         rowFilter: []
	//     }]
	// },
	'H3K4me3 (v1)': {
		viewConfig: hgDemoViewConfig3K4Revision,
		options: [
			{
				viewId: 'cistrome-view-3k4-revision',
				trackId: 'cistrome-track-3k4-revision',
				rowInfoAttributes: [
					{ field: 'Cell Type', type: 'nominal', position: 'right', width: 200 },
					{ field: 'Clustering', type: 'tree', position: 'right', width: 200 },
					{ field: 'PEAKS_TOTAL', type: 'quantitative', title: 'Total peaks', position: 'right', width: 100 },
					{
						field: 'SEQUENCE_LENGTH',
						type: 'quantitative',
						title: 'Sequence length',
						position: 'right',
						width: 100
					},
					{
						field: 'READS_MAPPED',
						type: 'quantitative',
						title: 'Reads mapped',
						position: 'right',
						width: 100
					}
				],
				rowAggregate: [],
				rowSort: [
					{ field: 'Clustering', type: 'tree', order: 'ascending' }
					// {field: "Cell Type", type: "nominal", order: "ascending"}
				],
				rowFilter: []
			}
		]
	}
};

export const miraDemos = {
	e18_mouse_brain_10x_top_500_smooth: {
		viewConfig: hgDemoViewConfigMiraMouse500Smooth,
		options: [
			{
				viewId: 'mira-view-mouse-500-smooth',
				trackId: 'mira-track-mouse-500-smooth',
				rowInfoAttributes: [
					// {
					//     field: "max_topic", type: "nominal", position: "right", width: 40,
					//         categories: [
					//             "atac_topic_11",
					//             "atac_topic_5",
					//             "atac_topic_6",
					//             "atac_topic_2",
					//             "atac_topic_8",
					//             "atac_topic_7",
					//             "atac_topic_0",
					//             "atac_topic_9",
					//             "atac_topic_10",
					//             "atac_topic_1",
					//             "atac_topic_12",
					//             "atac_topic_4",
					//             "atac_topic_3",
					//         ],
					//         colors: TOPIC_COLORS
					// },
					// { ...streamBase, width: 100 },
					// { ...streamBase, only: { field: 'tree_states', value: 'Astrocyte' }, title: 'Astrocyte' },
					// { ...streamBase, only: { field: 'tree_states', value: 'Astrocyte, Excitatory, Inhibitory' }, title: 'Astrocyte, Excitatory, Inhibitory' },
					// { ...streamBase, only: { field: 'tree_states', value: 'Inhibitory' }, title: 'Inhibitory' },
					// { ...streamBase, only: { field: 'tree_states', value: 'Excitatory, Inhibitory' }, title: 'Excitatory, Inhibitory' },
					// { ...streamBase, only: { field: 'tree_states', value: 'Excitatory' }, title: 'Excitatory' },
					// {
					//     field: "tree_states", type: "nominal", position: "right", width: 200, categories: STATES, colors: STATE_COLORS
					// },
					// { field: "Astrocyte", type: "nominal", position: "right", width: 30, categories: STATES, colors: STATE_COLORS, title: 'A'},
					// { field: "Astrocyte, Excitatory, Inhibitory", type: "nominal", position: "right", width: 30, categories: STATES, colors: STATE_COLORS, title: 'AEI'},
					// { field: "Inhibitory", type: "nominal", position: "right", width: 30, categories: STATES, colors: STATE_COLORS, title: 'I'},
					// { field: "Excitatory, Inhibitory", type: "nominal", position: "right", width: 30, categories: STATES, colors: STATE_COLORS, title: 'EI'},
					// { field: "Excitatory", type: "nominal", position: "right", width: 30, categories: STATES, colors: STATE_COLORS, title: 'E'},
					{
						field: 'mira_pseudotime',
						type: 'quantitative',
						position: 'right',
						width: 80,
						title: 'Pseudotime',
						color: 'grey'
					},
					{
						title: 'Tree State',
						field: 'tree_states',
						type: 'nominal',
						position: 'right',
						width: 80
					},
					{
						title: 'Lineage',
						field: 'tree_states',
						type: 'lineage',
						position: 'right',
						width: 180,
					},
					{title: "Topic 11", field: "atac_topic_11", type: 'quantitative', position: 'right', width: 60, color: '#FF9896' },
					{title: "Topic 5",  field: "atac_topic_5", type: 'quantitative', position: 'right', width: 60, color: '#2077B4' },
					{title: "Topic 6",  field: "atac_topic_6", type: 'quantitative', position: 'right', width: 60, color: '#7F7F7F' },
					{title: "Topic 2",  field: "atac_topic_2", type: 'quantitative', position: 'right', width: 60, color: '#AEC7E8' },
					{title: "Topic 8",  field: "atac_topic_8", type: 'quantitative', position: 'right', width: 60, color: '#BCBD22' },
					{title: "Topic 7",  field: "atac_topic_7", type: 'quantitative', position: 'right', width: 60, color: '#F7B6D2' },
					{title: "Topic 0",  field: "atac_topic_0", type: 'quantitative', position: 'right', width: 60, color: '#DBDB8E' },
					{title: "Topic 9",  field: "atac_topic_9", type: 'quantitative', position: 'right', width: 60, color: '#E377C2' },
					{title: "Topic 10", field: "atac_topic_10", type: 'quantitative', position: 'right', width: 60, color: '#C5B0D5' },
					{title: "Topic 1",  field: "atac_topic_1", type: 'quantitative', position: 'right', width: 60, color: '#D62828' },
					{title: "Topic 12", field: "atac_topic_12", type: 'quantitative', position: 'right', width: 60, color: '#2CA02C' },
					{title: "Topic 4",  field: "atac_topic_4", type: 'quantitative', position: 'right', width: 60, color: 'grey' },
					{title: "Topic 3",  field: "atac_topic_3", type: 'quantitative', position: 'right', width: 60, color: '#9FDBE5' },
					{
						field: 'index',
						url: 'https://s3.amazonaws.com/gosling-lang.org/data/cistrome/e18_mouse_brain_10x_rna_main.zarr/',
						type: 'expression',
						position: 'right',
						width: 80,
						title: 'XKR4'
					}
				],
				rowAggregate: [],
				rowSort: [
					{ field: 'tree_states', type: 'nominal', order: 'descending', sort: [
						'Astrocyte, Excitatory, Inhibitory',
						'Astrocyte',
						'Excitatory, Inhibitory',
						'Inhibitory',
						'Excitatory',
						'mira_pseudotime'
					] },
					{ field: 'mira_pseudotime', type: 'quantitative', order: 'ascending' }
				],
				rowFilter: []
			}
		]
	},
	e18_mouse_brain_10x_top_500: {
		viewConfig: hgDemoViewConfigMiraMouse,
		options: [
			{
				viewId: 'mira-view-mouse',
				trackId: 'mira-track-mouse',
				rowInfoAttributes: [
					{field: [
					    'topic_0',
					    'topic_1',
					    'topic_2',
					    'topic_3',
					    'topic_4',
					    'topic_5',
					    'topic_6',
					    'topic_7',
					    'topic_8',
					    'topic_9',
					    'topic_10',
					    'topic_11',
					    'topic_12',
					    'topic_13',
					    'topic_14',
					    'topic_15',
					    'topic_16',
					    'topic_17',
					    'topic_18',
					    'topic_19',
					    'topic_20',
					    'topic_21',
					], type: "quantitative", position: "right", width: 60},
					// { field: 'max_topic', type: 'nominal', position: 'right', width: 100 },
					// {field: "cluster_by_topic_0", type: "nominal", position: "right", width: 100},
					// {field: "cluster_by_topic_9", type: "nominal", position: "right", width: 100},
					{
						field: 'mira_pseudotime',
						type: 'quantitative',
						position: 'right',
						width: 80,
						title: 'Pseudotime',
						color: 'grey'
					},
					// { field: 'topic_0', type: 'quantitative', position: 'right', width: 60, color: '#FF9896' },
					// { field: 'topic_1', type: 'quantitative', position: 'right', width: 60, color: '#2077B4' },
					// { field: 'topic_2', type: 'quantitative', position: 'right', width: 60, color: '#7F7F7F' },
					// { field: 'topic_3', type: 'quantitative', position: 'right', width: 60, color: '#AEC7E8' },
					// { field: 'topic_4', type: 'quantitative', position: 'right', width: 60, color: '#BCBD22' },
					// { field: 'topic_5', type: 'quantitative', position: 'right', width: 60, color: '#F7B6D2' },
					// { field: 'topic_6', type: 'quantitative', position: 'right', width: 60, color: '#DBDB8E' },
					// { field: 'topic_7', type: 'quantitative', position: 'right', width: 60, color: '#E377C2' },
					// { field: 'topic_8', type: 'quantitative', position: 'right', width: 60, color: '#C5B0D5' },
					// { field: 'topic_9', type: 'quantitative', position: 'right', width: 60, color: '#D62828' },
					// { field: 'topic_10', type: 'quantitative', position: 'right', width: 60, color: '#2CA02C' },
					// // {field: "topic_11", type: "quantitative", position: "right", width: 60, color: 'grey'},
					// { field: 'topic_12', type: 'quantitative', position: 'right', width: 60, color: '#9FDBE5' }

					// {field: "atac_topic_11", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_5", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_6", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_2", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_8", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_7", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_0", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_9", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_10", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_1", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_12", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_4", type: 'quantitative', position: 'right', width: 60 },
					// {field: "atac_topic_3", type: 'quantitative', position: 'right', width: 60 },
				],
				rowAggregate: [],
				rowSort: [
					{ field: 'max_topic', type: 'nominal', order: 'ascending' },
					// {field: "topic_9", type: "quantitative", order: "descending"},
					{ field: 'mira_pseudotime', type: 'quantitative', order: 'ascending' }
				],
				rowFilter: []
			}
		]
	},
	e18_mouse_brain_10x_top_4000: {
		viewConfig: hgDemoViewConfigMiraMouse4000,
		options: [
			{
				viewId: 'mira-view-mouse-4000',
				trackId: 'mira-track-mouse-4000',
				rowInfoAttributes: [
					{
						field: 'mira_pseudotime',
						type: 'quantitative',
						position: 'right',
						width: 80,
						title: 'Pseudotime',
						color: 'grey'
					}
					// {field: ["topic_0", "topic_1", "topic_2", "topic_3", "topic_4", "topic_5", "topic_6", "topic_7", "topic_8", "topic_9", "topic_10", "topic_11", "topic_12"], type: "quantitative", position: "right", width: 60},
					// {field: "max_topic", type: "nominal", position: "right", width: 100},
					// {field: "cluster_by_topic_0", type: "nominal", position: "right", width: 100},
					// {field: "cluster_by_topic_9", type: "nominal", position: "right", width: 100},
					// {field: "topic_0", type: "quantitative", position: "right", width: 60},
					// {field: "topic_1", type: "quantitative", position: "right", width: 60},
					// {field: "topic_2", type: "quantitative", position: "right", width: 60},
					// {field: "topic_3", type: "quantitative", position: "right", width: 60},
					// {field: "topic_4", type: "quantitative", position: "right", width: 60},
					// {field: "topic_5", type: "quantitative", position: "right", width: 60},
					// {field: "topic_6", type: "quantitative", position: "right", width: 60},
					// {field: "topic_7", type: "quantitative", position: "right", width: 60},
					// {field: "topic_8", type: "quantitative", position: "right", width: 60},
					// {field: "topic_9", type: "quantitative", position: "right", width: 60},
					// {field: "topic_10", type: "quantitative", position: "right", width: 60},
					// {field: "topic_11", type: "quantitative", position: "right", width: 60},
					// {field: "topic_12", type: "quantitative", position: "right", width: 60},
				],
				rowAggregate: [
					//  {field: "cluster_by_topic_9", type: "nominal", notOneOf: []},
				],
				rowSort: [
					// {field: "max_topic", type: "nominal", order: "ascending"}
					// {field: "topic_9", type: "quantitative", order: "descending"},
					{ field: 'mira_pseudotime', type: 'quantitative', order: 'ascending' }
				],
				rowFilter: []
			}
		]
	}
};
