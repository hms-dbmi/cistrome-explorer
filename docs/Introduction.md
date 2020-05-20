[Cistrome Explorer](http://cisvis.gehlenborglab.org) is an interactive visualization tool to explore and analyze [cistrome data](http://cistrome.org/db/#/).
Cistrome Explorer uses [HiGlass](https://higlass.io) to show cistrome data in heatmaps while providing HiGlass wrapper components for showing metadata, such as dendromgram for hierarchical clustering results and bar charts for quality scores. These wrapper components are highly interactive which allow ones to rearrange, filter, and highlight rows in HiGlass heatmaps. Take a look at the [demo](http://cisvis.gehlenborglab.org) or [source codes](https://github.com/hms-dbmi/cistrome-higlass-wrapper) in GitHub.

[Cistrome Explorer](http://cisvis.gehlenborglab.org) support two types of wrapper components: *sample-wise* and *genome-wise* components.

Sample-wise components: 

1. **Nominal Bar Chart**: A nominal data field, such as a cell type, can be visualized using a nominal bar chart that shows text labels with different colors for individual categories.
1. **Quantitative Bar Chart**: One or more quantitative data fields, such as quality scores for individual samples, can be represented as a horizontal (stacked) bar chart.
1. **Dendrogram**: Tree-structured data, such as hierarchical clustering results, can be visualized using a Dendrogram track. You can filter or highlight rows by selecting a branch or moving a minimum similarity bar using the mouse cursor.

1. **Link**: Useful external links for individual samples can be shown with a Link track where clicking on a link on this track opens a linked website.

Genome-wise component: 

1. **Genomic Interval Selection Bar**: To support querying for [Cistrome DB](http://cistrome.org/db/#/), a genomic interval selection bar can be displayed on the bottom of each HiGlass heatmap. On this component, you can select multiple regions of interest and query for a data table that shows factors bound in these regions. To learn more about the query, please refer to [Cistrome ToolKit](http://dbtoolkit.cistrome.org).

You can set up the wrapper components using the JSON format specification. 
For example, the following specification generates four tracks on the left and right of a HiGlass heatmap with two sorting and one filtering options:

```javascript
{
    rowInfoAttributes: [
        {field: "Hierarchical Clustering (Average)", type: "tree", position: "left"},
        {field: "qc_frip", type: "quantitative", position: "left"},
        {field: "qc_fastqc", type: "quantitative", position: "left"},
        {field: "Metadata URL", type: "url", position: "left", title: "cid"},
        {field: "Hierarchical Clustering (Ward)", type: "tree", position: "right"},
        {field: "Cell Type", type: "nominal", position: "right"},
        {field: "Tissue Type", type: "nominal", position: "right"},
        {field: "Species", type: "nominal", position: "right"}
    ],
    rowSort: [
        {field: "Tissue Type", type: "nominal", order: "ascending"},
        {field: "qc_frip", type: "quantitative", order: "descending"}
    ],
    rowFilter: [
        {field: "Tissue Type", type: "nominal", notOneOf: ["None"]}
    ]
}
```



Please refer to the definition of the [JSON schema](https://github.com/hms-dbmi/cistrome-higlass-wrapper/blob/e1f9d2e83fa7af684f6cb827b8f7aae92a6f6b8a/src/utils/options.js#L9).