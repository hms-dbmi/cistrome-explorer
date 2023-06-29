import { useEffect, useState } from "react";
import { AnnDataSource, ObsFeatureMatrixAnndataLoader } from '@vitessce/zarr';
import './GeneExpressionSelection.scss';

export default function GeneExpressionSelection(props) {
    const {
        left,
        top,
        width,
        height,
        url = 'https://s3.amazonaws.com/gosling-lang.org/data/cistrome/e18_mouse_brain_10x_rna_main.zarr/',
        genes,
        onGeneSelection = (genes) => { }
    } = props;

    const [metadata, setMetadata] = useState({});
    const [keyword, setKeyword] = useState('');
    const [selectedGenes, setSelectedGenes] = useState([...genes]);

    async function getMetadata() {
		const source = new AnnDataSource({ url });
		const config = {
			url,
			fileType: 'obsFeatureMatrix.mudata.zarr',
			options: {
				path: 'X'
			}
		};
		const loader = new ObsFeatureMatrixAnndataLoader(source, config);

		// obsIndex is cell IDs. varIndex is gene IDs.
		const {
			data: { rows: obsIndex, cols: varIndex }
		} = await loader.loadAttrs();

        setMetadata({
            obsIndex,
            varIndex,
            loader
        });
		// We can load the data for a subset of genes by selecting an array of gene IDs.
		// const { data } = await loader.loadGeneSelection({ selection: ['XKR4'] });
		// const expression = obsIndex.map((cellId, i) => ({ cellId, normalizedExpression: data[0][i] / 256 }));
		// callback(expression);
	}

    useEffect(() => {
        getMetadata();
    }, [url]);

    useEffect(() => {
        onGeneSelection(selectedGenes);
    }, [selectedGenes]);

    return (
        <div
            style={{
                position: 'absolute',
                top: `${top}px`,
                left: `${left}px`,
                width: `${width}px`,
                height: `${height}px`
            }}
        >
            <input
                id='ges-available-search'
                style={{ 
                    width: `${width/2.0}px`,
                    height: `30px`,
                    border: '1px solid black',
                }}
                onChange={e => setKeyword(e.target.value)}
                placeholder="Search for genes"
            />
            <div
                id='ges-available-genes'
                style={{ 
                    width: `${width/2.0}px`,
                    height: `calc(100% - ${30}px)`,
                    border: '1px solid black',
                    overflowY: 'scroll',
                }}
            >
                <ul>
                    {metadata.varIndex?.filter(gene => gene.toUpperCase().includes(keyword.toUpperCase())).sort().slice(0, 30).map(gene => {
                        return (
                            <li 
                                key={gene}
                                onClick={() => setSelectedGenes([...selectedGenes, gene])}
                            >
                                {gene}
                            </li>
                        );
                    })}
                </ul>
            </div>
            <textarea
                id='ges-selected-genes'
                style={{
                    width: `${width/2.0 - 6}px`,
                    height: `calc(100%)`,
                }}
                value={selectedGenes.join('\n')}
                onChange={e => setSelectedGenes(e.target.value.split('\n'))}
            />
        </div>
    )
}