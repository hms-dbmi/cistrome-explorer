#%%
import numpy as np
import pandas as pd
import json
import urllib
from io import StringIO
# %%
data_folders = [
    ('table-H3K27ac (v1)', 'cistrome-track-3k27', 'rowInfoRevision.json'),
    ('table-H3K4me3 (v1)', 'cistrome-track-3k4', 'rowInfoRevision.json'),
    ('table-ATAC (v1)', 'cistrome-track-atac', 'rowInfoRevision.json'),
    ('table-ATAC (v0)', 'cistrome-track-atac', 'rowInfo.json')
]
for (name, data_folder, file_path) in data_folders:
    print(data_folder)
    row_info = pd.read_json('../../src/demo/fakedata/' + data_folder + '/' + file_path)
    external_ids = ','.join(row_info.GSM.tolist())
    
    # get cids
    url = f'http://develop.cistrome.org/cistrome/samples?external_ids={external_ids}&format=json&limit=99999'
    samples = pd.read_json(url)
    sample_data = pd.DataFrame.from_dict(samples.samples.to_dict(), orient='index')
    cids = sample_data.id.tolist()

    # get metadata
    isFirst = True
    i = 0
    for cid in cids:
        url = f'http://dc2.cistrome.org/api/inspector?id={cid}'
        # i += 1
        # if i == 3:
        #     break
        try:
            row = pd.read_json(url, lines=True)
            if isFirst:
                metadata = row
                isFirst = False
            else:
                metadata = pd.concat([metadata, row])
        except:
            print(url + ' not working ')

    metadata.treats = metadata.treats.apply(lambda x: x[0])
    metadata = pd.DataFrame.from_dict(metadata.treats.reset_index().treats.to_dict(), orient='index')
    
    # Remove unused
    metadata = metadata.drop(columns=[
        'other_ids',
        'cell_line__name',
        'is_correcting',
        'strain__name',
        'cell_pop__name',
        'paper__journal__name',
        'name',
        'disease_state__name',
        'link',
        'paper__lab'
    ])

    # sort rows
    metadata = metadata.sort_values(by=['unique_id'])

    # sort columns
    metadata = metadata[['unique_id', 'species__name', 'factor__name', 'cell_type__name', 'tissue_type__name', 'paper__reference', 'paper__pmid']]

    # readable columns
    metadata = metadata.rename(columns={
        'unique_id': 'ID',
        'species__name': 'Species',
        'factor__name': 'Factor',
        'cell_type__name': 'Cell Type',
        'tissue_type__name': 'Tissue Type',
        'paper__reference': 'Paper Reference',
        'paper__pmid': 'Paper PMID'
    })

    # metadata.to_excel(f'./{data_folder}.xlsx') # No module named 'openpyxl'
    metadata.to_csv(f'./{name}.csv', index=False)
# metadata
# %%
