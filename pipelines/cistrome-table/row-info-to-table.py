#%%
import numpy as np
import pandas as pd
import json
import urllib
from io import StringIO
# %%
data_folders = [
    'cistrome-track-3k27',
    'cistrome-track-3k4',
    'cistrome-track-atac'
]
for data_folder in data_folders:
    print(data_folder)
    row_info = pd.read_json('../../src/demo/fakedata/' + data_folder + '/rowInfo.json')
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
    metadata = metadata.drop(columns=['other_ids'])
    # metadata.to_excel(f'./{data_folder}.xlsx') # No module named 'openpyxl'
    metadata.to_csv(f'./{data_folder}.csv', index=False)
# metadata
# %%
