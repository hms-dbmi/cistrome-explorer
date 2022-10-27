#%%
import numpy as np
import pandas as pd
import json
import urllib
from io import StringIO
import operator

# %%
data_folders = [
    ('table-H3K27ac-v1', 'cistrome-track-3k27', 'rowInfoRevision.json'),
    ('table-H3K4me3-v1', 'cistrome-track-3k4', 'rowInfoRevision.json'),
    # ('table-ATAC-v1', 'cistrome-track-atac', 'rowInfoRevision.json'),
    # ('table-ATAC-v0', 'cistrome-track-atac', 'rowInfo.json')
]
for (name, data_folder, file_path) in data_folders:
    print(data_folder)
    row_info = pd.read_json('../../src/demo/fakedata/' + data_folder + '/' + file_path)
    external_ids = ','.join(row_info.GSM.tolist())
    gsm_to_cell = dict(zip(row_info.GSM.tolist(), row_info['Cell Type'].tolist()))
    
    # get cids
    url = f'http://develop.cistrome.org/cistrome/samples?external_ids={external_ids}&format=json&limit=99999'
    samples = pd.read_json(url)
    sample_data = pd.DataFrame.from_dict(samples.samples.to_dict(), orient='index')
    cids = sample_data.id.tolist()

    # get metadata
    isFirst = True
    i = 0
    if 'ATAC' in name:
        for cid in cids:
            url = f'http://develop.cistrome.org/cistrome/samples/{cid}?fields=ontologies'
            row = pd.read_json(url, lines=True)
            samples = row.samples.apply(lambda x: x[0])
            row = pd.DataFrame.from_dict(samples.to_dict(), orient='index')
            celltype = pd.DataFrame.from_dict(row.ontologies.apply(lambda x: sorted(x, key=lambda k: k['value'])[-1] if len(x) >= 1 else None).to_dict(), orient='index')
            row['paper__pmid'] = 0
            row['paper__reference'] = ''
            try:
                row['ontologies'] = celltype['term']
            except:
                print(url + ' no cell type ')

            # publication
            url = f'http://dc2.cistrome.org/api/inspector?id={cid}'
            try:
                inspector = pd.read_json(url, lines=True)
                inspector.treats = inspector.treats.apply(lambda x: x[0])
                inspector = pd.DataFrame.from_dict(inspector.treats.reset_index().treats.to_dict(), orient='index')
                row['paper__reference'] = inspector['paper__reference']
                row['paper__pmid'] = inspector['paper__pmid']
            except:
                print(url + ' not working ')

            if isFirst:
                metadata = row
                isFirst = False
                print(cid, row)
            else:
                # break
                metadata = pd.concat([metadata, row])
        
         # Remove unused
        #  ID,Species,Factor,Cell Type,Tissue Type,Paper Reference,Paper PMID
        metadata = metadata.drop(columns=[
            'assembly',
            'experiment_type',
            'external_id_type',
            'sample_type',
            'title',
            'released',
            'factors',
            'qcs',
            'series',
            'id'
        ])
        
        # readable columns
        metadata = metadata.rename(columns={
            'species': 'Species',
            # 'factor__name': 'Factor',
            'external_id': 'ID',
            'ontologies': 'Cell Type',
            'paper__reference': 'Paper Reference',
            'paper__pmid': 'PMID'
        })
        metadata['PMID'] = metadata['PMID'].apply(lambda x: '0' if x == None else x)
        metadata['PMID'] = metadata['PMID'].astype(int)
        metadata['Factor'] = 'ATAC'
        metadata['Species'] = 'Homo sapiens'

        metadata = metadata[['ID', 'Species', 'Factor', 'Cell Type', 'Paper Reference', 'PMID']]
        metadata['Cell Type'] = metadata['ID']
        metadata['Cell Type'] = metadata['Cell Type'].apply(lambda x: gsm_to_cell[x])
        metadata.to_csv(f'./{name}.csv', index=False)
    else:
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
            'paper__pmid': 'PMID'
        })
        # print(metadata.PMID)
        metadata['PMID'] = metadata['PMID'].apply(lambda x: x if x == x else 0).astype(int)

        # metadata.to_excel(f'./{data_folder}.xlsx') # No module named 'openpyxl'
        metadata['Cell Type'] = metadata['ID']
        metadata['Cell Type'] = metadata['Cell Type'].apply(lambda x: gsm_to_cell[x])
        metadata.to_csv(f'./{name}.csv', index=False)
# metadata
# %%
import requests
import re

result = ''
for gsm in ['GSM2544227',
'GSM2544222',
'GSM2695562',
'GSM2695561',
'GSM2679879',
'GSM2679876',
'GSM2679870',
'GSM2371314',
'GSM2371311',
'GSM2679897',
'GSM2679891',
'GSM2679887',
'GSM2679885',
'GSM2679889',
'GSM2679867',
'GSM2371313',
'GSM2702716',
'GSM2702713',
'GSM2714170',
'GSM3452726',
'GSM3331038',
'GSM3331039',
'GSM3331034',
'GSM3331036',
'GSM3331033',
'GSM3320297',
'GSM3320295',
'GSM2829008',
'GSM3391515',
'GSM3391514',
'GSM3391517',
'GSM3391516',
'GSM3057381',
'GSM2817576',
'GSM3320256',
'GSM3320254',
'GSM3320255',
'GSM3320251',
'GSM3331032',
'GSM3320344',
'GSM3320340',
'GSM3320323',
'GSM3320328',
'GSM3320329',
'GSM2714162',
'GSM2714163',
'GSM2714168',
'GSM3320300',
'GSM3320302',
'GSM2688936',
'GSM2688934',
'GSM3320289',
'GSM3320283',
'GSM3320282',
'GSM3320343',
'GSM3320380',
'GSM3320381',
'GSM2829010',
'GSM2689786',
'GSM3320357',
'GSM3320258',
'GSM3320259',
'GSM2688944',
'GSM2688941',
'GSM3320249',
'GSM3320248',
'GSM3320241',
'GSM3320240',
'GSM3320322',
'GSM3320261',
'GSM3320371',
'GSM3320358',
'GSM3320334',
'GSM3320331',
'GSM3320330',
'GSM3320339',
'GSM2564056',
'GSM2564055',
'GSM2564045',
'GSM2564038',
'GSM2564037',
'GSM2471276',
'GSM2471294',
'GSM3137778',
'GSM2471255',
'GSM3146466',
'GSM3701911',
'GSM2945825',
'GSM3032066',
'GSM3032065',
'GSM3032064',
'GSM3020373',
'GSM3020336',
'GSM3020332',
'GSM3020327',
'GSM3020325',
'GSM3020323',
'GSM3020298',
'GSM3020280',
'GSM3020259',
'GSM3020226',
'GSM3716514',
'GSM3716513',
'GSM3597647',
'GSM3597646',
'GSM3106831',
'GSM3312795',
'GSM3770763',
'GSM3770758',
'GSM3693107',
'GSM3444160',
'GSM3444162',
'GSM3444163',
'GSM3770776',
'GSM3873133',
'GSM4025770',
'GSM4025765',
'GSM4051648',
'GSM4191897',
'GSM4191896',
'GSM4191895',
'GSM3231157',
'GSM3321011',
'GSM3462844',
'GSM3067777',
'GSM3067775',
'GSM3067769',
'GSM4399738',
'GSM4332570',
'GSM4332548',
'GSM4332560',
'GSM4332558',
'GSM4332553',
'GSM4332549',
'GSM4332572',
'GSM4332574',
'GSM4332556',
'GSM4332562',
'GSM4332575',
'GSM4332557',
'GSM4674818',
'GSM4767215',
'GSM4767222',
'GSM4767225',
'GSM4426262',
'GSM4570071',
'GSM3738858',
'GSM2650673',
'GSM4565972',
'GSM4565974',
'GSM4565965',
'GSM4486191',
'GSM4162161',
'GSM4162164',
'GSM4162165',
'GSM4162166',
'GSM4520713',
'GSM4314811',
'GSM4648801',
'GSM4648806',
'GSM4648809',
'GSM4648810',
'GSM4648811',
'GSM4648812',
'GSM4648813',
'GSM4648816',
'GSM4648817',
'GSM4648818',
'GSM4648819',
'GSM4648820',
'GSM4648821',
'GSM4837479',
'GSM4299096',
'GSM4299098',
'GSM4299082',
'GSM4299092',
'GSM4299107',
'GSM3900498',
'GSM4764093',
'GSM4764100',
'GSM4764102',
'GSM4764103',
'GSM4764104',
'GSM4764106',
'GSM4291167',
'GSM4289909',
'GSM4289910',
'GSM4289912',
'GSM4289915',
'GSM4289918',
'GSM4289919',
'GSM4138221',
'GSM4138222',
'GSM4139429',
'GSM4139469',
'GSM4259632',
'GSM4259633',
'GSM4259634',
'GSM4259635',
'GSM4259645',
'GSM3962425',
'GSM3813957',
'GSM3630261',
'GSM4904185',
'GSM4904183',
'GSM4904184',
'GSM4904182',
'GSM4904187',
'GSM4904181',
'GSM4904180',
'GSM4542936',
'GSM4542937',
'GSM4542938',
'GSM3632980',
'GSM3632981',
'GSM3632982',
'GSM3632985',
'GSM4801596',
'GSM4610440',
'GSM4625592',
'GSM4625610',
'GSM3584949',
'GSM4306144',
'GSM4306142',
'GSM4306149',
'GSM4306153',
'GSM4306145',
'GSM4306152',
'GSM4306148',
'GSM4306143',
'GSM4306151',
'GSM4306140',
'GSM4306141',
'GSM3302065',
'GSM4043955',
'GSM4409114',
'GSM4409124',
'GSM4340223',
'GSM4340228',
'GSM4340224',
'GSM4340227',
'GSM4139397',
'GSM4139399',
'GSM4139412',
'GSM4139413',
'GSM4139417',
'GSM4139433',
'GSM4139460',
'GSM4139473',
'GSM4565966',
'GSM3738846',
'GSM3738857',
'GSM3738859',
'GSM3738871',
'GSM3738872',
'GSM4570072',
'GSM4110238',
'GSM4520712',
'GSM4674823',
'GSM4289908',
'GSM4764099',
'GSM4054640',
'GSM3738852',
'GSM3738854',
'GSM3738860',
'GSM3738864',
'GSM3738879',
'GSM4570064',
'GSM4674819',
'GSM4767218',
'GSM4767220',
'GSM4767223',
'GSM4299086',
'GSM4299108',
'GSM4299100',
'GSM4332567',
'GSM4837481',
'GSM4837487',
'GSM4156591',
'GSM4156595',
'GSM4295137',
'GSM4295141',
'GSM4301065',
'GSM4289916',
'GSM4764101',
'GSM4064961',
'GSM3962404',
'GSM3962407',
'GSM3962412',
'GSM3962414',
'GSM3962418',
'GSM4625594',
'GSM4625587',
'GSM4625578',
'GSM3738372',
'GSM3738373',
'GSM3738371',
'GSM4054634',
'GSM4054638',
'GSM3584961',
'GSM3780773',
'GSM3780775',
'GSM3780772',
'GSM3780768',
'GSM3780771',
'GSM3780770',
'GSM3780767',
'GSM4818345',
'GSM4818343',
'GSM4818346',
'GSM4409122',
'GSM3632987',
'GSM4340221']:
    try:
        gsm_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"
        gsm_page = requests.get(gsm_url)
        gse = re.search('acc=GSE(.*)" onmouseout', gsm_page.text).group(1)
        # print(gse)

        gse_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE{gse}"
        gse_page = requests.get(gse_url)
        pmid = re.search('pubmed_id" id="(.*)"><a', gse_page.text).group(1)
        # print(pmid)

        # url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        # page = requests.get(url)
        # citation_authors = re.search('citation_authors" content="(.*)">', page.text).group(1)
        # citation_authors = citation_authors.replace(';', ', ')
        # citation_title = re.search('citation_title" content="(.*)">', page.text).group(1)
        # citation_journal_title = re.search('citation_journal_title" content="(.*)">', page.text).group(1)
        # citation = citation_authors + citation_title + ', ' + citation_journal_title

        result += pmid + '\n'
    except:
        result += 'Missing' + '\n'
print(result)
# %%
