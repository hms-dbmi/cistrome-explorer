# %%
import json
import math
# %%
path = './rowInfo4000.json'
# %%
with open(path, 'r') as f:
    rowInfo = json.load(f)
    
    for row in rowInfo:
        maxValue = row['topic_0']
        maxIndex = 0
        for i in range(13):
            curValue = row[f'topic_{i}']
            if maxValue < curValue:
                maxIndex = i
                maxValue = curValue
        row['max_topic'] = f'topic_{maxIndex}' if maxIndex >= 10 else f'topic_0{maxIndex}'

        for c in ['topic_0', 'topic_9']:
            clusterValue = row[c]
            row[f'cluster_by_{c}'] = f'Cluster {math.floor(clusterValue * 2 / 100) + 1}'

with open(path, 'w') as f:
    json.dump(rowInfo, f)
# %%
