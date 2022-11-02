# %%
import json
# %%
with open('./rowInfo.json', 'r') as f:
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

with open('./rowInfoWithCategory.json', 'w') as f:
    json.dump(rowInfo, f)
# %%
