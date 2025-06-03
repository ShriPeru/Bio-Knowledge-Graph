import json
from itertools import chain

file_name = 'rawMedicalEntities.json'

with open(file_name, 'r') as file:
    rawMedicalEntities = json.load(file)
# print(paper_json)
store = rawMedicalEntities['publications']

concatenated = []
for publications in store:
    prev_word = None
    prev_end = None
    prev_start = None
    prev_score = None
    merged_num = 1
    prev_entity_group = None
    cur_entities = []
    for entity in publications['entities']:
        cur_end = entity['end']

        cur_word = entity['word']
        cur_start = entity['start']
        cur_score = entity['score']
        cur_entity_group = entity['entity_group']
        if prev_word != None and prev_end == cur_start and prev_entity_group == cur_entity_group:
            prev_word_iter = iter(prev_word)
            print("JOINED")
            combine_word = prev_word + cur_word.lstrip('#')
            combine_word = combine_word.replace('#','')
            merged_num+=1
            print(combine_word)
            prev_score += float(cur_score)
            prev_end = cur_end
            prev_word = combine_word
        else:
            if prev_word is not None:
                new_score = float(prev_score) / merged_num
                entity = {
                    "entity_group": prev_entity_group,
                    "score": str(new_score),
                    "word": prev_word,
                    "start": prev_start,
                    "end": prev_end
                }
                cur_entities.append(entity)
                
            prev_word = cur_word
            prev_end = cur_end
            prev_start = cur_start
            prev_score = float(cur_score)
            prev_entity_group = cur_entity_group
            merged_num = 1

    if prev_word is not None:
        new_score = float(prev_score) / merged_num
        cur_entities.append({
            "entity_group": prev_entity_group,
            "score": str(new_score),
            "word": prev_word,
            "start": prev_start,
            "end": prev_end
    })
    new_pmd_entities = {
        "pmid": publications['pmid'],
        "title": publications['title'],
        "abstract": publications['abstract'],
        "entities": cur_entities
    }
    concatenated.append(new_pmd_entities)


print(concatenated)
    
            

