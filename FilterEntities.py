import json
from itertools import chain

file_name = 'medicalEntities.json'

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
            if prev_word is None and cur_word.startswith("##"):
                cur_word = cur_word.lstrip("#")
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

merged_words_entity_groups = []
threshold_gap = 4
for publications in concatenated:
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
        if prev_word != None and prev_end >= cur_start - threshold_gap and prev_end <= cur_start:
            print("JOINED")
            combine_word = prev_word + " " + cur_word
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
    merged_words_entity_groups.append(new_pmd_entities)
    
from collections import defaultdict

entity_aggregate = defaultdict(lambda: {"count": 0, "score_sum": 0.0, "entity_group": None})

for pub in merged_words_entity_groups:
    for e in pub['entities']:
        key = (e['word'], e['entity_group'])
        entity_aggregate[key]['count'] += 1
        entity_aggregate[key]['score_sum'] += float(e['score'])
        entity_aggregate[key]['entity_group'] = e['entity_group']

final_kg_entities = []
for (word, entity_group), data in entity_aggregate.items():
    avg_score = data['score_sum'] / data['count']
    final_kg_entities.append({
        "word": word,
        "entity_group": entity_group,
        "avg_score": avg_score,
        "count": data['count']
    })

print("distinct kg entities")
print(final_kg_entities)

            

