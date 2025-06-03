import json
import string

const_file_name = "NormalizedEntities.json"

with open(const_file_name, 'r') as file:
    file_normalized_entities = json.load(file)

publication_normalized_entities = file_normalized_entities['publications']

publication_pairs_list = []

for publications in publication_normalized_entities:
    pairs_list = []
    for i in range(len(publications['entities'])):
        for j in range(i+1, len(publications['entities'])):
            entity_1 = publications['entities'][i]
            entity_2 = publications['entities'][j]
            if (entity_1['word'] == entity_2['word'])  and (entity_1['entity_group'] == entity_2['entity_group']):
                continue
            else:
                pair = {
                        "word1": entity_1['word'],
                        "entity_group1": entity_1['entity_group'],
                        "score1": entity_1['score'],
                        "word2": entity_2['word'],
                        "entity_group2": entity_2['entity_group'],
                        "score2": entity_2['score']
                    }
                pairs_list.append(pair)
    publication_pairs = {
        "pmid": publications['pmid'],
        "title": publications['title'],
        "abstract": publications['abstract'],
        "entity_pairs": pairs_list
    }
    publication_pairs_list.append(publication_pairs)
    

medical_entity_pairs = {
    "topic": file_normalized_entities['topic'],
    "publications": publication_pairs_list
}

with open("UnorderedPairs.json", mode="w", encoding="utf-8") as write_file:
    json.dump(medical_entity_pairs, write_file, indent = 2)
# print(medical_entity_pairs)
    

                