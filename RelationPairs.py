import json

const_file_name = "AssignedEntities.json"

with open(const_file_name, 'r') as file:
    publications_file = json.load(file)

publications = publications_file['publications']
publication_pairs_list = []

for publication in publications:
    pairs_list = []

    for sentence in publication.get('sentence', []):
        entities = sentence.get('entities', [])
        if len(entities) < 2:
            continue  # skip if fewer than 2 entities

        for i in range(len(entities)):
            for j in range(i + 1, len(entities)):
                entity_1 = entities[i]
                entity_2 = entities[j]

                if entity_1['cui'] == entity_2['cui']:
                    continue  # skip identical CUIs

                pair = {
                    "entity1": entity_1,
                    "entity2": entity_2,
                    "sentence": sentence['sentence']
                }
                pairs_list.append(pair)

    publication_pairs = {
        "pmid": publication['pmid'],
        "abstract": publication['abstract'],
        "entity_pairs": pairs_list
    }

    publication_pairs_list.append(publication_pairs)

# Wrap final output under topic
medical_entity_pairs = {
    "topic": publications_file.get('topic', ""),
    "publications": publication_pairs_list
}

# Save output
with open("UnorderedAssignedPairs.json", mode="w", encoding="utf-8") as write_file:
    json.dump(medical_entity_pairs, write_file, indent=2)
