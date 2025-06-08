import json

const_file_name = "AssignedEntities.json"

with open(const_file_name, 'r') as file:
    publications_file = json.load(file)

publications = publications_file['publications']
publication_pairs_list = []

for publication in publications:
    list_sentence_pairs = []

    for sentence in publication["sentence"]:
        pairs_list = []

        # Generate all unordered entity pairs inside the sentence
        entities = sentence['entities']
        for i in range(len(entities)):
            for j in range(i+1, len(entities)):
                entity_1 = entities[i]
                entity_2 = entities[j]

                # Skip exact duplicates
                if (entity_1['word'] == entity_2['word']) and (entity_1['entity_group'] == entity_2['entity_group']):
                    continue

                pair = {
                    "word1": entity_1['word'],
                    "entity_group1": entity_1['entity_group'],
                    "score1": entity_1['score'],
                    "word2": entity_2['word'],
                    "entity_group2": entity_2['entity_group'],
                    "score2": entity_2['score']
                }
                pairs_list.append(pair)

        sentence_pair = {
            "sentence": sentence['sentence'],
            "entity_pairs": pairs_list
        }
        list_sentence_pairs.append(sentence_pair)

    publication_pairs = {
        "pmid": publication['pmid'],
        "abstract": publication['abstract'],
        "sentence_pairs": list_sentence_pairs  # <-- keep sentence-level structure
    }
    publication_pairs_list.append(publication_pairs)

# final output structure
publications_assigned = {
    "topic": publications_file["topic"],
    "publications": publication_pairs_list
}

with open("Pairs.json", mode="w", encoding="utf-8") as write_file:
    json.dump(publications_assigned, write_file, indent=2)
