import json

with open("FilteredMEDCAT.json", 'r') as file:
    file_normalized_entities = json.load(file)

publication_normalized_entities = file_normalized_entities['publications']
publication_pairs_list = []

for publication in publication_normalized_entities:
    entities = publication['entities']
    pairs_list = []
    for i in range(len(entities)):
        for j in range(i + 1, len(entities)):
            e1 = entities[i]
            e2 = entities[j]

            # Skip identical concepts (same CUI and type)
            if (e1['cui'] == e2['cui']) and (e1['types'] == e2['types']):
                continue

            pair = {
                "entity1": e1["pretty_name"],
                "cui1": e1["cui"],
                "type1": e1["types"][0] if e1["types"] else "",
                "acc1": e1["acc"],

                "entity2": e2["pretty_name"],
                "cui2": e2["cui"],
                "type2": e2["types"][0] if e2["types"] else "",
                "acc2": e2["acc"]
            }

            pairs_list.append(pair)

    publication_pairs = {
        "pmid": publication['pmid'],
        "title": publication['title'],
        "abstract": publication['abstract'],
        "entity_pairs": pairs_list
    }

    publication_pairs_list.append(publication_pairs)

medical_entity_pairs = {
    "topic": file_normalized_entities['topic'],
    "publications": publication_pairs_list
}

with open("UnorderedPairs.json", mode="w", encoding="utf-8") as write_file:
    json.dump(medical_entity_pairs, write_file, indent=2)


                