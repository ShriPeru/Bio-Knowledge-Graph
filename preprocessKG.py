import json

const_file_name = "ExtractedRelations.json"

with open(const_file_name, 'r') as file:
    publications_file = json.load(file)

pairs = publications_file

filtered_list = []
seen_triples = set()
conf_threshold = 0.3

for pair in pairs:
    if float(pair["confidence"]) >= conf_threshold:
        triple_tuple = (pair["entity1"], pair["relation"], pair["entity2"])
        if triple_tuple not in seen_triples:
            seen_triples.add(triple_tuple)
            filtered_list.append({
                "entity1": pair["entity1"],
                "relation": pair["relation"],
                "entity2": pair["entity2"]
            })

with open("KG_filtered.json", mode="w", encoding="utf-8") as write_file:
    json.dump(filtered_list, write_file, indent=2)

import pandas as pd

with open('KG_filtered.json', encoding='utf-8') as inputfile:
    df = pd.read_json(inputfile)

df.to_csv('KG_relations.csv', encoding='utf-8', index=False)
