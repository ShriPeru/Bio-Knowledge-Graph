import json

const_file_name = "ExtractedRelations.json"

with open(const_file_name, 'r') as file:
    publications_file = json.load(file)

publications = publications_file

filtered_list = []
conf_threshold = 0.5
for relation in publications:
    if float(relation["confidence"]) >= conf_threshold:
        triple = {
            "entity1": relation["entity1"],
            "relation": relation["relation"],
            "entity2": relation["entity2"]
        }
        filtered_list.append(triple)



with open("KG_filtered.json", mode="w", encoding="utf-8") as write_file:
    json.dump(filtered_list, write_file, indent=2)
