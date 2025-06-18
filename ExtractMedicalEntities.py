import json
from medcat.cat import CAT

# ✅ Load MedMentions model zip
cat = CAT.load_model_pack(r"C:\Users\ShriP\Downloads\medmen_wstatus_2021_oct.zip")

with open("publications.json", 'r') as file:
    papers = json.load(file)['publications']

output = []

for pub in papers:
    result = cat.get_entities(pub["abstract"])
    entities = []
    for k, v in result.items():
        # entities.append({
            #json formatted with the stats
        # })
        print(v)
        print()
        entities.append(v)
    output.append({
        "pmid": pub["pmid"],
        "title": pub["title"],
        "abstract": pub["abstract"],
        "entities": entities
    })

with open("MedCAT_Entities.json", "w") as f:
    json.dump({"topic": "sleep improves memory", "publications": output}, f, indent=2)
