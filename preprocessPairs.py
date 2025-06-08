import json

# Load Pairs.json
with open('Pairs.json', 'r') as f:
    data = json.load(f)

extracted_pairs = []

for publication in data['publications']:
    pmid = publication['pmid']

    for sentence_data in publication['sentence_pairs']:
        sentence_text = sentence_data['sentence']

        for pair in sentence_data['entity_pairs']:
            e1_text = pair['word1']
            e1_type = pair['entity_group1']

            e2_text = pair['word2']
            e2_type = pair['entity_group2']

            extracted_pairs.append({
                "pmid": pmid,
                "sentence": sentence_text,
                "entity1": e1_text,
                "entity1_type": e1_type,
                "entity2": e2_text,
                "entity2_type": e2_type
            })

# Save the formatted pairs to feed into BioRED model
with open("FormattedForBioRED.json", "w") as f:
    json.dump(extracted_pairs, f, indent=2)
