import json
from transformers import pipeline

# Load UnorderedAssignedPairs.json
with open("UnorderedAssignedPairs.json", "r", encoding="utf-8") as f:
    data = json.load(f)

# Zero-shot classifier with BioLinkBERT
classifier = pipeline("zero-shot-classification", model="facebook/bart-large-mnli")

labels = ["Cause-Effect", "Treats", "Prevents", "Association", "No-Relation"]

results = []
for pub in data["publications"]:
    pmid = pub["pmid"]
    title = pub.get("title", "")
    for pair in pub["entity_pairs"]:
        sentence = pair.get("sentence", pub["abstract"])
        e1 = pair["entity1"]["pretty_name"]
        e2 = pair["entity2"]["pretty_name"]

        hypothesis = f"There is a {e1}–{e2} relation in this sentence."

        output = classifier(sentence, labels)
        label, score = output["labels"][0], output["scores"][0]

        results.append({
            "pmid": pmid,
            "title": title,
            "sentence": sentence,
            "entity1": e1,
            "entity2": e2,
            "relation": label,
            "confidence": round(score, 4)
        })

with open("ExtractedRelations.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2)
