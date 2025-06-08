import json
from transformers import pipeline

# Load your Pairs.json file
with open("Pairs.json", "r") as f:
    data = json.load(f)

# Initialize zero-shot classifier
classifier = pipeline("zero-shot-classification", model="facebook/bart-large-mnli")

# Define some simple relation labels (we keep this very simple)
labels = ["Cause-Effect", "Treatment", "Association", "No-Relation"]

results = []

for publication in data['publications']:
    pmid = publication['pmid']

    for sentence_data in publication['sentence_pairs']:
        sentence_text = sentence_data['sentence']

        for pair in sentence_data['entity_pairs']:
            e1 = pair['word1']
            e2 = pair['word2']

            # Create a simple hypothesis for zero-shot learning
            hypothesis = f"There is a relation between {e1} and {e2}."

            # Run classification
            output = classifier(sentence_text, labels)

            # Take top prediction
            pred_label = output['labels'][0]
            score = output['scores'][0]
            results.append({
                "pmid": pmid,
                "sentence": sentence_text,
                "entity1": e1,
                "entity2": e2,
                "relation": pred_label,
                "confidence": score
            })

# Save results
with open("ExtractedRelations.json", "w") as f:
    json.dump(results, f, indent=2)
