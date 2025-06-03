import re
import string
import json
import unicodedata
from spellchecker import SpellChecker

# Use the pyspellchecker for general English correction
spell = SpellChecker(distance=1)

# Add biomedical-specific terms (this list should be expanded as you process more data)
biomedical_vocab = [
    "hypoxia", "intermittent", "sleep", "apnea", "physical", "activity", 
    "cognitive", "deficit", "brain", "memory", "exercise", "estrogen",
    "processing", "speed", "gaba", "older", "adults"
]

# Add biomedical terms to dictionary
spell.word_frequency.load_words(biomedical_vocab)

def normalize_text(word):
    # Step 1: Lowercase
    word = word.lower()

    # Step 2: Remove any weird unicode characters
    word = unicodedata.normalize('NFKD', word).encode('ascii', 'ignore').decode('utf-8', 'ignore')

    # Step 3: Remove non-ascii or non-alpha numeric chars except hyphen
    word = re.sub(r'[^a-z0-9\- ]+', '', word)

    # Step 4: Add space between lowercase-uppercase boundaries (handles merged words)
    word = re.sub(r'([a-z])([A-Z])', r'\1 \2', word)

    # Step 5: Add space between word-number boundaries
    word = re.sub(r'([a-zA-Z])(\d)', r'\1 \2', word)
    word = re.sub(r'(\d)([a-zA-Z])', r'\1 \2', word)

    # Step 6: Split fully merged words using known biomedical words
    for vocab_word in biomedical_vocab:
        word = re.sub(r'(?<!\w)(' + vocab_word + r')(?=[a-z])', vocab_word + ' ', word)

    # Step 7: Remove multiple spaces
    word = re.sub(r'\s+', ' ', word).strip()

    # Step 8: Correct spelling (word by word)
    tokens = word.split()
    corrected_tokens = []
    for token in tokens:
        corrected = spell.correction(token)
        if corrected is None:
            corrected_tokens.append(token)
        else:
            corrected_tokens.append(corrected)

    word = ' '.join(corrected_tokens)

    return word

# Example usage:
if __name__ == "__main__":
    # Load your filtered entities file:
    with open("FilteredPerPublication.json", "r") as f:
        data = json.load(f)

    # Apply normalization to all words
    for publication in data['publications']:
        for entity in publication['entities']:
            entity['word'] = normalize_text(entity['word'])

    # Save back
    with open("NormalizedEntities.json", "w") as f:
        json.dump(data, f, indent=2)

    print("Normalization complete!")
