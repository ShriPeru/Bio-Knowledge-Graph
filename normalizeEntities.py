import re
import string
import json
import unicodedata
from spellchecker import SpellChecker
from flashtext import KeywordProcessor

# Use the pyspellchecker for general English correction (backup spell checker)
spell = SpellChecker(distance=1)

# Expanded biomedical vocabulary (add more phrases as you go)
biomedical_vocab = [
    "hypoxia", "intermittent hypoxia", "sleep", "apnea", "physical activity", 
    "moderate vigorous physical activity", "light intensity physical activity",
    "cognitive", "cognitive deficit", "memory", "exercise", "estrogen",
    "processing speed", "gaba", "older adults"
]

# For flashtext phrase matching
keyword_processor = KeywordProcessor(case_sensitive=False)
for phrase in biomedical_vocab:
    keyword_processor.add_keyword(phrase)

# Add biomedical individual words to spellchecker
for word in ' '.join(biomedical_vocab).split():
    spell.word_frequency.add(word.lower())

def normalize_text(word):
    # Step 1: Lowercase
    word = word.lower()

    # Step 2: Remove weird unicode characters
    word = unicodedata.normalize('NFKD', word).encode('ascii', 'ignore').decode('utf-8', 'ignore')

    # Step 3: Remove non-alphanumeric except hyphen
    word = re.sub(r'[^a-z0-9\- ]+', '', word)

    # Step 4: Handle camel-case / digit splits (basic prep)
    word = re.sub(r'([a-z])([A-Z])', r'\1 \2', word)
    word = re.sub(r'([a-zA-Z])(\d)', r'\1 \2', word)
    word = re.sub(r'(\d)([a-zA-Z])', r'\1 \2', word)

    # Step 5: Run phrase-based matching using flashtext
    word = keyword_processor.replace_keywords(word)

    # Step 6: Clean extra spaces
    word = re.sub(r'\s+', ' ', word).strip()

    # Step 7: Backup spell correction (word-by-word)
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

# Actual pipeline
if __name__ == "__main__":
    with open("FilteredPerPublication.json", "r") as f:
        data = json.load(f)

    for publication in data['publications']:
        for entity in publication['entities']:
            entity['word'] = normalize_text(entity['word'])

    with open("NormalizedEntities.json", "w") as f:
        json.dump(data, f, indent=2)

    print("Normalization complete!")
