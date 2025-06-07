import json
import re
from collections import defaultdict


def recompute_offsets(publications):
    for pub in publications:
        abstract = pub['abstract']
        for entity in pub['entities']:
            word = entity['word']
            original_start = entity['start']  # this is from your merged data, used as approximate anchor

            # Find all occurrences of the word in the abstract
            occurrences = [m.start() for m in re.finditer(re.escape(word), abstract)]

            if not occurrences:
                print(f"Warning: '{word}' not found in abstract for PMID {pub['pmid']}")
                continue

            # Find the occurrence closest to original_start
            best_match = min(occurrences, key=lambda x: abs(x - original_start))

            entity['start'] = best_match
            entity['end'] = best_match + len(word)


file_name = 'medicalEntities.json'

# load json
with open(file_name, 'r') as file:
    rawMedicalEntities = json.load(file)
store = rawMedicalEntities['publications']

# First pass: merge subwords (## based merging)
concatenated = []
for publications in store:
    prev_word = None
    prev_end = None
    prev_start = None
    prev_score = None
    merged_num = 1
    prev_entity_group = None
    prev_is_subword = False
    cur_entities = []

    for entity in publications['entities']:
        cur_end = entity['end']
        cur_start = entity['start']
        cur_score = float(entity['score'])
        cur_entity_group = entity['entity_group']

        word_raw = entity['word']
        cur_is_subword = word_raw.startswith("##")
        cur_word = word_raw.replace("#", "") if cur_is_subword else word_raw

        if prev_word is not None and prev_end == cur_start and prev_entity_group == cur_entity_group:
            if prev_is_subword or cur_is_subword:
                combine_word = (prev_word + cur_word).replace(" ", "")
            else:
                combine_word = prev_word + " " + cur_word

            merged_num += 1
            prev_score += cur_score
            prev_end = cur_end
            prev_word = combine_word
            prev_is_subword = prev_is_subword or cur_is_subword
        else:
            if prev_word is not None:
                new_score = prev_score / merged_num
                cur_entities.append({
                    "entity_group": prev_entity_group,
                    "score": str(new_score),
                    "word": prev_word,
                    "start": prev_start,
                    "end": prev_end
                })

            prev_word = cur_word
            prev_end = cur_end
            prev_start = cur_start
            prev_score = cur_score
            prev_entity_group = cur_entity_group
            prev_is_subword = cur_is_subword
            merged_num = 1

    if prev_word is not None:
        new_score = prev_score / merged_num
        cur_entities.append({
            "entity_group": prev_entity_group,
            "score": str(new_score),
            "word": prev_word,
            "start": prev_start,
            "end": prev_end
        })

    concatenated.append({
        "pmid": publications['pmid'],
        "title": publications['title'],
        "abstract": publications['abstract'],
        "entities": cur_entities
    })

# Second pass: merge broken tokens even when entity_group differs
merged_words_entity_groups = []
threshold_gap = 4  # adjustable hyperparameter

for publications in concatenated:
    prev_word = None
    prev_end = None
    prev_start = None
    prev_score = None
    prev_entity_group = None
    merged_num = 1
    cur_entities = []

    for entity in publications['entities']:
        cur_end = entity['end']
        cur_start = entity['start']
        cur_score = float(entity['score'])
        cur_entity_group = entity['entity_group']
        cur_word = entity['word']

        should_merge = False

        # Merge if tokens are close enough regardless of entity group
        if prev_word is not None:
            if (prev_end >= cur_start - threshold_gap and prev_end <= cur_start):
                should_merge = True

        if should_merge:
            combine_word = prev_word.rstrip() + " " + cur_word.lstrip()

            merged_num += 1
            prev_score += cur_score
            prev_end = cur_end
            prev_word = combine_word

            # Pick entity_group of higher confidence
            if cur_score > prev_score / merged_num:
                prev_entity_group = cur_entity_group

        else:
            if prev_word is not None:
                new_score = prev_score / merged_num
                cur_entities.append({
                    "entity_group": prev_entity_group,
                    "score": str(new_score),
                    "word": prev_word,
                    "start": prev_start,
                    "end": prev_end
                })

            prev_word = cur_word
            prev_end = cur_end
            prev_start = cur_start
            prev_score = cur_score
            prev_entity_group = cur_entity_group
            merged_num = 1

    if prev_word is not None:
        new_score = prev_score / merged_num
        cur_entities.append({
            "entity_group": prev_entity_group,
            "score": str(new_score),
            "word": prev_word,
            "start": prev_start,
            "end": prev_end
        })

    merged_words_entity_groups.append({
        "pmid": publications['pmid'],
        "title": publications['title'],
        "abstract": publications['abstract'],
        "entities": cur_entities
    })
recompute_offsets(merged_words_entity_groups)

# Output 1 — Aggregated entity vocabulary
entity_aggregate = defaultdict(lambda: {"count": 0, "score_sum": 0.0, "entity_group": None})

for pub in merged_words_entity_groups:
    for e in pub['entities']:
        key = (e['word'], e['entity_group'])
        entity_aggregate[key]['count'] += 1
        entity_aggregate[key]['score_sum'] += float(e['score'])
        entity_aggregate[key]['entity_group'] = e['entity_group']

final_kg_entities = []
for (word, entity_group), data in entity_aggregate.items():
    avg_score = data['score_sum'] / data['count']
    final_kg_entities.append({
        "word": word,
        "entity_group": entity_group,
        "avg_score": avg_score,
        "count": data['count']
    })

filtered_medical_entity = {
    "topic": rawMedicalEntities['topic'],
    "publications": final_kg_entities
}

with open("FilteredMedicals.json", mode="w", encoding="utf-8") as write_file:
    json.dump(filtered_medical_entity, write_file, indent = 2)

# Output 2 — Fully merged entities per publication for relationship extraction
filtered_publication_level = {
    "topic": rawMedicalEntities['topic'],
    "publications": merged_words_entity_groups
}

with open("FilteredPerPublication.json", mode="w", encoding="utf-8") as write_file:
    json.dump(filtered_publication_level, write_file, indent = 2)

print("Done. Smart merged files generated.")
