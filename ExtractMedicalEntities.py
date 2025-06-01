import json

from transformers import pipeline
from transformers import AutoTokenizer, AutoModelForTokenClassification


class EncodeMedicalEntity:
    medical_entities = []
    tokenizer = AutoTokenizer.from_pretrained("d4data/biomedical-ner-all")
    model = AutoModelForTokenClassification.from_pretrained("d4data/biomedical-ner-all")

    pipe = pipeline("ner", model=model, tokenizer=tokenizer, aggregation_strategy="simple") # pass device=0 if using gpu
    results = pipe("""The patient reported no recurrence of palpitations at follow-up 6 months after the ablation.""")
    print(results)

    file_name = 'publications.json'

    with open(file_name, 'r') as file:
        paper_json = json.load(file)
    # print(paper_json)
    store = paper_json['publications']

    for publication in paper_json['publications']:
        medical_entity = {
            "pmid": publication['pmid'],
            "title": publication['title'],
            "abstract": publication['abstract'],
            "entities": [{**item, 'score': str(item['score'])} for item in pipe(publication['abstract'] )] 
            }

        medical_entities.append(medical_entity)

    final_medical_entity = {
        "topic": paper_json['topic'],
        "publications": medical_entities
    }
    with open("medicalEntities.json", mode="w", encoding="utf-8") as write_file:
        json.dump(final_medical_entity, write_file, indent = 2)