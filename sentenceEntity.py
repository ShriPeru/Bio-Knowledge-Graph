import scispacy
import spacy
import json

class mapSentences:
    def __init__(self, start, end, sentence):
        self.start = start
        self.end = end
        self.sentence = sentence
        self.entity_list = []

    def placeEntity(self, entity, tolerance=5):
        if (entity['start'] >= self.start - tolerance) and (entity['end'] <= self.end + tolerance):
            self.entity_list.append(entity)

    def convertToJson(self):
        return {
            "sentence": self.sentence,
            "start": self.start,
            "end": self.end,
            "entities": self.entity_list
        }

class SentenceFormation:
    nlp = spacy.load("en_core_sci_sm")

    const_file_name = "FilteredPerPublication.json"

    with open(const_file_name, 'r') as file:
        publications_file = json.load(file)

    publications = publications_file['publications']
    publications_with_assigned_entities = []

    for publication in publications:
        doc = nlp(publication['abstract'])
        cur_sent_map = []

        for sent in doc.sents:
            cur_sent_map.append(mapSentences(sent.start_char, sent.end_char, sent.text))

        for entity in publication['entities']:
            assigned = False
            for sent_map in cur_sent_map:
                sent_map.placeEntity(entity)
                assigned = True  # Whether placed or not is inside placeEntity()
            
            # If entity not placed anywhere, you can optionally log it
            # if not assigned:
            #     print("Unplaced entity:", entity['word'])

        sentences_list = [mapping.convertToJson() for mapping in cur_sent_map]

        publication_sentence = {
            "pmid": publication["pmid"],
            "abstract": publication['abstract'],
            "sentence": sentences_list
        }

        publications_with_assigned_entities.append(publication_sentence)

    publications_assigned = {
        "topic": publications_file["topic"],
        "publications": publications_with_assigned_entities
    }

    with open("AssignedEntities.json", mode="w", encoding="utf-8") as write_file:
        json.dump(publications_assigned, write_file, indent=2)
