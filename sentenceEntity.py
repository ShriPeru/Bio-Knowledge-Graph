import scispacy
import spacy
import json

class mapSentences:
    def __init__(self, start, end, sentence, entity_list = []):
        print('sentence: ', sentence)
        print('start ', start)
        print('end ', end)
        self.start = start
        self.end = end
        self.sentence = sentence
        self.entity_list = entity_list

    def placeEntity(self, entity):
        if (entity['start'] >= self.start) and (entity['end'] <= self.end):
            self.entity_list.append(entity)
            return True
        else:
            return False
    def convertToJson(self):
        sentence_obj = {
            "sentence": self.sentence,
            "start": self.start,
            "end": self.end,
            "entities": self.entity_list
        }
        return sentence_obj
    

class SentenceFormation:
    nlp = spacy.load("en_core_sci_sm")
    # abstract = "Declining spatial working memory (WM) is an early hallmark of Alzheimer's disease (AD). Sleep disturbance exacerbates spatial WM and increases AD risk. The GABAergic system, crucial for sleep regulation, may mediate this link. We thus investigate the relationship between spatial WM and hippocampal GABAergic signaling during rapid eye movement sleep deprivation (REM-SD) in AD model mice. We assessed spatial and non-spatial WM, locomotor activity, and anxiety-like behavior in 6-month-old triple transgenic (3xTg) AD mice and wild-type (WT) controls, with and without REM-SD (5 days, 4 h/day). We then used immunofluorescence to quantify GABA<sub>A</sub>\u03b11, GABA<sub>B</sub>R1, GAD67, and GABA levels in the prefrontal cortex (PFC) and hippocampus and analyze the correlations with behavioral outcomes. REM-SD increased locomotor activity, reduced anxiety-like behavior, and improved non-spatial WM in 3xTg-AD mice. Conversely, REM-SD impaired spatial WM in WT mice, which was also demonstrated in 3xTg-AD mice. Increased hippocampal GABA levels are correlated with improved non-spatial WM in 3xTg+SD mice. In contrast, impaired spatial WM in WT+SD mice was associated with elevated hippocampal GABA and GABA<sub>B</sub>R1, decreased hippocampal GAD67, and reduced PFC GABA levels. Notably, spatial WM in 3xTg+SD and 3xTg control mice related to increased GABA<sub>A</sub>\u03b11 in the PFC and hippocampus and GAD67 in hippocampal CA1, along with decreased GABA<sub>B</sub>R1 and GAD67 in the dentate gyrus. REM-SD-induced alterations in WM performance are linked to GABAergic signaling changes in the PFC and hippocampus, with distinct patterns in WT and 3xTg-AD mice. This study provides insight into AD pathologies and potential therapeutic targets for sleep-related cognitive impairments."
    # doc = nlp(abstract)
    # for sent in doc.sents:
    #     print(sent.start_char, sent.end_char, sent.text)


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
            it = iter(cur_sent_map)
            map_sentence_obj = next(it)
            try:
                while (map_sentence_obj.placeEntity(entity) == False):
                    map_sentence_obj = next(it)
            except StopIteration:
                print("couldn't place in sentence object")
                print(entity['word'])
        sentences_list = []
        for mappings in cur_sent_map:
            sentences_list.append(mappings.convertToJson())

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
        json.dump(publications_assigned, write_file, indent = 2)

            




    