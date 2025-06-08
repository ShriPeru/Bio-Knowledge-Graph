from Bio import Entrez
import pprint
import json
Entrez.email = "your.email@example.com"

class ExtractData:
    papers = []
    topic = "sleep improves memory"
    handle = Entrez.esearch(db="pubmed", term=topic, retmax=5)
    record = Entrez.read(handle)
    handle.close()

    id_list = record["IdList"]
    handle = Entrez.efetch(db="pubmed", id=",".join(id_list), rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    for article in records['PubmedArticle']:
        abstract_data = article['MedlineCitation']['Article'].get('Abstract')
        title = article['MedlineCitation']['Article']['ArticleTitle']
        id = article['MedlineCitation'].get('PMID')
        if abstract_data:
            abstract_text = ' '.join(abstract_data['AbstractText'])
            data = {
                "pmid": str(id),
                "title": str(title),
                "abstract": str(abstract_text)

            }
            papers.append(data)
        else:
            print("No abstract available")
    final_data = {
        "topic": topic,
        "publications": papers
    }
    with open("publications.json", mode="w", encoding="utf-8") as write_file:
        json.dump(final_data, write_file, indent = 2)