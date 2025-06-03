from Bio import Entrez
import pprint
import json
Entrez.email = "your.email@example.com"

class ExtractData:
# Search for papers (example: keyword "cancer")
    papers = []
    # topic = "Mitochondrial dysfunction in neurodegenerative diseases"
    topic = "magnesium"
    handle = Entrez.esearch(db="pubmed", term=topic, retmax=100)
    record = Entrez.read(handle)
    handle.close()

    id_list = record["IdList"]

    # Fetch full records in XML
    handle = Entrez.efetch(db="pubmed", id=",".join(id_list), rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    # print(records)
    # pprint.pprint(records)
    handle.close()
    # print("record title: \n")
    
    # print(records['PubmedArticle']['MedlineCitation']['Article'].get('ArticleTitle'))
    # for article in records['PubmedArticle']:
    #     abstract_title = article['MedlineCitation']['Article'].get('ArticleTitle')
    # for article in records['PubmedArticle']:
    #     title = article['MedlineCitation']['Article']['ArticleTitle']
    #     print(f"Title: {title}\n")
    # # Extract only the abstracts
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
            # print(title)
            # print()
            # print(pmid)
            # print()
            # print(abstract_text)
            # print()
        else:
            print("No abstract available")

    # papers is a list of all json objects 
    # papers_json = json.dumps(papers)
    # print(papers_json)

    final_data = {
        "topic": topic,
        "publications": papers
    }
    with open("publications.json", mode="w", encoding="utf-8") as write_file:
        json.dump(final_data, write_file, indent = 2)