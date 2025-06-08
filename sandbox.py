import pandas as pd

with open('KG_filtered.json', encoding='utf-8') as inputfile:
    df = pd.read_json(inputfile)

df.to_csv('KG_converted.csv', encoding='utf-8', index=False)