import requests
import pandas as pd
from Bio.KEGG import REST
import io


# Consultando los organismos disponibles en Kegg
result = REST.kegg_list("organism").read()

# Transformando los resultados a un dataframe
df = pd.read_table(io.StringIO(result), header=None, sep="\t")
print(df)

filtro = df[df[2].str.contains("babesia", case=False)]
print(filtro[[1,2]])

organism_codes = filtro[1].tolist()

# Consultando los genes asociados a B. bovis (bbo)
result = REST.kegg_list("bbo").read()

# Transformando los resultados a un dataframe
df = pd.read_table(io.StringIO(result), header=None, sep="\t")
print(df)

genes_bbo = df.iloc[:,0].tolist()
print(len(genes_bbo))

print(REST.kegg_get("dme:Dmel_CG17636").read())

# Rutas metabolicas (pathways), presentes en B. bovis
result = REST.kegg_list("pathway", "bbo").read()

# Transformando los resultados a un dataframe
df = pd.read_table(io.StringIO(result), header=None)
print(df)

def pathways_for_gene(gene_id):
    url = f"https://rest.kegg.jp/link/pathway/{gene_id}"
    r = requests.get(url)

    if r.status_code != 200:
        print("HTTP error:", r.status_code)
        return None

    lines = r.text.strip().split("\n")
    if not lines or lines == ['']:
        print(f"No pathways found for {gene_id}")
        return []

    pathways = [line.split("\t")[1] for line in lines]
    return pathways

genes = [
    # B. bovis
    "bbo:BBOV_I001300",
    "bbo:BBOV_I001310",
    "bbo:BBOV_I001320",
    "bbo:BBOV_I001330",
    "bbo:BBOV_I001340",
    "bbo:BBOV_IV011130",
    # Drosofila
    "dme:Dmel_CG17636",
    "dme:Dmel_CG40494",
    "dme:Dmel_CR43552",
    "dme:Dmel_CG17131"]

n = 0 
i = 0
genes_withp = [] # Empty
while n < 15:
    print(genes_bbo[i], paths := pathways_for_gene(genes_bbo[i]))
    i += 1
    if len(paths) > 0:
        genes_withp.append(paths)
        n += 1
    print(n)