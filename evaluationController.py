from evaluationMethods import investigate_proteom, investigate_tfclass, divide, evaluate_results, investigate_orts
from database import RESULT_PATH
import csv


# Checks one of the lists supplied by "divide" if a uniprotID or genesymbol is contained
def checklist(uniprot, genesymbol, list):
    for ck in list:
        if ck[0] == uniprot or ck[1] == genesymbol:
            return True
    return False


res_w = open("results_prot_hs2.txt", 'w')
        # proteome_uniprot_list, proteome_gene_symbol_list = investigate_proteom()
        # Parse the HMMER output
        # pl = divide('/sybig/home/ttu/BSc/ScanResults/26_08_2019 13_52_57/Tbl', -2, 1.0)  # Homo sapiens
        # pl = divide('/sybig/home/ttu/BSc/ScanResults/02_09_2019 16_38_33/Tbl', -2, 1.0)  # Droso
        # pl = divide('/sybig/home/ttu/BSc/ScanResults/14_09_2019 15_27_17/Tbl', -2, actualj)  # native TF seqs
pl = divide('/sybig/home/ttu/BSc/ScanResults/12_08_2019 14_12_40/Tbl', -2, 1.0)
with open('/sybig/home/ttu/BSc/test.csv', "w") as test:
    test_writer = csv.writer(test, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for proteins in pl:
        if len(proteins) == 5:
            test_writer.writerow([proteins[0], proteins[1], proteins[2], proteins[3], proteins[4]])
# investigate_tfclass()
investigate_orts()
retl = evaluate_results("test2ort.csv")
for el in retl:
    res_w.write(str(el) + ", ")
    res_w.write("\n")
res_w.close()


# Sort the HMMER output by queries for each protein, get rid of queries with too low e-values (higher than TH)
# Sort the queries by subfamilies, take the three to five best ones and put them into "blocks"
# Check the mean distance for each block, then check for significant distances inside the block (if blocks are big
#   enough). Consider checking the mean distance over the top X matches aswell.
# Categorize the query into 1: Sequence clearly matches Genus, 2: Sequence quite clearly a new genus
#   3: Distances are so low that it cant rightly be categorized, but it matches to a subfamily/geni block
#   4: Cant be categorized or family/class is the highest match
# Check how good categorization is when dealing with sample sequences
