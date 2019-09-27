from hmmerMethods import parse_hmmer_output
import pandas as pd
from getMethods import alignment_list, tfc_uri_to_label, get_children_nodes
from database import PROTEOME_PATH, FUSEKI_URL, FCREPO_URL
import csv
from queryCreator import createquery
import requests as rq
from parseMethods import parseturtle
import sqlite3


class Statichits:
    gen_hits = 0
    newgen_hits = 0
    grp_hits = 0
    fc_hits = 0
    maybe_hits = 0
    hits = 0
    ng_hits = 0


# Checks number of genera connected to a subfamily
def check_genera(tfc_id):
    nr_genera = 0
    al_list = alignment_list(5)
    for al in al_list:
        check_this = tfc_uri_to_label(al)
        if match_id(tfc_id, check_this) == -1:
            nr_genera += 1
    return nr_genera


# Checks if two TFC IDs belong to the same branch
# Returns the point where the two IDs diverge, or -1 if they dont diverge at all
def match_id(id1, id2):
    id1s = id1.split('.')
    id2s = id2.split('.')
    to_here = min(len(id1s), len(id2s))
    for i in range(0, to_here):
        if id1s[i] != id2s[i]:
            return i
    return -1


# Returns the classification level of a given TFCID
def level(tfc_id):
    tfc_id = tfc_id.split('.')
    if len(tfc_id) == 1:
        return "Superclass"
    elif len(tfc_id) == 2:
        return "Class"
    elif len(tfc_id) == 3:
        return "Family"
    elif len(tfc_id) == 4:
        return "Subfamily"
    else:
        return "Genus"


# Builds "hit blocks" (Pandas dataframes) out of a hit list divided by subfamilies
def build_blocks(hitlist):
    blocks = []
    block_id = 0
    classes = []

    # First, identify subfamilies
    for line in hitlist:
        if level(line[1]) == "Subfamily":
            if line[1] not in classes:
                classes.append(line[1])
                block_id += 1
        elif level(line[1]) == "Genus":
            splitter = line[1].split('.')
            sub = splitter[0] + "." + splitter[1] + "." + splitter[2] + "." + splitter[3]
            if sub not in classes:
                classes.append(sub)
                block_id += 1

    # Match all hits to respective subfamilies
    for line in hitlist:
        for i in range(0, block_id):
            if match_id(line[1], classes[i]) == -1:
                blocks.append([i, line[1], line[2], float(line[3]), float(line[4]), line[5], float(line[6]), float(line[7])])

    block_df = pd.DataFrame(columns=['BlockID', 'HitDescriptor', 'E-Value', 'Score', 'Bias', 'E-Value Best Domain',
                                     'Score Best Domain', 'Bias Best Domain'], data=blocks)
    block_df = block_df.infer_objects()
    block_df_score = block_df.sort_values(by='Score', ascending=False)
    block_df_sorted = block_df.sort_values(by=['BlockID', 'Score'], ascending=[True, False])
    return block_df_score, block_df_sorted


# Evaluates a given list ("block") with the matches that a protein got with the hmmer profile(s)
def evaluate_protein_hits(block_dataframe_score, block_dataframe_sorted, prot_id, resultfile, dist_fct):
    if len(block_dataframe_score) == 0:
        return None

    # Gets lots of statistical values for the block
    def evaluate_block(block, pid, out):
        if len(block) == 0:
            return None
        mean_dist_block = mean_score_block = mean_score_dom = best_hit_score = best_hit_dscore = best_hit_distance = 0.0
        distances = []
        ranking = []
        best_hit = ""
        f_counter = 0
        f_current_score = -1
        max_dist = -1.0
        sfm_rank = fam_rank = cls_rank = -1

        for f_index, f_row in block.iterrows():
            ranking.append([f_row['HitDescriptor'], level(f_row['HitDescriptor'])])
            if level(f_row['HitDescriptor']) == "Subfamily":
                sfm_rank = f_counter
            if level(f_row['HitDescriptor']) == "Family":
                fam_rank = f_counter
            if level(f_row['HitDescriptor']) == "Class":
                cls_rank = f_counter

            if f_counter == 0:
                best_hit = f_row['HitDescriptor']
                best_hit_score = f_row['Score']
                best_hit_dscore = f_row['Score Best Domain']

            mean_score_block += f_row['Score']
            mean_score_dom += f_row['Score Best Domain']

            if f_counter != 0:
                mean_dist_block += float((f_row['Score'] / f_current_score))
                distances.append(float(f_row['Score'] / f_current_score))
                if f_counter == 1:
                    best_hit_distance = round(mean_dist_block, 2)
                max_dist = round(max(max_dist, f_row['Score'] / f_current_score), 2)
            if f_counter != len(block)-1:
                f_current_score = f_row['Score']

            f_counter += 1

        if mean_dist_block != 0 and len(block)-1 != 0:
            mean_dist_block = round(mean_dist_block/(len(block)-1), 2)
        mean_score_block = round(mean_score_block/len(block), 2)
        mean_score_dom = round(mean_score_dom/len(block), 2)

        # Output to file if output variable is set
        if pid < 5 and out:
            resultfile.write("-------Evaluating a new block-------\n")
            resultfile.write(block.to_string())
            resultfile.write("\nBlockID:\t\t" + str(pid) + " with " + str(len(block)) + " elements.\n")
            resultfile.write("Best hit:\t\t" + best_hit + " with score " + str(best_hit_score) + " (Sequence), " +
                             str(best_hit_dscore) + " (Domain) and distance to second best hit " +
                             str(best_hit_distance) + "\n")
            resultfile.write("Mean score:\t\t" + str(mean_score_block) + " (FullSequence), " + str(mean_score_dom) +
                             " (Domain) and mean distance " + str(mean_dist_block) + "\n")
            resultfile.write("Subfamily at " + str(sfm_rank) + ", family at " + str(fam_rank) + " and class at " +
                             str(cls_rank) + "\n")
            resultfile.write("Max distance: " + str(max_dist) + "\n")

        # Categorizes the query according to the given block
        resultfile.write("-------------Categorization here--------------\n")
        significants = []
        for dist in distances:
            if dist == 0.0:
                significants.append(2)
            elif float(dist) < dist_fct * mean_dist_block:
                significants.append(1)
            else:
                significants.append(0)
        significants.append(-1)

        sig_group = [ranking[0]]
        rk = 1
        f2 = False
        for sig in significants:
            if sig == 0:
                sig_group.append(ranking[rk])
                rk += 1
            elif sig == 2:
                sig_group.append(ranking[rk])
                rk += 1
                f2 = True
            else:
                cat = wtf(sig_group, resultfile, f2)
                return cat, best_hit
        cat = wtf(sig_group, resultfile, f2)
        return cat, best_hit

    resultfile.write("--------------------------------------------------------------------\n")
    resultfile.write("Evaluating a new superblock - Query: " + prot_id + "\n")
    # print(block_dataframe_score)
    # Calculate mean distance for the entire protein sequence / the best domain aswell as mean scores
    current_score = current_domain_score = -1
    ms = ms10 = msb = msb10 = 0.0
    md = mdb = md10 = mdb10 = 0.0
    counter = 0
    for index, row in block_dataframe_score.iterrows():
        ms += row['Score']
        if counter != 0:
            md += (row['Score']/current_score)
        if counter != len(block_dataframe_score)-1:
            current_score = row['Score']
        counter += 1
        if counter == 10 or (counter < 10 and counter == len(block_dataframe_score)-1):
            md10 = round(md/counter, 2)
            ms10 = round(ms/counter, 2)
    if len(block_dataframe_score)-1 != 0:
        md = round(md/(len(block_dataframe_score)-1), 2)
    ms = round(ms/len(block_dataframe_score), 2)
    block_dataframe_score = block_dataframe_score.sort_values(by='Score Best Domain', ascending=False)
    counter = 0
    for index, row in block_dataframe_score.iterrows():
        msb += row['Score Best Domain']
        if counter != 0:
            mdb += (current_domain_score - row['Score Best Domain'])
        if counter != len(block_dataframe_score)-1:
            current_domain_score = row['Score Best Domain']
        counter += 1
        if counter == 10 or (counter == len(block_dataframe_score)-1 and counter < 10):
            mdb10 = round(mdb/counter, 2)
            msb10 = round(msb/counter, 2)
    if len(block_dataframe_score) - 1 != 0:
        mdb = round(mdb/(len(block_dataframe_score)-1), 2)
    msb = round(msb/len(block_dataframe_score), 2)
    resultfile.write("Mean distance FullSeq: " + str(md) + ". Mean distance best domain: " + str(mdb) + "\n")
    resultfile.write("Mean distance TopSeq: " + str(md10) + ". Mean distance Top best domains: " + str(mdb10) + "\n")
    resultfile.write("Mean score FullSeq: " + str(ms) + ". Mean score best domain: " + str(msb) + "\n")
    resultfile.write("Mean score TopSeq: " + str(ms10) + ". Mean score Top best domains: " + str(msb10) + "\n")

    # Secondly, evaluate the divided blocks on their own
    best = "None"
    category = 0
    current_block_id = 0
    cutoff = 0
    counter = -1
    cat_list = []
    for index, row in block_dataframe_sorted.iterrows():
        counter += 1  # Current row
        if current_block_id == 0:
            if row['BlockID'] != current_block_id:
                category, best = evaluate_block(block_dataframe_sorted[cutoff:counter], current_block_id, True)
                cat_list.append(category)
                cutoff = counter
                current_block_id = row['BlockID']
            elif counter == len(block_dataframe_sorted)-1:
                category, best = evaluate_block(block_dataframe_sorted[cutoff:counter+1], current_block_id, True)
                cat_list.append(category)
    if not cat_list:
        cat_list.append(-1)
    return cat_list, best


# Divides HMMER Tbl file into parsable hitlists
def divide(tbl, thres, dist_f):

    Statichits.gen_hits = 0
    Statichits.grp_hits = 0
    Statichits.newgen_hits = 0
    Statichits.fc_hits = 0
    Statichits.maybe_hits = 0
    Statichits.hits = 0
    Statichits.ng_hits = 0

    with open('/sybig/home/ttu/BSc/category_out.txt', 'w') as result_out:
        returnlist = []
        queries = 0
        hmmer_list = parse_hmmer_output(tbl)
        current_prot = ''
        cutoff = 0
        index = -1
        for line in hmmer_list:
            index += 1
            if current_prot == '':
                current_prot = line[0]
                queries = 1
            if line[0] != current_prot or index == len(hmmer_list)-1:
                # uniprot_symbol = current_prot.split('|')[1]
                # gene_symbol = current_prot.split('|')[2]
                # gene_symbol = gene_symbol.split('_')[0]
                uniprot_symbol = "None"
                gene_symbol = "None"
                # gene_symbol = current_prot.split("_")
                # gene_symbol = gene_symbol[len(gene_symbol)-1].split("---")[0]
                queries += 1
                best_til_here = hmmer_list[cutoff:index][0][1]
                protein_query = clean_queries(hmmer_list[cutoff:index], thres)
                if not protein_query:
                    result_out.write("--------------------------------------------------------------------\n")
                    result_out.write("No fitting match in database found for " + str(current_prot) + ".\n")
                    returnlist.append([uniprot_symbol, gene_symbol, -1, best_til_here])
                    Statichits.maybe_hits += 1
                else:
                    newblock_score, newblock_sorted = build_blocks(protein_query)
                    if not newblock_score.empty and not newblock_sorted.empty:
                        category, best = evaluate_protein_hits(newblock_score, newblock_sorted, current_prot, result_out, dist_f)
                        if category:
                            category = category[0]
                        else:
                            category = 666
                        returnlist.append([uniprot_symbol, gene_symbol, category, best, current_prot.split("---")[1]])

                cutoff = index
                current_prot = line[0]
        print(str(queries) + " queries meet HMMer requirements")
        print(str(Statichits.hits) + " hits and " + str(Statichits.maybe_hits) + " theoretical extra hits (E-Value"
                                                                                 " not significant)")
        print(str(Statichits.gen_hits) + " direct Genus hits and " + str(Statichits.newgen_hits) + " hits that cant be"
                                                                                                   "categorized (might"
                                                                                                   "be a new Genus)")
        print(str(Statichits.ng_hits) + " hit a subfamily with no genera.")
        print(str(Statichits.fc_hits) + " hits that matched families and classes first")
        print(str(Statichits.grp_hits) + " hits match a group of similar nodes and cant be directly categorized")
    return returnlist


# Throws out queries that dont meet the e-value threshhold
def clean_queries(query_list, threshold):
    cleaned_list = []
    for line in query_list:
        if "e" in line[2] and "e" in line[5]:
            ev1s = line[2].split('e')
            ev2s = line[5].split('e')
            if -int(ev1s[1].lstrip("-").lstrip("0")) <= threshold and -int(ev2s[1].lstrip("-").lstrip("0")) <= threshold:
                cleaned_list.append(line)
    return cleaned_list


# Writes categorization results to file
def wtf(result_list, result_out, f2):
    category = 0
    if len(result_list) == 1:
        if result_list[0][1] == "Genus":
            Statichits.gen_hits += 1
            Statichits.hits += 1
            result_out.write("Sequence is quite clearly matched to existing Genus " + result_list[0][0] + "\n")
            category = 1
        elif result_list[0][1] == "Subfamily":
            Statichits.hits += 1
            if check_genera(result_list[0][0]) == 0:
                result_out.write("Sequence matches a subfamily with significant distance, but there are no genera\n")
                Statichits.ng_hits += 1
                category = 0
            else:
                result_out.write("Sequence matches a subfamily (" + result_list[0][0] + ") first - "
                                                                                        "Might be a new Genus\n")
                Statichits.newgen_hits += 1
                category = 2
        elif result_list[0][1] == "Family" or result_list[0][1] == "Class":
            Statichits.fc_hits += 1
            Statichits.hits += 1
            result_out.write("Sequence matches family or class node - Impossible to categorize\n")
            category = 3
    else:
        Statichits.hits += 1
        if f2:
            result_out.write("Sequence matches a group of nodes, among them distances of 0:\n")
            category = 4
            if result_list[0][1] == "Genus":
                result_out.write("Since first match is a genus, it is attributed to that genus node.\n")
                category = 1
                Statichits.gen_hits += 1
            else:
                Statichits.grp_hits += 1
        else:
            category = 5
            result_out.write("Sequence matches a group of nodes with fairly similar scores:\n")
            if result_list[0][1] == "Genus":
                Statichits.gen_hits += 1
                result_out.write("Since first match is a genus, it is attributed to that genus node.\n")
                category = 1
            else:
                Statichits.grp_hits += 1
        for res in result_list:
            result_out.write(str(res[0])+"  \t\t"+str(res[1])+"\n")
    return category


# Gathers UNIPROT IDs from the temporary sequence file, and gene symbols
def investigate_tsf():
    uniprot_list, gene_symbol_list = [], []
    with open('/sybig/home/ttu/BSc/Sequences/temp_seqfile', "r") as tfs:
        tfs_text = tfs.readlines()
        for line in tfs_text:
            if line[0] == '>':
                line_split = line.split('_')
                line_split2 = line_split[len(line_split)-1].split("---")
                uniprot_id = "-"
                gene_symbol = line_split2[0]
                uniprot_list.append(uniprot_id)
                gene_symbol_list.append(gene_symbol)
    return uniprot_list, gene_symbol_list


# Gathers UNIPROT IDs out of the proteome so we can later cross-check which of those made it through HMMer
# Also gathers "gene symbols" like ALX1
def investigate_proteom():
    uniprot_list, gene_symbol_list = [], []
    with open(PROTEOME_PATH, "r") as proteome:
        proteome_text = proteome.readlines()
        for line in proteome_text:
            if line[0] == '>':
                line_split = line.split('|')
                uniprot_id = line_split[1]
                gene_symbol = line_split[2].split('_')
                gene_symbol = gene_symbol[0]
                uniprot_list.append(uniprot_id)
                gene_symbol_list.append(gene_symbol)
    return uniprot_list, gene_symbol_list


# Gathers a list of all human transcription factors, their uniprot IDs and gene symbols (if available)
# Uses the Fedora repository
def repo_gather_human():
    out_list = []
    ch_list = get_children_nodes(FCREPO_URL + 'TF_protein/tf_9606')
    for child in ch_list:
        add1, add2 = "-", "-"
        out = rq.get(child.rstrip(".")).text
        parsed = parseturtle(out, ['rdfs:label', 'sybig:xref'], 0)
        if not parsed:
            print("No turtlefile results found")
        else:
            for element in parsed:
                if "UNIPROT" in element:
                    add1 = element.split(':')[len(element.split(':'))-1]
            add2 = parsed[len(parsed)-1]
            out_list.append([child, add1, add2])
    # print(out_list)
    return out_list


# Compares the categorization results with TFClass to find true and false negatives/positives
# Output: Uniprot / GS / Category / Best Hit //// TFClass ID / Has DBD? / Has Profile?
def investigate_tfclass():
    profile_list = []
    conn = sqlite3.connect("Sequences.db")
    cu = conn.cursor()
    tfclass_list = []
    found_list_uniprot, found_list_gs = [], []
    ctr = 0
    tfclass_checklist = repo_gather_human()
    print(str(len(tfclass_checklist)) + " human TFs found in the TFClass repository")
    al_ch_list_full, al_ch_list = alignment_list(5), []
    for al in al_ch_list_full:
        alx = al.split('/')[len(al.split('/'))-1]
        al_ch_list.append(alx)
    with open('/sybig/home/ttu/BSc/test.csv', "r") as proteome_out:
        log = open('/sybig/home/ttu/BSc/testlog.csv', 'w')
        for strng in tfclass_checklist:
            log.write(str(strng[0]) + " " + str(strng[1]) + " " + str(strng[2]) + "\n")
        reader = csv.reader(proteome_out, delimiter="\t")
        for row in reader:
            conn_profile, best_hit, conn_dbd, conn_gs, conn_tfc, conn_p = "No Profile", "None", "No DBD", "-", "-", "-"
            last_gs, ent = "-", "-"
            query_uniprot = row[0]
            query_gs = row[1]
            for ch in tfclass_checklist:
                if ch[1] == query_uniprot or ch[2] == query_gs or ch[2].replace('ZNF', 'ZN') == query_gs or\
                        ch[2].replace('ZNF', 'Z') == query_gs:  # Specifically look at proteins that are in TFClass
                    ctr += 1
                    out = rq.get(ch[0]).text
                    conn_p = tfc_uri_to_label(ch[0])
                    conn_tfc = parseturtle(out, ['sybig:belongs_to'], 0)
                    if conn_tfc:
                        conn_tfc = tfc_uri_to_label(conn_tfc[0])
                    dbdl = parseturtle(out, ['sybig:dbd'], 0)
                    if dbdl:
                        conn_dbd = "Has DBD"
                    if conn_gs:
                        conn_gs = parseturtle(out, ['rdfs:label'], 0)
                        if conn_gs:
                            conn_gs = conn_gs[0]
                    if conn_tfc in al_ch_list:
                        profile_list.append(conn_tfc)
                        conn_profile = "Has Profile"
                        for res in cu.execute("SELECT Entropy FROM ALIGNMENTS WHERE TfcID = \"" + str(conn_tfc) + "\""):
                            ent = res[0]
                    log.write(
                        "Adding " + query_uniprot + " with query GS " + query_gs + " and category " + row[2] + ".\n")
                    log.write("Best hit HMMer descriptor is " + row[3] + ".\n")
                    log.write("Connected TFClass protein is " + conn_p + ".\n")
                    log.write(conn_dbd + ", " + conn_profile + ".\n")
                    last_gs = conn_gs
                    tfclass_list.append([query_uniprot, query_gs, row[2], row[3], conn_p, conn_tfc, conn_gs, conn_dbd,
                                         conn_profile, ent])
                    if ch[1] == query_uniprot:
                        found_list_uniprot.append(query_uniprot)
                    if ch[2] == query_gs:
                        found_list_gs.append(query_gs)
            if last_gs == "-":
                tfclass_list.append([query_uniprot, query_gs, row[2], row[3], conn_p, conn_tfc, conn_gs, conn_dbd,
                                    conn_profile, ent])

    log.close()
    with open('/sybig/home/ttu/BSc/test2ort.csv', "w") as csv_in:
        test_writer = csv.writer(csv_in, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL)
        test_writer.writerow(["Query", "Gene Symbol", "Category", "HMMer HitDescriptor", "TFClass Node",
                              "TFClass Node ID", "TFClass Node GS", "TFClass node has DBD?", "In Profiles", "Entropy"])
        for e in tfclass_list:
            test_writer.writerow([e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9]])
    print(str(ctr) + " human transcription factors with either UNIPROT ID or fitting gene symbol detected.")
    print(str(len(profile_list)) + " profile hits")
    print("Transcription factors that have not been found, but should be according to profiles:")
    for tcl in tfclass_checklist:
        if tcl[1] not in found_list_uniprot and tcl[2] not in found_list_gs:
            cand_out = rq.get(tcl[0].rstrip(".")).text
            ctl = parseturtle(cand_out, ['sybig:belongs_to'], 0)
            if ctl:
                if tfc_uri_to_label(ctl[0]) in al_ch_list:
                    # print(tfc_uri_to_label(ctl[0]))
                    for res in cu.execute("SELECT Entropy FROM ALIGNMENTS WHERE TfcID =\"" +
                                          str(tfc_uri_to_label(ctl[0]))+"\""):
                        ent = res[0]
                    print(tcl[0] + ", " + tcl[1] + ", " + tcl[2] + " " + str(ent))
    conn.close()
    return tfclass_list


# Evaluate the resulting CSV according to whether or not the hits have profiles etc.
def evaluate_results(file):
    row_counter = -1
    with open('/sybig/home/ttu/BSc/' + file, 'r') as results:
        total = 0
        true_pos, fp, false_pos, tnp = 0, 0, [0, 0, 0, 0, 0], 0
        bullshit, bullshit_c45, similar_c45 = 0, 0, 0
        sfg_matches, sfg_matches_c45 = 0, 0
        no1cares = 0
        should_find_this, this_too = 0, 0
        subfams_matched, subfams_wrong = 0, 0
        res_reader = csv.reader(results, delimiter="\t")
        for row in res_reader:
            row_counter += 1
            if row_counter == 0:
                continue
            total += 1
            category = int(row[2])
            res_hmmer = row[3]
            if "_" in res_hmmer:
                res_hmmer = res_hmmer.split("_")[0]
            res_tfc = row[5]
            prof = False
            if row[8] == "Has Profile":
                prof = True

            m_ret = match_id(res_tfc, res_hmmer)
            if category == 1:
                if prof:
                    if res_hmmer == res_tfc:
                        # Correct genus matched
                        true_pos += 1
                    else:
                        false_pos[m_ret] += 1
                        if m_ret < 4:
                            print(str(m_ret) + " - " + res_hmmer + " / " + res_tfc)
                else:

                    if m_ret == -1:
                        # Found a genus completely
                        tnp += 1
                    elif m_ret < 3:
                        # Found a genus
                        sfg_matches += 1
                    else:
                        bullshit += 1

            elif category == 2:
                if prof:
                    if m_ret == -1:
                        subfams_matched += 1
                    else:
                        subfams_wrong += 1
                else:
                    no1cares += 1

            elif category == 0:
                no1cares += 1

            elif category == -1:
                if prof:
                    should_find_this += 1
                else:
                    no1cares += 1

            elif category == 3:
                if prof:
                    this_too += 1
                else:
                    no1cares += 1

            elif category == 4 or category == 5:
                if prof:
                    if m_ret == -1:
                        similar_c45 += 1
                    elif m_ret < 3:
                        bullshit_c45 += 1
                    else:
                        sfg_matches_c45 += 1
                else:
                    no1cares += 1

    print(str(total) + "queries in total")
    print(str(no1cares) + " queries were not interesting enough to actually consider.")
    print(str(should_find_this) + " queries should have been considered, but weren't.\n")
    print(str(true_pos) + " genera were correctly identified by HMMer and matched correctly.")
    print(str(fp) + " genera were identified by HMMer, but matched to a wrong genus.")
    print(str(false_pos[0]) + ", " + str(false_pos[1]) + ", " + str(false_pos[2]) + ", " + str(false_pos[3]) + ", " +
                                                                                    str(false_pos[4]))
    print(str(tnp) + " query sequences were not found in TFClass, but still matched to a perfect genus")
    print(str(sfg_matches) + " query sequences were not found in TFClass, but still matched to a genus and matched to "
                             "a fitting subfamily.")
    print(str(bullshit) + " query sequences were not found in TFClass, but still matched to a genus and matched to a "
                          "wrong subfamily.\n ")
    print(str(subfams_matched) + " query sequences were matched to the correct subfamily despite having a "
                                 "genus profile.")
    print(str(subfams_wrong) + " query sequences were matched to a wrong subfamily despite having a genus profile.\n")
    print(str(this_too) + " query sequences were matched to a family profile despite having a genus profile.\n")
    print(str(similar_c45) + " query sequences were matched to groups, with the correct subfamily as best hit, "
                             "despite a genus profile existing.")
    print(str(sfg_matches_c45) + " query sequences were matched to groups, with a dubious best hit, despite a "
                                 "genus profile existing.")
    print(str(bullshit_c45) + " query sequences were matched to groups, with a wrong best hit, despite a genus "
                              "profile existing.")

    return (no1cares, should_find_this, true_pos, false_pos, sfg_matches, bullshit,
            subfams_matched, subfams_wrong, this_too, similar_c45, sfg_matches_c45, bullshit_c45)


def median_entropy(file):
    row_counter = -1
    with open('/sybig/home/ttu/BSc/' + file, 'r') as results:
        cat1, cat2, cat3, cat4, cat5 = 0.0, 0.0, 0.0, 0.0, 0.0
        cat1c, cat2c, cat3c, cat4c, cat5c = 0, 0, 0, 0, 0
        reader = csv.reader(results, delimiter="\t")
        for row in reader:
            row_counter += 1
            if row_counter == 0:
                continue
            if row[2] == "-1":
                cat5c += 1
                cat5 += float(row[9])
            elif row[2] == "3" or row[2] == "0":
                cat4c += 1
                cat4 += float(row[9])
            elif row[2] == "1":
                cat1c += 1
                cat1 += float(row[9])
            elif row[2] == "2":
                cat3c += 1
                cat3 += float(row[9])
            else:
                cat2c += 1
                cat2 += float(row[9])
        cat1 /= cat1c
        cat2 /= cat2c
        cat3 /= cat3c
        cat4 /= cat4c
        cat5 /= cat5c
    return cat1, cat2, cat3, cat4, cat5


def investigate_orts():
    profile_list = []
    conn = sqlite3.connect("Sequences.db")
    cu = conn.cursor()
    tfclass_list = []
    al_ch_list_full, al_ch_list = alignment_list(5), []
    for al in al_ch_list_full:
        alx = al.split('/')[len(al.split('/')) - 1]
        al_ch_list.append(alx)
    with open('/sybig/home/ttu/BSc/test.csv', "r") as proteome_out:
        log = open('/sybig/home/ttu/BSc/testlog.csv', 'w')
        reader = csv.reader(proteome_out, delimiter="\t")
        for row in reader:
            conn_profile, best_hit, conn_dbd, conn_gs, conn_p = "No Profile", "None", "No DBD", "-", "-"
            ent = "-"
            conn_tfc = row[4]
            query_uniprot = row[0]
            query_gs = row[1]
            for res in cu.execute("SELECT Entropy FROM ALIGNMENTS WHERE TfcID = \"" + str(conn_tfc) + "\""):
                    ent = res[0]
            tfclass_list.append([query_uniprot, query_gs, row[2], row[3], conn_p, conn_tfc, conn_gs, conn_dbd,
                                conn_profile, ent])
    with open('/sybig/home/ttu/BSc/test2ort.csv', "w") as csv_in:
        test_writer = csv.writer(csv_in, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL)
        test_writer.writerow(["Query", "Gene Symbol", "Category", "HMMer HitDescriptor", "TFClass Node",
                              "TFClass Node ID", "TFClass Node GS", "TFClass node has DBD?", "In Profiles", "Entropy"])
        for e in tfclass_list:
            test_writer.writerow([e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9]])
