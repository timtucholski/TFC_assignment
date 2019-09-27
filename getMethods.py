# ------------------------------------TOOLBOX FILE FOR FEDORA-----------------------------------------#
import requests as rq
import subprocess as sp
import os
import sqlite3
import numpy as np
from database import FUSEKI_URL, ALIGN_PATH, SEQ_PATH, amino_dict, DBD_TREE_PATH, FCREPO_URL
from parseMethods import parseturtle
from queryCreator import createquery


# Take a TFC Uri and return just the label
def tfc_uri_to_label(tfc_id):
    tfc_id = tfc_id.split('/')
    if tfc_id[len(tfc_id)-1] == '':
        tfc_id = tfc_id[len(tfc_id)-2]
    else:
        tfc_id = tfc_id[len(tfc_id)-1]
    return tfc_id


# Create a query that GETs JSON data concerning a single TFClass node -> Mostly useful for getting URIs pertaining
# to certain identifiers or relations
def get_node(tfc_id):
    query = createquery('SELECT ?subject ?pred ?object WHERE {\n' + '?subject sybig:tfclass_id "' + tfc_id + '" .\n'
                        + '?subject ?pred ?object}\n')
    res = rq.get(FUSEKI_URL, params={'format': 'json', 'query': query})
    data = res.json()
    # Return the URI to given TFClass ID aswell
    tmp = data['results']['bindings'][0]
    subject = tmp['subject']['value']

    return data, subject


# For a given TFClassID, download alignment (if it exists)
def get_alignment(tfc_id):
    if "http://" not in tfc_id:
        data, subject = get_node(tfc_id)
        # subject = subject.replace("localhost", "coyote")
    else:
        subject = tfc_id
        # subject = subject.replace("localhost", "coyote")
    tmp = rq.get(subject).text
    # Since ldp:contains is not passed to Fuseki, parse turtle data for it
    res = parseturtle(tmp, ["ldp:contains"], 0)
    for arguments in res:
        if "fasta.txt" in arguments:
            arguments = arguments.replace(" .", '')
            fnull = open(os.devnull, 'w')
            sp.call(['wget', '-P', ALIGN_PATH, arguments], stdout=fnull, stderr=sp.STDOUT, close_fds=True)
            arg_splitter = arguments.split('/')
            arg_splitter = arg_splitter[len(arg_splitter)-1]
            id_splitter = tfc_id.split('/')
            id_splitter = id_splitter[len(id_splitter)-1]
            os.rename(ALIGN_PATH + arg_splitter, ALIGN_PATH + id_splitter)
            return 1  # No need to investigate further
    return 1


# Takes an URI and returns all relations as a list (done through turtle, since Fuseki has holes)
def get_triples(uri):
    tmp = rq.get(uri).text
    res = parseturtle(tmp, [], 1)
    return res


# Get URIs of all children nodes of a given URI, return as list
def get_children_nodes(uri):
    returnlist = []
    res = get_triples(uri)
    for x in range(len(res)):
        if res[x][0] == "ldp:contains":
            returnlist.append(res[x][1])
    returnlist.sort()
    return returnlist


# Takes a bottom level TFClass ID and returns the TF label aswell as a list of related TF_protein URIs
def get_related_proteins(tfc_id):
    data, subject = get_node(tfc_id)
    tmp = rq.get(subject).text
    label = parseturtle(tmp, ["rdfs:label"], 0)
    res = parseturtle(tmp, ["sybig:contains"], 0)
    res.sort()
    return label, res


# Takes a bottom level TFProtein URI and returns the TFClass ID connected to it (ID, URI and label)
def get_related_tfc(protein_uri):
    tmp = rq.get(protein_uri).text
    res1 = parseturtle(tmp, ["sybig:belongs_to"], 0)
    if res1:
        res1 = res1[0]
        tmp = rq.get(res1).text
        if tmp:
            res2 = parseturtle(tmp, ["rdfs:label", "sybig:tfclass_id"], 0)
            if len(res2) >= 2:
                res2.sort()
                res = [res1, res2[0], res2[1]]
            elif len(res2) == 1:
                res = [res1, res2[0]]
            else:
                id_splitter = res1.split('\t')
                id_splitter = id_splitter[len(id_splitter)-1]
                res = [res1, id_splitter]
            return res
        else:
            return [res1, "None", "None"]
    else:
        return []


# Collects all TFClass URIs (mode 0) or IDs (mode 1) in the repository
def collect_tfc_ids(mode, level):
    tf_list = []
    if level == "All":
        query = createquery('SELECT ?subject WHERE {?subject sybig:tfclass_id ?o}')
    else:
        query = createquery('SELECT ?subject WHERE {?subject sybig:tfclass_id ?o ; sybig:level "' + str(level) + '" .}')
    tf_list_unparsed = rq.get(FUSEKI_URL, params={'format': 'json', 'query': query})
    tf_list_json = tf_list_unparsed.json()
    for x in tf_list_json['results']['bindings']:
        if mode == 0:
            element = x['subject']['value']  # .replace("localhost", "coyote")
            tf_list.append(element)
        else:
            temp_ids = x['subject']['value']
            id_splitter = temp_ids.split('/')
            id_splitter = id_splitter[len(id_splitter) - 1]
            tf_list.append(id_splitter)
    tf_list.sort()
    return tf_list


# Collects all TFProtein URIs through a query
def collect_tf_proteins():
    prot_list = []
    query = createquery('SELECT ?subject WHERE { ?subject sybig:level ?o FILTER ('
                        ' !EXISTS {	?subject sybig:tfclass_id ?o2})}')
    prot_list_unparsed = rq.get(FUSEKI_URL, params={'format': 'json', 'query': query})
    prot_list_json = prot_list_unparsed.json()
    for x in prot_list_json['results']['bindings']:
        prot_list.append(x['subject']['value'])
    prot_list.sort()
    return prot_list


# Collects all dbds and sequences from Fedora (must use turtle files for parsing)
def collect_dbd_seqs(mode):
    # Runs for a while (20+ minutes)
    # Mode 0: Just return. Mode 1: Write to file. Mode 2: Save to DB (recommended)
    big_seq_list = []
    big_dbd_list = []
    if mode == 2:
        conn = sqlite3.connect('Sequences.db')
        cu = conn.cursor()
    p_list = collect_tf_proteins()
    dbd_path = SEQ_PATH + "Fedora_DBDs/"
    fseq_path = SEQ_PATH + "Fedora_Seqs/"
    for nodes in p_list:
        # nodes = nodes.replace("localhost", "coyote")
        related_tfc = get_related_tfc(nodes)
        if not related_tfc:
            related_tfc = ["-", "None found", "-"]
        if "http" in related_tfc:
            print(related_tfc)
        tmp = rq.get(nodes).text
        dbds = parseturtle(tmp, ["sybig:dbd"], 0)
        seqs = parseturtle(tmp, ["sybig:sequence"], 0)
        nodesplitter = nodes.split('/')
        nodesplitter = nodesplitter[len(nodesplitter)-1]
        # print(nodesplitter)
        if mode < 2:
            for x in dbds:
                big_dbd_list.append([nodesplitter, x])
            for y in seqs:
                big_seq_list.append([nodesplitter, y])
        else:
            if "http" in related_tfc[1]:
                related_tfc[1] = "NONE"
            if dbds and not seqs:
                cu.execute("INSERT INTO SEQUENCES VALUES('" + str(nodesplitter) + "', '" + str(related_tfc[1]) + "', '"
                           + str(dbds[0]) + "', '" + "None" + "')")
            elif seqs and not dbds:
                cu.execute("INSERT INTO SEQUENCES VALUES('" + str(nodesplitter) + "', '" + str(related_tfc[1]) + "', '"
                           + "None" + "', '" + str(seqs[0]) + "')")
            elif seqs and dbds:
                cu.execute("INSERT INTO SEQUENCES VALUES('" + str(nodesplitter) + "', '" + str(related_tfc[1]) + "', '"
                           + str(dbds[0]) + "', '" + str(seqs[0]) + "')")

    if mode == 1:
        with open("DBDs_current", "w") as dbd:
            for d in big_dbd_list:
                dbd.write(d[0] + "\t" + d[1] + "\n")
        os.rename("DBDs_current", dbd_path + "DBDs_current")
        with open("SEQs_current", "w") as seq:
            for s in big_seq_list:
                seq.write(s[0] + "\t" + s[1] + "\n")
        os.rename("SEQs_current", fseq_path + "SEQs_current")
    if mode == 2:
        conn.commit()
        conn.close()
    return big_dbd_list, big_seq_list


# Returns a "normalized" descriptor for a TFClass protein name, without addendums such as "DBD_ma5"
def normalize_prot_name(line):
    protein = line[1:len(line) - 1]
    if protein.find('-annot') != -1:
        pos = protein.find('-annot')
        protein = protein[0:pos]
    else:
        protein_tmp = protein.split('_')
        protein = ''
        for j in range(0, len(protein_tmp) - 1):
            protein += protein_tmp[j]
            if j != len(protein_tmp) - 2:
                protein += '_'
    protein = protein.replace("-DBD", "")
    protein = protein.replace("_ma", "")
    if protein.find('annot') != -1:
        pos = protein.find('annot')
        protein = protein[0:pos]
    return protein


# Returns a list of alignments currently in the Alignments folder
def alignment_list(level):
    return_list = []
    for root, dirs, files in os.walk(ALIGN_PATH):
        for file in files:
            if level == 0:
                return_list.append(file)
            else:
                filex = file.split('.')
                if len(filex) == level:
                    return_list.append(file)
    return_list.sort()
    return return_list


# Calculates entropy according to Shannon for all alignments, saves to DB
# Currently counts gaps
def entropy_to_db():
    al_list = []
    al_nr = 0

    def calc_entropy(pi):
        col_entropy = 0.0
        for probability in pi:
            if probability != 0:
                col_entropy += (probability * np.log2(probability))
        return col_entropy

    conn = sqlite3.connect("Sequences.db")
    cu = conn.cursor()
    for alignment in alignment_list(0):
        with open(ALIGN_PATH + alignment, "r") as align_file:
            alignment_entropy = 0.0
            # Get alignment length
            o_align_file = align_file.readlines()
            shape = (21, len(o_align_file[1]))
            entropy_matrix = np.zeros(shape)
            line_nr = 0
            for line in o_align_file:
                line = line.replace('\n', '')
                if line[0] != '>':
                    line_nr += 1
                    for pos in range(0, len(line)):
                        entropy_matrix[amino_dict[line[pos]], pos] += 1
            entropy_matrix /= line_nr
            for i in range(0, len(o_align_file[1])):
                alignment_entropy += calc_entropy(entropy_matrix[:, i])
                # print("Entropy at position " + str(i) + " is " + str(calc_entropy(entropy_matrix[:, i])))
            alignment_entropy = -alignment_entropy / len(o_align_file[1])
            print(alignment + " -> " + str(alignment_entropy))
            al_list.append([alignment, alignment_entropy])

    al_list.sort()
    for al in al_list:
        al_nr += 1
        cu.execute("INSERT INTO ALIGNMENTS VALUES ('" + str(al_nr) + "', '" + str(al[0]) + "',  '" + str(al[1]) + "')")
    conn.commit()
    conn.close()


# Saves alignments to DBs, gaps are ignored in entropy calculation
def entropy_to_db_ignoregaps():
    al_list = []
    al_nr = 0

    def calc_entropy(pi):
        col_entropy = 0.0
        for probability in pi:
            if probability != 0:
                col_entropy += (probability * np.log2(probability))
        return col_entropy

    conn = sqlite3.connect("Sequences.db")
    cu = conn.cursor()
    for alignment in alignment_list(0):
        with open(ALIGN_PATH + alignment, "r") as align_file:
            alignment_entropy = 0.0
            # Get alignment length
            o_align_file = align_file.readlines()
            shape = (21, len(o_align_file[1])-1)
            entropy_matrix = np.zeros(shape)
            line_nr = 0
            for line in o_align_file:
                line = line.replace('\n', '')
                if line[0] != '>':
                    line_nr += 1
                    for pos in range(0, len(line)):
                        if line[pos] != 'X' and line[pos] != '-':
                            entropy_matrix[amino_dict[line[pos]], pos] += 1
            sums = entropy_matrix.sum(axis=0)
            # print(sums)
            for x in range(0, len(sums)):
                # print(entropy_matrix[:, x])
                counted = sums[x]
                if counted != 0:
                    entropy_matrix[:, x] /= counted

            for i in range(0, len(o_align_file[1])-1):
                alignment_entropy += calc_entropy(entropy_matrix[:, i])
                # print("Entropy at position " + str(i) + " is " + str(calc_entropy(entropy_matrix[:, i])))
            alignment_entropy = -alignment_entropy / len(o_align_file[1])
            print(alignment + " -> " + str(alignment_entropy))
            al_list.append([alignment, alignment_entropy])

    al_list.sort()
    cu.execute("DELETE FROM ALIGNMENTS")
    for al in al_list:
        al_nr += 1
        cu.execute("INSERT INTO ALIGNMENTS VALUES ('" + str(al_nr) + "', '" + str(al[0]) + "',  '" + str(al[1]) + "')")
    conn.commit()
    conn.close()


# Calculates entropy for one alignment (does not work with varying lenghts in a file)
def entropy_alignment(alignment_id):

    def calc_entropy(pi):
        col_entropy = 0.0
        for probability in pi:
            if probability != 0:
                col_entropy += (probability * np.log2(probability))
        return col_entropy

    with open(ALIGN_PATH + alignment_id, "r") as align_file:
        alignment_entropy = 0.0
        # Get alignment length
        o_align_file = align_file.readlines()
        shape = (21, len(o_align_file[1]))
        entropy_matrix = np.zeros(shape)
        line_nr = 0
        for line in o_align_file:
            line = line.replace('\n', '')
            if line[0] != '>':
                line_nr += 1
                for pos in range(0, len(line)):
                    entropy_matrix[amino_dict[line[pos]], pos] += 1
        entropy_matrix /= line_nr
        for i in range(0, len(o_align_file[1])):
            alignment_entropy += calc_entropy(entropy_matrix[:, i])
            # print("Entropy at position " + str(i) + " is " + str(calc_entropy(entropy_matrix[:, i])))
        alignment_entropy = -alignment_entropy / len(o_align_file[1])
    return alignment_entropy


# Creates an unaligned FASTA file containing DBDs for each genus using the repository (deprecated)
def create_genera_alignments_repo():
    all_tfs = collect_tfc_ids(0, "Genus")
    for tfc in all_tfs:
        current_tfc = rq.get(tfc).text
        child_list = parseturtle(current_tfc, ["sybig:contains"], 0)
        alignlist = []
        for child in child_list:
            current_child = rq.get(child).text
            protname = parseturtle(current_child, ["sybig:name"], 0)[0]
            dbd = parseturtle(current_child, ["sybig:dbd"], 0)
            if dbd:
                alignlist.append('>'+str(protname)+'\n')
                alignlist.append(dbd[0] + '\n')
        if alignlist:
            with open(ALIGN_PATH + "Unaligned_" + tfc_uri_to_label(tfc), 'w') as current_align:
                for line in alignlist:
                    current_align.write(line)


# Reduces TFC ID by one level
def reduce_tfc_id(tfc):
    if tfc == '':
        return ''
    reduced = False
    cutoff = len(tfc) - 1
    for i in range(len(tfc) - 1, 0, -1):
        if reduced:
            return tfc[0:cutoff]
        else:
            if tfc[i] != '.':
                cutoff -= 1
            else:
                reduced = True


# Creates genera alignments from subfamily alignments
def create_genera_alignments():

    def protein_query(protname):
        protein_to_search = tfc_label = ''
        query = createquery('SELECT ?subject WHERE {?subject sybig:name "' + str(protname) + '"}')
        # print("Checking where " + protname + "belongs to")
        uri_unparsed = rq.get(FUSEKI_URL, params={'format': 'json', 'query': query})
        uri = uri_unparsed.json()
        for x in uri['results']['bindings']:
            protein_to_search = x['subject']['value']
            # print("Found " + protein_to_search)
            query = createquery('SELECT ?subject WHERE {<' + str(protein_to_search) + '> sybig:belongs_to ?subject}')
            tfc_unparsed = rq.get(FUSEKI_URL, params={'format': 'json', 'query': query})
            tfc = tfc_unparsed.json()
            for y in tfc['results']['bindings']:
                tfc_label = tfc_uri_to_label(y['subject']['value'])
        return protein_to_search, tfc_label

    # Check which alignment files we need to look at to create genera alignments - Start with family level
    alignlist = alignment_list(4)
    alignlist_family = alignment_list(3)
    # First, add all sequences from subfamily alignment files to their respective files
    for alignfile in alignlist:
        with open(ALIGN_PATH + alignfile, 'r') as potential:
            deviation = "_1"
            currentlines = potential.readlines()
            for line in currentlines:
                just_created = False
                if line[0] == '>':
                    protein_name = line[1:len(line)-1]
                    protein_name = normalize_prot_name(">"+protein_name)
                else:
                    dbd = line
                    pts, tfcl = protein_query(protein_name)
                    if pts != '' and tfcl != '' and tfcl != 'rest':
                        if os.path.exists(ALIGN_PATH + tfcl):
                            with open(ALIGN_PATH + tfcl, 'r') as checklength:
                                cl = checklength.readlines()
                                if len(cl) >= 2:
                                    alignlength = len(cl[1])
                                else:
                                    alignlength = -1
                        else:
                            alignlength = -1
                        with open(ALIGN_PATH + tfcl, 'a') as add_to_this:
                            if len(dbd) == alignlength or alignlength == -1:
                                add_to_this.write(">" + protein_name + "\n")
                                add_to_this.write(dbd)
                            else:
                                if not os.path.exists(ALIGN_PATH + tfcl + deviation):
                                    create = open(ALIGN_PATH + tfcl + deviation, 'w')
                                    create.close()
                                with open(ALIGN_PATH + tfcl + deviation, 'a') as add_dev:
                                    add_dev.write(">" + protein_name + "\n")
                                    add_dev.write(dbd)

    # Then, fill in gaps by checking the family alignment files
    for alignfile in alignlist_family:
        with open(ALIGN_PATH + alignfile) as current:
            currentlines = current.readlines()
            for line in currentlines:
                if line[0] == '>':
                    protein_name = line[1:len(line)-1]
                    protein_name = normalize_prot_name(">"+protein_name)
                else:
                    dbd = line
                    pts, tfcl = protein_query(protein_name)
                    if tfcl != '' and tfcl != 'rest':
                        if not os.path.exists(ALIGN_PATH + tfcl):
                            create = open(ALIGN_PATH + tfcl, 'w')
                            create.close()
                        with open(ALIGN_PATH + tfcl, 'r') as check:
                            # print("Checking file " + tfcl)
                            cl = check.readlines()
                            # print(len(cl))
                            if len(cl) >= 2:
                                dbd_len = len(cl[1])  # Sequence length of alignment file
                            else:
                                dbd_len = len(dbd)
                            found = False
                            for line2 in cl:
                                if protein_name in line2:
                                    found = True
                            if not found:
                                if tfcl != '' and pts != '' and tfcl != 'rest' and len(dbd) == dbd_len:
                                    with open(ALIGN_PATH + tfcl, 'a') as add_to_this:
                                        print("Adding a DBD of length " + str(len(dbd)) + " to alignment " + tfcl + " with length " + str(dbd_len))
                                        add_to_this.write(">" + protein_name + "\n")
                                        add_to_this.write(dbd)


# Adds DBDs to all corresponding nodes in Fedora
def add_dbds():

    def add_single_dbd(dbd):
        dbd = dbd.rstrip()
        query = createquery('SELECT ?subject WHERE { ?subject <http://sybig.de/tfclass#name> "' +
                            protein.strip() + '"}')
        query_out = rq.get(FUSEKI_URL, params={'format': 'json', 'query': query})
        query_out_json = query_out.json()
        for element in query_out_json['results']['bindings']:
            protein_uri = element['subject']['value']
            # protein_uri = protein_uri.replace("localhost", "coyote")
            turtlefile = rq.get(protein_uri).text
            print(protein_uri + " ----- " + dbd)
            with open("temp_turtle", "w") as update:
                update.write(turtlefile)
                update.write('\n<' + protein_uri + '> sybig:dbd "' + dbd + '" .')
            os.system('curl -X PUT -H "Content-type: text/turtle" --data-binary "@temp_turtle" "' + protein_uri + '"')

    fasta_l = open("fasta_dbd_list", "r")
    fasta_l_text = fasta_l.readlines()
    fasta_l_text.sort()
    fasta_l.close()

    for file in fasta_l_text:
        file = file.strip()
        file = file.replace("/sybig/home/cservices/tfclass3/dbd/dbd_tree/", DBD_TREE_PATH)
        with open(file, "r") as current_dbd_file:
            protein = ""
            current_dbd = current_dbd_file.readlines()
            for line in current_dbd:
                if line == '':
                    continue
                # print(line)
                if line[0] == '>':
                    line = line.split('\t')
                    if len(line) == 1:
                        protein = normalize_prot_name(line[0])
                    else:
                        protein = normalize_prot_name(line[0])
                        add_single_dbd(line[1])
                else:
                    add_single_dbd(line)


def count_genera_alignments():
    al = alignment_list(5)
    print(len(al))


def add_uniprots(file):
    with open(file, 'r') as uplist:
        upl = uplist.readlines()
        for row in upl:
            row_s = row.split(';')
            os.system('curl -H "Accept:text/turtle" -X GET "' + str(row_s[0]) + '" >node.rdf')
            with open("node.rdf", 'a') as addhere:
                addhere.write('<' + str(row_s[0]) + '> sybig:xref "UNIPROT:' + str(row_s[1].strip()) + '" .\n')
            os.system('curl -X PUT -H "Content-Type: text/turtle" --data-binary "@node.rdf" "' + str(row_s[0]) + '"')


# ali_list = collect_tfc_ids(0, 'All')
# for element in ali_list:
#    get_alignment(element)
# create_genera_alignments()
# entropy_to_db_ignoregaps()
# count_genera_alignments()
# add_uniprots("Missing_Uniprots")








