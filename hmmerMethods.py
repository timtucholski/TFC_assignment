from database import HMMER_PATH, PROFILE_PATH, ALIGN_PATH, SEQ_PATH, RESULT_PATH
from getMethods import tfc_uri_to_label
from Visualizer import vis_best_evalues
from datetime import datetime
from DBhandler import get_sequences_from_db, get_profiled_sequences
import os
import subprocess as sp
fnull = open(os.devnull, 'w')

MASTER_PROFILE = PROFILE_PATH + "MasterProfile"
PRESSED_MASTER_PROFILE = PROFILE_PATH + "MasterProfile_pressed/"


# A method that builds a singularprofile from a given alignment file using hmmbuild
def build_profile(file):
    sp.call([HMMER_PATH + 'hmmbuild', PROFILE_PATH + 'TestFile_out', ALIGN_PATH + file], stdout=fnull,
            stderr=sp.STDOUT, close_fds=True)


# A method that builds a master profile database from all current alignment files in ALIGN_PATH
def build_master():
    alignment_list = []
    for root, dirs, files in os.walk(ALIGN_PATH):
        for file in files:
            alignment_list.append(os.path.join(root, file))
    alignment_list = sorted(alignment_list)
    master_profile = open(PROFILE_PATH + "MasterProfile", 'w')
    for alignment in alignment_list:
        al_nr = alignment.split('/')
        al_nr = al_nr[len(al_nr)-1]
        sp.call([HMMER_PATH + 'hmmbuild', '-n', al_nr,  PROFILE_PATH + al_nr, alignment], stdout=fnull,
                stderr=sp.STDOUT, close_fds=True)
        current_rl = open(PROFILE_PATH + al_nr, 'r')
        current = current_rl.readlines()
        for line in current:
            master_profile.write(line)
        current_rl.close()
    master_profile.close()


# A method that combines all orthologues files into one and tags them accordingly for testing
def orthologue_master():
    fasta_list = []
    for root, dirs, files in os.walk(SEQ_PATH + "Orthologs/"):
        for file in files:
            fasta_list.append(os.path.join(root, file))
    master_file_lines = []
    for file in fasta_list:
        if "Master" in file:
            continue
        with open(file, 'r') as current:
            cl = current.readlines()
            for line in cl:
                if line[0] == '>':
                    line_split = line.split('|')
                    hitdescriptor = line_split[3].replace('[', '').replace(']', '').replace('\n', '')
                    hitdescriptor = hitdescriptor.strip()
                    hitdescriptor = hitdescriptor.replace(' ', '_')
                    new_line = ">" + hitdescriptor + "---" + tfc_uri_to_label(file)
                    new_line = new_line.strip()
                    new_line = new_line + "\n"
                    master_file_lines.append(new_line)
                else:
                    master_file_lines.append(line)
    with open(SEQ_PATH + "Orthologs/Master", 'w') as master:
        for line in master_file_lines:
            master.write(line)


# A method that prepares a profile database for hmmscan (hmmpress)
def press_profile(name):
    sp.call([HMMER_PATH + 'hmmpress', '-f', PROFILE_PATH + name], stdout=fnull, stderr=sp.STDOUT, close_fds=True)
    # os.mkdir(PROFILE_PATH + name + '_pressed/')
    os.rename(PROFILE_PATH + name + '.h3i', PROFILE_PATH + name + '_pressed/' + name + '.h3i')
    os.rename(PROFILE_PATH + name + '.h3p', PROFILE_PATH + name + '_pressed/' + name + '.h3p')
    os.rename(PROFILE_PATH + name + '.h3f', PROFILE_PATH + name + '_pressed/' + name + '.h3f')
    os.rename(PROFILE_PATH + name + '.h3m', PROFILE_PATH + name + '_pressed/' + name + '.h3m')


# A method that searches all sequences from file X against profile
def match_sequence(sequence_file, profile_name):
    now = datetime.now()
    dt_string = now.strftime("%d_%m_%Y %H_%M_%S")
    tpath = RESULT_PATH + dt_string + "/"
    os.mkdir(tpath)
    sp.call([HMMER_PATH + 'hmmscan', '-o', tpath + "MainOutput", '--tblout', tpath + "Tbl", '--domtblout', tpath + "DomTbl", '--pfamtblout', tpath + "PfamTbl", PRESSED_MASTER_PROFILE + profile_name, SEQ_PATH + sequence_file], stdout=fnull, stderr=sp.STDOUT, close_fds=True)
    # evaluate_best_hits(dt_string)


# A method that searches all sequences from file X against profile, also outputs errors to console
def match_sequence_with_output(sequence_file, profile_name):
    now = datetime.now()
    dt_string = now.strftime("%d_%m_%Y %H_%M_%S")
    tpath = RESULT_PATH + dt_string + "/"
    os.mkdir(tpath)
    sp.call([HMMER_PATH + 'hmmscan', '-o', tpath + "MainOutput", '--tblout', tpath + "Tbl", '--domtblout', tpath + "DomTbl", '--pfamtblout', tpath + "PfamTbl", PRESSED_MASTER_PROFILE + profile_name, SEQ_PATH + sequence_file])
    # evaluate_best_hits(dt_string)


# A method that parses HMMER output (Tbl)
def parse_hmmer_output(tbl_file):
    result_list = []
    with open(tbl_file, 'r') as hmmer_results:
        hmmer_results_rl = hmmer_results.readlines()
        for line in hmmer_results_rl:
            if line[0] != '#':
                current_line = line.strip().split(' ')
                current_line = [x for x in current_line if (x != ';' and x != '')]
                result_list.append([current_line[2], current_line[0], current_line[4], current_line[5], current_line[6],
                                    current_line[7], current_line[8], current_line[9]])
    return result_list


# A method that writes a list of DB sequences to a file (for testing)
def list_to_file(seq_list, name):
    detect = len(seq_list[0])
    with open(name, 'w') as ts:
        for row in seq_list:
            if "http" in row[1]:
                row[1] = "NONE"
            if detect == 3:
                ts.write(">" + row[0] + "---" + row[1] + "\n")
                ts.write(row[2] + "\n")
            elif detect == 4:
                ts.write(">" + row[0] + "---" + row[1] + "(DBD)" + "\n")
                ts.write(row[2] + "\n")
                ts.write(">" + row[0] + "---" + row[1] + "(SEQ)" + "\n")
                ts.write(row[3] + "\n")
    os.rename(name, SEQ_PATH + name)


# A method that builds a table out of the best hits from a HMMER Tbl output (deprecated)
def best_hit_list():
    bhl = []
    reslist = parse_hmmer_output(RESULT_PATH + "Tbl")
    current = ""
    current_best_hit = ""
    current_best_hit_second = ""
    sec = True
    index = 0
    ev = -1.0
    dist = -1.0
    minv = -1.0
    category = 0
    oneliner = False
    writeout = open(RESULT_PATH + "currentout.tsv", "w")
    writeout.write("SeqDescriptor\tHitDescriptor\tDistFactor\thmmer Score\tCategory\n")
    for x in reslist:
        if x[0] != current:
            if (ev != -1.0 and dist != -1.0 and not oneliner) or (index == len(reslist)-1 and not oneliner):
                bhl.append([current, current_best_hit, dist, ev, category])
                writeout.write(current+"\t"+current_best_hit + " / " + current_best_hit_second + "\t"+str(dist)+"\t"+str(ev)+"\t"+str(category)+"\n")
                # print("Appending: " + current + " which hit " + current_best_hit + " with an ev of " + str(ev) + " and distance " + str(dist) + " and category " + str(category))
            current = x[0]
            current_best_hit = x[1]
            oneliner = True
            # Get Score
            ev = float(x[3])
            sec = True
            # Get category
            csp = current.split("---")
            csp = csp[1].split("(")
            csp = csp[0]
            s_hits = x[1].split(".")
            hits = csp.split(".")
            limit = min(len(s_hits), len(hits))
            category = 0
            for i in range(0, limit):
                if hits[i] != s_hits[i]:
                    category += 1
                    i = limit
        else:
            if sec:
                dist = float(x[3])/ev
                sec = False
                oneliner = False
                current_best_hit_second = x[1]
            if index != len(reslist)-1 and reslist[index+1][0] != current:
                minv = float(x[3])
            if index == len(reslist)-2:
                minv = float(reslist[index+1][3])
        index += 1
    writeout.close()
    vis_best_evalues(bhl)
    return bhl


# New method to evaluate hits (deprecated)
def evaluate_best_hits(dat):

    hitlist = []

    # Matches two TFC IDs
    def match_tfc_ids(id1, id2):
        category = 0
        id1s = id1.split('.')
        id2s = id2.split('.')
        match_length = min(len(id1s), len(id2s))
        for i in range(0, match_length):
            if id1s[i] == id2s[i]:
                category += 1
            else:
                if category == 0:
                    print(id1 + " " + id2 + " - Nothing hit")
                    return "Nothing hit"
                elif category == 1:
                    print(id1 + " " + id2 + " - Only superclass hit")
                    return "Superclass hit"
                elif category == 2:
                    print(id1 + " " + id2 + " - Class hit")
                    return "Class hit"
                elif category == 3:
                    return "Family hit"
                elif category == 4:
                    return "Subfamily hit"
                elif category == 5:
                    return "Genus hit"
        return "Genus hit"  # Only reached when genus hit

    # Convert HMMER output to a python list
    result_list = parse_hmmer_output(RESULT_PATH + dat + "/" + "Tbl")
    # Sort hits into categories: Genus hit, Subfamily hit, Family hit, Class hit, Superclass hit, nothing hit
    # Save e-value, bias, score for both the best sequence and the best domain

    current_prot = ''
    score = 0.0
    dom_score = 0.0
    sec = False
    prot_hit = ''
    current_prot_desc = ''

    for line in result_list:
        if sec:
            sec = False
            dist = float(line[3]) / score
            domdist = float(line[6]) / dom_score
            fdist = score - float(line[3])
            domfdist = dom_score - float(line[6])
            if current_prot_desc != "NONE":
                hitlist.append([current_prot, prot_hit, dist, fdist, domdist, domfdist, match_tfc_ids(prot_hit,
                                                                                                      current_prot_desc)])
        if current_prot != line[0]:  # New protein
            current_prot = line[0]
            sec = True
            prot_hit = line[1]
            score = float(line[3])
            dom_score = float(line[6])
            current_prot_desc = current_prot.split('---')
            current_prot_desc = current_prot_desc[1]

    vis_best_evalues(hitlist)


# build_master()
# press_profile("MasterProfile")
# match_sequence('Orthologs/1/1.1.1.2.2', "MasterProfile")
# parse_hmmer_output('/sybig/home/ttu/BSc/ScanResults/Tbl')
# best_hit_list()
# orthologue_master()


