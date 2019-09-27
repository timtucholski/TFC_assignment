import sqlite3
import os
from database import ALIGN_PATH
from getMethods import normalize_prot_name


# Return sequences from the DB
def get_sequences_from_db(mode):
    # Mode 0: Return both DBDs and sequences. Mode 1: Return DBDs. Mode 2: Return sequences.
    conn = sqlite3.connect("Sequences.db")
    cu = conn.cursor()
    return_list = []
    if mode == 0:
        for row in cu.execute("SELECT * FROM SEQUENCES"):
            return_list.append([row[0], row[1], row[2], row[3]])
    elif mode == 1:
        for row in cu.execute("SELECT ProtID, TfcID, DBD FROM SEQUENCES WHERE DBD != 'None'"):
            return_list.append([row[0], row[1], row[2]])
    elif mode == 2:
        for row in cu.execute("SELECT ProtID, TfcID, SEQ FROM SEQUENCES WHERE SEQ != 'None'"):
            return_list.append([row[0], row[1], row[2]])
    return return_list


# Only returns sequences from which a profile was built
def get_profiled_sequences(mode):
    # Get sequences that fit to one of the alignments
    seq_list = []
    for root, dirs, files in os.walk(ALIGN_PATH):
        for file in files:
            with open(ALIGN_PATH + file, 'r') as current_alignment:
                for line in current_alignment:
                    if line[0] == '>':
                        protein = normalize_prot_name(line)
                        seq_list.append(protein)
    conn = sqlite3.connect("Sequences.db")
    cu = conn.cursor()
    return_list = []
    for candidate in seq_list:
        if mode == 0:
            for row in cu.execute("SELECT * FROM SEQUENCES WHERE ProtID = '" + candidate +"'"):
                return_list.append([row[0], row[1], row[2], row[3]])
        elif mode == 1:
            for row in cu.execute("SELECT ProtID, TfcID, DBD FROM SEQUENCES WHERE DBD != 'None' AND ProtID = '" + candidate +"'"):
                return_list.append([row[0], row[1], row[2]])
        elif mode == 2:
            for row in cu.execute("SELECT ProtID, TfcID, SEQ FROM SEQUENCES WHERE SEQ != 'None' AND ProtID = '" + candidate +"'"):
                return_list.append([row[0], row[1], row[2]])

    return return_list


# Return entropy connected to a TFClass identifier
def return_entropy(tfc_id):
    conn = sqlite3.connect("Sequences.db")
    cu = conn.cursor()
    cu.execute("SELECT Entropy FROM ALIGNMENTS WHERE TfcID = '" + str(tfc_id) + "'")
    retv = cu.fetchone()[0]
    return retv


# --------TESTING----------
# print(get_sequences_from_db(2))
# get_profiled_sequences()
# print(return_entropy("1.2.6.4"))
