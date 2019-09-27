from DBhandler import get_sequences_from_db, get_profiled_sequences
from hmmerMethods import list_to_file, match_sequence, match_sequence_with_output
from database import SEQ_PATH
# from getMethods import collect_dbd_seqs

# collect_dbd_seqs(2)
# test_sequences = get_profiled_sequences()
# list_to_file(test_sequences)
match_sequence_with_output("seq_prof_dbd", "MasterProfile")

