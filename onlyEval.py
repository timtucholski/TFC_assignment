from evaluationMethods import evaluate_results, median_entropy

retl = evaluate_results("proteome_out_onlyprof.csv")
print(median_entropy("proteome_out_onlyprof.csv"))
