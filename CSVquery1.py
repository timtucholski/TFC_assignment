import requests as rq
from queryCreator import createquery
from database import FUSEKI_URL

with open('TF_ids.tsv', "w") as tsvfile:
    tf_list = []
    query = createquery('SELECT ?subject WHERE {?subject sybig:tfclass_id ?o}')
    tf_list_unparsed = rq.get(FUSEKI_URL, params={'format': 'json', 'query': query})
    tf_list_json = tf_list_unparsed.json()
    # print(tf_list_json)
    for x in tf_list_json['results']['bindings']:
        tf_list.append(x['subject']['value'])

    tf_list.sort()
    tsvfile.write("TFClass_ID\tTF_Label\tConsensus Sequence\n")
    res_list = []
    last_triple = []
    for tfc in tf_list:
        # print("Found " + tfc)
        query = createquery('SELECT ?n ?o ?i WHERE\n {<' + tfc + '> sybig:consensus_sequence ?i . \n'
                            + '<' + tfc + '> rdfs:label ?o . \n'
                            + '<' + tfc + '> sybig:tfclass_id ?n}')
        # print("Created query: " + query)
        res = rq.get(FUSEKI_URL, params={'format': 'json', 'query': query})
        data = res.json()
        if data['results']['bindings']:  # Non-empty result
            for tmp in data['results']['bindings']:
                if "http://" not in tmp['n']['value'] and tmp != last_triple:
                    last_triple = tmp
                    res_list.append([tmp['n']['value'], tmp['o']['value'], tmp['i']['value']])

        else:
            query = createquery('SELECT ?n ?o  WHERE\n {<' + tfc + '> sybig:tfclass_id ?n .\n <' + tfc + '> rdfs:label ?o}')
            res = rq.get(FUSEKI_URL, params={'format': 'json', 'query': query})
            data = res.json()
            for newc in data['results']['bindings']:
                tf_name = newc['n']['value']
                tf_label = newc['o']['value']
                if "http://" not in tf_name and newc != last_triple:
                    last_triple = newc
                    res_list.append([tf_name, tf_label, ''])
    for trp in res_list:
        print(trp)
        current_line = trp[0]+'\t'+trp[1]+'\t'+trp[2]+'\n'
        tsvfile.write(current_line)
