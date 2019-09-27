def createquery(querystring):
    return ('PREFIX fedora: <http://fedora.info/definitions/v4/repository#>\n'
            + 'PREFIX sybig: <http://sybig.de/tfclass#>\n'
            + 'PREFIX ldp: <http://www.w3.org/ns/ldp#>\n'
            + 'PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>\n'
            + querystring)
