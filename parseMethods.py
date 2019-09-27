

def parseturtle(turtlefile, predicates, mode):
    # Parse a turtle file, returning all objects that correspond to any of the given predicates
    returnlist = []
    for lines in turtlefile.splitlines():
        args = lines.split(' ')
        args = [x for x in args if (x != ';' and x != '')]

        # Keep strings intact
        if len(args) >= 3:
            for tmp in args[2:len(args)]:
                args[1] = args[1] + ' ' + tmp
        args = args[:2]

        # Remove incomplete triples (should only be the subject line) and prefixes
        if not args or len(args) == 1 or "prefix" in args[0]:
            continue

        # Split off data type
        temp = args[1].split("^^")
        args[1] = temp[0]

        # Remove parentheses
        args[1] = args[1].replace('<', '')
        args[1] = args[1].replace('>', '')

        # Remove unneeded "
        temp = args[1]
        if temp[0] == '"':
            temp = temp[1:len(temp)]
        if temp[len(temp)-1] == '"':
            temp = temp[:len(temp)-1]
        args[1] = temp

        # Append results to list according to mode
        # mode 0: Only append objects. mode 1: Append predicates and objects. mode 2: only return line number
        if not predicates:
            if mode == 0:
                returnlist.append(args[1])
            elif mode == 1:
                returnlist.append([args[0], args[1]])
            else:
                returnlist.append(turtlefile.splitlines().index(lines))
        else:
            for preds in predicates:
                if args[0] == preds:
                    if not (args[0] == "rdfs:label" and "http://" in args[1]):  # No shitty IDs
                        if mode == 0:
                            returnlist.append(args[1])
                        elif mode == 1:
                            returnlist.append([args[0], args[1]])
                        else:
                            returnlist.append(turtlefile.splitlines().index(lines))
    return returnlist
