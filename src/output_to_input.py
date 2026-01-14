def out_to_in(folder, filename, new_file):
    count = 0
    empty_count = 0
    full_count = 1
    i = 0
    counts = []
    empties = [[]]
    fulls = [[]]
    g = open(folder+'/out_copy.vt3d', "w")
    with open('{}/{}'.format(folder, filename), 'r') as f:
        for line in f:
            l = line.split()
            if l:
                if l[0] != '0':
                    count += 1
                    empties[i].append(empty_count)
                    fulls[i].append(empty_count + full_count)
                    if l[0] == '1':
                        l.pop(0)
                    if len(counts) > 3:
                        g.writelines([str(len(l))+'\t'])
                    g.writelines([x+'\t' for x in l])
                    g.writelines(['\n'])
                    full_count += 1
                else:
                    empty_count += 1
            else:
                counts.append(count)
                empties.append([])
                fulls.append([])
                i += 1
                count = 0
                empty_count = 0
                full_count = 1
                g.writelines(['\n'])
    counts.append(count)
    g.close()

    counts.pop(0)
    counts.pop(0)
    counts = [str(m) for m in counts]

    entity_count = 0
    k = open('{}/{}'.format(folder, new_file), "w")
    with open(folder+'/out_copy.vt3d', 'r') as g:
        for line in g:
            l = line.split()
            if l:
                if entity_count == 0:
                    k.writelines([x+'\t' for x in counts])
                    k.writelines(['\n'])
                if entity_count == 1 or entity_count == 2:
                    k.writelines(line)
                if entity_count == 3:
                    l_num = [int(x) for x in l]
                    l_new = []
                    for m in range(len(l_num)):
                        for n in range(len(fulls[2])):
                            if l_num[m] == fulls[2][n]:
                                l_new.append(l_num[m] - empties[2][n])
                    k.writelines([str(x)+'\t' for x in l_new])
                    k.writelines(['\n'])
                if entity_count == 4 or entity_count == 5:
                    l_num = [int(x) for x in l]
                    l_new = [l_num[0]]
                    for m in range(1, len(l_num)):
                        for n in range(len(fulls[entity_count-1])):
                            if abs(l_num[m]) == fulls[entity_count-1][n]:
                                if l_num[m] > 0:
                                    l_new.append(l_num[m] - empties[entity_count-1][n])
                                if l_num[m] < 0:
                                    l_new.append(l_num[m] + empties[entity_count-1][n])
                    k.writelines([str(x)+'\t' for x in l_new])
                    k.writelines(['\n'])
            else:
                entity_count += 1
                k.writelines(['\n'])
    k.close()

