import networkx as nx
####################################################################################################
############################ USEFUL QUERY TEMPLATES ################################################
####################################################################################################
def qFIND(G, n1, n2):
    try:
        hos_of_n1 = list(G.successors(n1))
        hos_of_n2 = list(G.successors(n2))
    except nx.exception.NetworkXError:
        return None
    try:
        ho_common = list(set(hos_of_n1) & set(hos_of_n2))[0]
        return ho_common
    except IndexError:
        return None

def qVCEE_E(G, v, c, ea, eb):
    es_of_v = list(G.successors(v))
    ps_of_c = list(G.predecessors(c))
    es_of_c = []
    for p_of_c in ps_of_c:
        es_of_c += list(G.predecessors(p_of_c))
    es_common = list(set(es_of_c) & set(es_of_v))
    for e_common in es_common:
        if e_common != ea and e_common != eb:
            e = e_common
            break
    try:
        return e
    except UnboundLocalError:
        return None

def qVEEE_E(G, v, ea, eb, ec):
    es_of_v = list(G.successors(v))
    for e_of_v in es_of_v:
        if e_of_v != ea and e_of_v != eb and e_of_v != ec:
            e = e_of_v
            break
    try:
        return e
    except UnboundLocalError:
        return None
    

def qVEP_E(G, v, ea, p):
    es_of_p = list(G.predecessors(p))
    es_of_v = list(G.successors(v))
    ee = list(set(es_of_p) & set(es_of_v))
    if ee[0] == ea:
        e = ee[1]
    elif ee[1] == ea:
        e = ee[0]
    try:
        return e
    except UnboundLocalError:
        return None

def qFREE(G, label):
    free_nodes = list(G.successors("Free"))
    for free_node in free_nodes:
        if free_node[0] == label:
            G.remove_edge("Free", free_node)
            break
    return free_node

def qSIGNe(G, s, t, v, e):
    r_vs = list(G.get_edge_data(v, s))[0][1]
    r_ve = list(G.get_edge_data(v, e))[0][1]
    r_et = list(G.get_edge_data(e, t))[0][1]
    r = -r_vs*r_ve*r_et
    return r

####################################################################################################
#######################PATTERN-MATCHING & GRAPH-TRANSFORMATION QUERIES##############################
####################################################################################################
# ET
# This funciton performs pattern matching for ET transition
def patternMatching_ET(G, e_id):
    e7 = ("Edge", e_id)
    v1v2 = list(G.predecessors(e7))
    v1 = v1v2[0]
    v2 = v1v2[1]
    ps_of_e7 = list(G.successors(e7))
    cs_of_e7 = []
    for p_of_e7 in ps_of_e7:
        cs_of_e7 += list(G.successors(p_of_e7))
    c1 = cs_of_e7[0] # c1
    ps_of_c1 = list(G.predecessors(c1))
    p1p2 = list(set(ps_of_e7) & set(ps_of_c1))
    p1 = p1p2[0]
    p2 = p1p2[1]
    e1 = qVEP_E(G, v=v1, ea=e7, p=p1) # e1
    e2 = qVEP_E(G, v=v1, ea=e7, p=p2) # e2
    p4 = qFIND(G, n1=e1, n2=e2) # p4
    e3 = qVEP_E(G, v=v2, ea=e7, p=p1) # e3
    e4 = qVEP_E(G, v=v2, ea=e7, p=p2) # e4
    p5 = qFIND(G, n1=e3, n2=e4) # p5
    e5 = qVEEE_E(G, v=v1, ea=e1, eb=e2, ec=e7) # e5
    e6 = qVEEE_E(G, v=v2, ea=e3, eb=e4, ec=e7) # e6
    p3 = qFIND(G, n1=e5, n2=e6) # p3
    p6 = qFIND(G, n1=e2, n2=e5) # p6
    p8 = qFIND(G, n1=e1, n2=e5) # p8
    p9 = qFIND(G, n1=e3, n2=e6) # p9
    p7 = qFIND(G, n1=e4, n2=e6) # p7
    c2 = qFIND(G, n1=p8, n2=p9) # c2
    c3 = qFIND(G, n1=p6, n2=p7) # c3
    c4 = qFIND(G, n1=p6, n2=p8) # c4
    c5 = qFIND(G, n1=p7, n2=p9) # c5
    v3 = qFREE(G, label="Vertex") # v3
    e8 = qFREE(G, label="Edge") # e8
    e9 = qFREE(G, label="Edge") # e9
    p10 = qFREE(G, label="Polygon") # p10
    return v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5

# This funciton performs graph transformation for ET transition
def graphTransformation_ET(G, v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5):
    G.add_edge(v1, e3, list(G.get_edge_data(v2, e3))[0])
    G.add_edge(v1, e9, ("IS_PART_OF", -1))
    G.add_edge(v2, e2, list(G.get_edge_data(v1, e2))[0])
    G.add_edge(v2, e8, ("IS_PART_OF", 1))
    try:
        G.add_edge(v3, e5, list(G.get_edge_data(v1, e5))[0])
    except TypeError:
        pass
    try:
        G.add_edge(v3, e6, list(G.get_edge_data(v2, e6))[0])
    except TypeError:
        pass
    G.add_edge(v3, e8, ("IS_PART_OF", -1))
    G.add_edge(v3, e9, ("IS_PART_OF", 1))
    G.add_edge(e7, p4, ("IS_PART_OF", qSIGNe(G, s=e7 , t=p4, v=v2, e=e2)))
    G.add_edge(e7, p5, ("IS_PART_OF", qSIGNe(G, s=e7 , t=p5, v=v1, e=e3)))
    G.add_edge(e7, p10, list(G.get_edge_data(v1, e7))[0])
    try:
        G.add_edge(e8, p6, ("IS_PART_OF", qSIGNe(G, s=e8, t=p6, v=v2, e=e2)))
    except TypeError:
        pass
    try:
        G.add_edge(e8, p7, ("IS_PART_OF", qSIGNe(G, s=e8, t=p7, v=v3, e=e6)))
    except TypeError:
        pass
    G.add_edge(e8, p10, ("IS_PART_OF", 1))
    G.add_edge(e9, p10, ("IS_PART_OF", 1))
    try:
        G.add_edge(e9, p8, ("IS_PART_OF", qSIGNe(G, s=e9, t=p8, v=v3, e=e5)))
    except TypeError:
        pass
    try:
        G.add_edge(e9, p9, ("IS_PART_OF", qSIGNe(G, s=e9, t=p9, v=v1, e=e3)))
    except TypeError:
        pass
    try:
        G.add_edge(p10, c4, ("IS_PART_OF", qSIGNe(G, s=p10, t=c4, v=e7, e=p4)))
    except TypeError:
        pass
    try:
        G.add_edge(p10, c5, ("IS_PART_OF", qSIGNe(G, s=p10, t=c5, v=e7, e=p5)))
    except TypeError:
        pass
    G.remove_edge(v1, e2)
    try:
        G.remove_edge(v1, e5)
    except nx.exception.NetworkXError:
        pass
    G.remove_edge(v2, e3)
    try:
        G.remove_edge(v2, e6)
    except nx.exception.NetworkXError:
        pass
    G.remove_edge(e7, p1)
    G.remove_edge(e7, p2)
    try:
        G.remove_edge(e7, p3)
    except nx.exception.NetworkXError:
        pass

# TE
# This funciton performs pattern matching for TE transition
def patternMatching_TE(G, p_id):
    p10 = ("Polygon", p_id)
    es_of_p10 = list(G.predecessors(p10))
    e7 = es_of_p10[0]
    v1v2 = list(G.predecessors(e7))
    for v in v1v2:
        r_ve = list(G.get_edge_data(v, e7))[0][1]
        if r_ve == 1:
            v1 = v
        else:
            v2 = v
    for e_of_p10 in es_of_p10:
        vs_of_e = list(G.predecessors(e_of_p10))
        for v_of_e in vs_of_e:
            if v_of_e != v1 and v_of_e != v2:
                v3 = v_of_e
                break
    es_of_v2 = list(G.successors(v2))
    e7e8 = list(set(es_of_v2) & set(es_of_p10))
    for e in e7e8:
        if e != e7:
            e8 = e
            break
    for ep in es_of_p10:
        if ep != e7 and ep != e8:
            e9 = ep
            break
    c4 = list(G.successors(p10))[0]
    e1 = qVCEE_E(G, v=v1, c=c4, ea=e7, eb=e9) # e1
    e3 = qVEEE_E(G, v=v1, ea=e1, eb=e7, ec=e9) # e3
    e2 = qVCEE_E(G, v=v2, c=c4, ea=e7, eb=e8) # e2
    e4 = qVEEE_E(G, v=v2, ea=e2, eb=e7, ec=e8) # e4
    e5 = qVCEE_E(G, v=v3, c=c4, ea=e8, eb=e9) # e5
    e6 = qVEEE_E(G, v=v3, ea=e5, eb=e8, ec=e9) # e6
    p4 = qFIND(G, n1=e1, n2=e7) # p4
    p6 = qFIND(G, n1=e2, n2=e8) # p6
    p8 = qFIND(G, n1=e5, n2=e9) # p8
    p1 = qFIND(G, n1=e1, n2=e3) # p1
    p2 = qFIND(G, n1=e2, n2=e4) # p2
    p5 = qFIND(G, n1=e3, n2=e4) # p5
    p3 = qFIND(G, n1=e5, n2=e6) # p3
    p7 = qFIND(G, n1=e4, n2=e8) # p7
    p9 = qFIND(G, n1=e3, n2=e9) # p9
    c1 = qFIND(G, n1=p1, n2=p2) # c1
    c2 = qFIND(G, n1=p8, n2=p9) # c2
    c3 = qFIND(G, n1=p6, n2=p7) # c3
    c5 = qFIND(G, n1=p7, n2=p9) # c5
    return v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5

# This funciton performs graph transformation for TE transition
def graphTransformation_TE(G, v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5):  
    G.add_edge(v1, e2, list(G.get_edge_data(v2, e2))[0])
    G.add_edge(v1, e5, list(G.get_edge_data(v3, e5))[0])
    try:
        G.add_edge(v2, e3, list(G.get_edge_data(v1, e3))[0])
    except TypeError:
        pass
    try:
        G.add_edge(v2, e6, list(G.get_edge_data(v3, e6))[0])
    except TypeError:
        pass
    try:
        G.add_edge(e7, p1, ("IS_PART_OF", qSIGNe(G, s=e7, t=p1, v=v2, e=e3)))
    except TypeError:
        pass
    try:
        G.add_edge(e7, p2, ("IS_PART_OF", qSIGNe(G, s=e7, t=p2, v=v1, e=e2)))
    except TypeError:
        pass
    try:
        G.add_edge(e7, p3, ("IS_PART_OF", qSIGNe(G, s=e7, t=p3, v=v2, e=e6)))
    except TypeError:
        pass
    try:
        G.remove_edge(v1, e3)
    except nx.exception.NetworkXError:
        pass
    G.remove_edge(v1, e9)
    G.remove_edge(v2, e2)
    G.remove_edge(v2, e8)
    G.remove_edge(v3, e5)
    try:
        G.remove_edge(v3, e6)
    except nx.exception.NetworkXError:
        pass
    G.remove_edge(v3, e8)
    G.remove_edge(v3, e9)
    G.remove_edge(e7, p4)
    try:
        G.remove_edge(e7, p5)
    except nx.exception.NetworkXError:
        pass
    G.remove_edge(e7, p10)
    G.remove_edge(e8, p6)
    try:
        G.remove_edge(e8, p7)
    except nx.exception.NetworkXError:
        pass
    G.remove_edge(e8, p10)
    G.remove_edge(e9, p8)
    try:
        G.remove_edge(e9, p9)
    except nx.exception.NetworkXError:
        pass
    G.remove_edge(e9, p10)
    G.remove_edge(p10, c4)
    try:
        G.remove_edge(p10, c5)
    except nx.exception.NetworkXError:
        pass
    G.add_edge("Free", v3, "CONTAINS")
    G.add_edge("Free", e8, "CONTAINS")
    G.add_edge("Free", e9, "CONTAINS")
    G.add_edge("Free", p10, "CONTAINS")

