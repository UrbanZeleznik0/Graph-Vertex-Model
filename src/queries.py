import networkx as nx
####################################################################################################
#################################### MISCELLANEOUS QUERIES #########################################
####################################################################################################
# This function counts common polygons between cells c4 and c5
def count_shared_polygons_c4c5(e_id, G):
    ps_of_e = list(G.successors(("Edge", e_id)))
    cs_of_ps = []
    for p_of_e in ps_of_e:
        cs_of_ps += list(G.successors(p_of_e))
    cs_e_distinct = list(set(cs_of_ps))
    cs_all = []
    vs_of_e = list(G.predecessors(("Edge", e_id)))
    for v_of_e in vs_of_e:
        es_of_v = list(G.successors(v_of_e))
        for e_of_v in es_of_v:
            ps_of_ev = list(G.successors(e_of_v))
            for p_of_ev in ps_of_ev:
                cs_of_pev = list(G.successors(p_of_ev))
                cs_all += cs_of_pev
    cs_all_distinct = list(set(cs_all))
    c4c5 = list(set(cs_all_distinct) - set(cs_e_distinct))
    try:
        c4 = c4c5[0]
        c5 = c4c5[1]
        ps_of_c4 = list(G.predecessors(c4))
        ps_of_c5 = list(G.predecessors(c5))
        num_shared_polygons_c4c5 = len(set(ps_of_c4) & set(ps_of_c5))
        return num_shared_polygons_c4c5
    except IndexError:
        return 0

def impossible_ET(G, p4, p5, p6, p7, p8, p9):
    es_of_p4 = list(G.predecessors(p4))
    es_of_p5 = list(G.predecessors(p5))
    num_e45 = len(set(es_of_p4) & set(es_of_p5))
    try:
        es_of_p6 = list(G.predecessors(p6))
        es_of_p7 = list(G.predecessors(p7))
        num_e67 = len(set(es_of_p6) & set(es_of_p7))
    except nx.exception.NetworkXError:
        num_e67 = 0
    try:
        es_of_p8 = list(G.predecessors(p8))
        es_of_p9 = list(G.predecessors(p9))
        num_e89 = len(set(es_of_p8) & set(es_of_p9))
    except nx.exception.NetworkXError:
        num_e89 = 0
    return [num_e45, num_e67, num_e89]

def impossible_TE(G, p1, p2, p3):
    try:
        es_of_p1 = list(G.predecessors(p1))
        es_of_p2 = list(G.predecessors(p2))
        num_e12 = len(set(es_of_p1) & set(es_of_p2))
    except nx.exception.NetworkXError:
        num_e12 = 0
    try:
        es_of_p1 = list(G.predecessors(p1))
        es_of_p3 = list(G.predecessors(p3))
        num_e13 = len(set(es_of_p1) & set(es_of_p3))
    except nx.exception.NetworkXError:
        num_e13 = 0
    try:
        es_of_p2 = list(G.predecessors(p2))
        es_of_p3 = list(G.predecessors(p3))
        num_e23 = len(set(es_of_p2) & set(es_of_p3))
    except nx.exception.NetworkXError:
        num_e23 = 0
    return [num_e12, num_e13, num_e23]

def hos_of_lo_hho(G, lo, hho):
    hos1 = list(G.successors(lo))
    hos2 = [ho2 for ho2 in list(G.predecessors(hho)) if ho2[0] != "Div"]
    hos_common = list(set(hos1) & set(hos2))
    return hos_common

def qSIGNe(G, s, t, v, e):
    r_vs = list(G.get_edge_data(v, s))[0][1]
    r_ve = list(G.get_edge_data(v, e))[0][1]
    r_et = list(G.get_edge_data(e, t))[0][1]
    r = -r_vs*r_ve*r_et
    return r
