####################################################################################################
############################ USEFUL QUERY TEMPLATES ################################################
####################################################################################################
def qFREE(G, label):
    free_nodes = list(G.successors("Free"))
    for free_node in free_nodes:
        if free_node[0] == label:
            G.remove_edge("Free", free_node)
            break
    return free_node

def qFREEd(G):
    free_nodes = list(G.successors("Free"))
    for free_node in free_nodes:
        if free_node[0] == "Div":
            G.remove_edge("Free", free_node)

def qSIGNe(G, s, t, v, e):
    r_vs = list(G.get_edge_data(v, s))[0][1]
    r_ve = list(G.get_edge_data(v, e))[0][1]
    r_et = list(G.get_edge_data(e, t))[0][1]
    r = -r_vs*r_ve*r_et
    return r

def next_edge_vertex(G, v, e, p):
    es_of_v = list(G.successors(v))
    es_of_p = [e for e in list(G.predecessors(p)) if e[0] != "Div"]
    es_common = list(set(es_of_v) & set(es_of_p))
    if len(es_common) < 2:
        print("mistake")
        return None, None
    for e_common in es_common:
        if e_common != e:
            e_next = e_common
            break
    vs_of_e = [v for v in list(G.predecessors(e_next)) if v[0] != "Div"]
    for v_of_e in vs_of_e:
        if v_of_e != v:
            v_next = v_of_e
            break
    return e_next, v_next

def side_e_v2(G, div, e_start, v1, v_first, v_second, p5):
    e = e_start
    v = v1
    num_side_edges = 0
    side_edges = []
    side_vertices = []
    wrong_div = False
    while v != v_first and v != v_second: # while we don't reach the vertices of cut 'e2'
        e_next, v_next = next_edge_vertex(G, v, e, p5) # find the next edge-vertex pair
        if e_next is None:
            wrong_div = True
            break
        e = e_next
        side_edges.append(e_next)
        G.add_edge(div, e_next, "e{}".format(6 + num_side_edges))
        v = v_next
        if v != v_first and v != v_second:
            side_vertices.append(v_next)
            G.add_edge(div, v_next, "v{}".format(5 + num_side_edges))
        num_side_edges += 1
    v2 = v # the last one is 'v2'
    G.add_edge(div, v2, "v2")
    return num_side_edges, side_edges, side_vertices, v2, wrong_div

def ho_of_lo12(G, lo1, lo2):
    hos1 = list(G.successors(lo1))
    hos2 = list(G.successors(lo2))
    ho_common = list(set(hos1) & set(hos2))[0]
    return ho_common

def lo_of_ho12(G, ho1, ho2):
    los1 = [lo1 for lo1 in list(G.predecessors(ho1)) if lo1[0] != "Div"]
    los2 = [lo2 for lo2 in list(G.predecessors(ho2)) if lo2[0] != "Div"]
    lo_common = list(set(los1) & set(los2))[0]
    return lo_common

def ho_of_lo_not_hho(G, lo, hho):
    hos = list(G.successors(lo))
    wrong_hos = [wh for wh in list(G.predecessors(hho)) if wh[0] != "Div"]
    for ho in hos:
        if ho not in wrong_hos:
            n = ho
            break
    return n

def exist_node(G, div, name):
    nodes_div = list(G.successors(div))
    for node_div in nodes_div:
        if list(G.get_edge_data(div, node_div))[0] == name:
            n = node_div
            break
    return n

def prev_div(G, div_id):
    v2 = None
    v4 = None
    nodes_div = list(G.successors(("Div", div_id)))
    for node_div in nodes_div:
        name = list(G.get_edge_data(("Div", div_id), node_div))[0]
        if name == "e2":
            e2 = node_div
        if name == "v4":
            v4 = node_div
        if name == "e4":
            e4 = node_div
        if name == "v2":
            v2 = node_div
        if name == "v1":
            v1 = node_div
    if v2 is None:
        v2 = v1
    if v4 is None:
        return None, None, None, None
    return e2, v4, e4, v2

def first_div(G):
    nodes_div = list(G.successors(("Div", 0)))
    for node_div in nodes_div:
        name = list(G.get_edge_data(("Div", 0), node_div))[0]
        if name == "e1":
            e1 = node_div
        if name == "v3":
            v3 = node_div
        if name == "e3":
            e3 = node_div
        if name == "v1":
            v1 = node_div
    return e1, v3, e3, v1

def prev_div_e2_p5(G, div_id):
    nodes_div = list(G.successors(("Div", div_id)))
    for node_div in nodes_div:
        name = list(G.get_edge_data(("Div", div_id), node_div))[0]
        if name == "e2":
            e2 = node_div
        if name == "p5":
            p5 = node_div
    return e2, p5

def other_ho_of_lo(G, lo, ho, hho):
    hos1 = list(G.successors(lo))
    hos2 = [h for h in list(G.predecessors(hho)) if h[0] != "Div"]
    hos_common = list(set(hos1) & set(hos2))
    for ho_common in hos_common:
        if ho_common != ho:
            n = ho_common
            break
    return n

def nodes_of_divs(G, num_cut_ps, target_name):
    node_list = []
    for i in range(num_cut_ps):
        nodes_div = list(G.successors(("Div", i)))
        for node_div in nodes_div:
            name = list(G.get_edge_data(("Div", i), node_div))[0]
            if name == target_name:
                node_list.append(node_div)
    return node_list

def next_lower(G, p, c, p5_list):
    neighbor_ps = []
    es_of_p = [e for e in list(G.predecessors(p)) if e[0] != "Div"]
    ps_of_c = [p for p in list(G.predecessors(c)) if p[0] != "Div"]
    for p_of_c in ps_of_c:
        if p_of_c not in p5_list:
            es_of_pc = [e for e in list(G.predecessors(p_of_c)) if e[0] != "Div"]
            es_common = list(set(es_of_p) & set(es_of_pc))
            if len(es_common) == 1:
                neighbor_ps.append(p_of_c)
    return neighbor_ps

def d_of_p(G, num_cut_ps, p):
    d_list = []
    for i in range(num_cut_ps):
        nodes_div = list(G.successors(("Div", i)))
        for node_div in nodes_div:
            if node_div == p:
                d_list.append(("Div", i))
                break
    return d_list

def neighbor_ps(G, p, c2):
    neighbors = []
    ps_of_c2 = [p for p in list(G.predecessors(c2)) if p[0] != "Div"]
    es_of_p = [e for e in list(G.predecessors(p)) if e[0] != "Div"]
    ps_of_es = []
    for e in es_of_p:
        ps_of_es += list(G.successors(e))
    ps_common = list(set(ps_of_es) & set(ps_of_c2))
    for p_common in ps_common:
        if p_common != p:
            neighbors.append(p_common)
    return neighbors


####################################################################################################
#######################PATTERN-MATCHING & GRAPH-TRANSFORMATION #####################################
####################################################################################################

# This funciton performs pattern matching for division
def patternMatching_D(G, div, dividing_cell, cut_e_first, cut_e_second, num_cut_ps):
    p5 = ho_of_lo12(G, cut_e_first, cut_e_second)
    G.add_edge(div, p5, "p5")
    c1 = dividing_cell
    G.add_edge(div, c1, "c1")
    
    cs_of_p5 = list(G.successors(p5))
    for c_of_p5 in cs_of_p5:
        if c_of_p5 != c1:
            c3 = c_of_p5
            break
    try:
        G.add_edge(div, c3, "c3")
    except UnboundLocalError:
        c3 = None
    
    e5 = qFREE(G, "Edge")
    G.add_edge(div, e5, "e5")
    p4 = qFREE(G, "Polygon")
    G.add_edge(div, p4, "p4")

    if div[1] == 0:
        p3 = qFREE(G, "Polygon")
        G.add_edge(div, p3, "p3")
        c2 = qFREE(G, "Cell")
        G.add_edge(div, c2, "c2")
        e1 = cut_e_first
        G.add_edge(div, e1, "e1")
        e2 = cut_e_second
        G.add_edge(div, e2, "e2")
        
        v1 = [v for v in list(G.predecessors(e1)) if v[0] != "Div"][0]
        G.add_edge(div, v1, "v1")

        vv_of_e2 = [v for v in list(G.predecessors(e2)) if v[0] != "Div"]
        v_first = vv_of_e2[0]
        v_second = vv_of_e2[1]

        num_side_edges, side_edges, side_vertices, v2, wrong_div = side_e_v2(G, div, e1, v1, v_first, v_second, p5)
        if wrong_div:
            return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None

        v3 = qFREE(G, "Vertex")
        G.add_edge(div, v3, "v3")
        e3 = qFREE(G, "Edge")
        G.add_edge(div, e3, "e3")
        v4 = qFREE(G, "Vertex")
        G.add_edge(div, v4, "v4")
        e4 = qFREE(G, "Edge")
        G.add_edge(div, e4, "e4")

    else:
        p3 = exist_node(G, ("Div", 0), "p3")
        G.add_edge(div, p3, "p3")
        c2 = exist_node(G, ("Div", 0), "c2")
        G.add_edge(div, c2, "c2")

        prev_e2, prev_v4, prev_e4, prev_v2 = prev_div(G, div[1]-1)
        if prev_e2 is None:
            return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None

        e1 = prev_e2
        G.add_edge(div, e1, "e1")
        v3 = prev_v4
        G.add_edge(div, v3, "v3")
        e3 = prev_e4
        G.add_edge(div, e3, "e3")
        G.add_edge(e3, p5)
        v1 = prev_v2
        G.add_edge(div, v1, "v1")
        
        if e1 == cut_e_first:
            e2 = cut_e_second
        elif e1 == cut_e_second:
            e2 = cut_e_first
        G.add_edge(div, e2, "e2")

        vv_of_e2 = [v for v in list(G.predecessors(e2)) if v[0] != "Div"]
        v_first = vv_of_e2[0]
        v_second = vv_of_e2[1]
        
        if div[1] == num_cut_ps-1:
            first_e1, first_v3, first_e3, first_v1 = first_div(G) # 'e1', 'v3', 'e3', 'v1' of the first polygon division
            num_side_edges, side_edges, side_vertices, v2, wrong_div = side_e_v2(G, div, e3, v1, first_v1, first_v1, p5) # v2 = first_v1
            if wrong_div:
                return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None
            v4 = first_v3
            G.add_edge(div, v4, "v4")
            e4 = first_e3
            G.add_edge(div, e4, "e4")
            G.add_edge(e4, p5)
        else:
            num_side_edges, side_edges, side_vertices, v2, wrong_div = side_e_v2(G, div, e3, v1, v_first, v_second, p5)
            if wrong_div:
                return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None
            v4 = qFREE(G, "Vertex")
            G.add_edge(div, v4, "v4")
            e4 = qFREE(G, "Edge")
            G.add_edge(div, e4, "e4")

    p1 = ho_of_lo_not_hho(G, e1, c1)
    G.add_edge(div, p1, "p1")
    p2 = ho_of_lo_not_hho(G, e2, c1)
    G.add_edge(div, p2, "p2")

    return num_side_edges, side_edges, side_vertices, v1, v2, v3, v4, e1, e2, e3, e4, e5, p1, p2, p3, p4, p5, c1, c2, c3


# This funciton performs graph transformation for division
def graph_transformations_D(G, div, num_cut_ps, num_side_edges, side_edges, side_vertices, v1, v2, v3, v4, e1, e2, e3, e4, e5, p1, p2, p3, p4, p5, c1, c2, c3):
    if div[1] != num_cut_ps-1:
        G.add_edge(v2, e4, ("IS_PART_OF", 1))
        G.add_edge(v4, e2, list(G.get_edge_data(v2, e2))[0])
        G.add_edge(v4, e4, ("IS_PART_OF", -1))
        G.add_edge(e4, p2, ("IS_PART_OF", qSIGNe(G, s=e4, t=p2, v=v4, e=e2)))
        G.remove_edge(v2, e2)
        if div[1] == 0:
            G.add_edge(v1, e3, ("IS_PART_OF", 1))
            G.add_edge(v3, e1, list(G.get_edge_data(v1, e1))[0])
            G.add_edge(v3, e3, ("IS_PART_OF", -1))
            G.add_edge(e3, p1, ("IS_PART_OF", qSIGNe(G, s=e3, t=p1, v=v3, e=e1)))
            G.remove_edge(v1, e1)
        else:
            G.remove_edge(e3, p5)
    else:
        G.remove_edge(e3, p5)
        G.remove_edge(e4, p5)
    
    G.add_edge(v3, e5, ("IS_PART_OF", 1))
    G.add_edge(v4, e5, ("IS_PART_OF", -1))
    G.add_edge(e3, p4, ("IS_PART_OF", 1))
    G.add_edge(e5, p3, ("IS_PART_OF", 1))
    G.add_edge(e5, p4, ("IS_PART_OF", qSIGNe(G, s=e5, t=p4, v=v3, e=e3)))
    G.add_edge(e4, p4, ("IS_PART_OF", qSIGNe(G, s=e4, t=p4, v=v4, e=e5)))
    G.add_edge(e5, p5, ("IS_PART_OF", qSIGNe(G, s=e5, t=p5, v=v3, e=e1)))
    G.add_edge(p4, c1, ("IS_PART_OF", qSIGNe(G, s=p4, t=c1, v=e5, e=p5)))
    try:
        G.add_edge(p4, c3, ("IS_PART_OF", qSIGNe(G, s=p4, t=c3, v=e5, e=p5)))
    except TypeError:
        pass
    G.add_edge(c2, ("Tissue", 0))

    for i in range(num_side_edges):
        if i == 0:
            G.add_edge(side_edges[0], p4, ("IS_PART_OF", qSIGNe(G, s=side_edges[0], t=p4, v=v1, e=e3)))
        else:
            G.add_edge(side_edges[i], p4, ("IS_PART_OF", qSIGNe(G, s=side_edges[i], t=p4, v=side_vertices[i-1], e=side_edges[i-1])))
        G.remove_edge(side_edges[i], p5)
        
    if div[1] == 0:
        G.add_edge(p3, c1, ("IS_PART_OF", qSIGNe(G, s=p3, t=c1, v=e5, e=p5)))
        G.add_edge(p3, c2, ("IS_PART_OF", -list(G.get_edge_data(p3, c1))[0][1]))


def PM_GT_D(G, dividing_cell, cut_ps, num_cut_ps, cut_es):
    qFREEd(G)
    for div_id in range(len(cut_ps)):
        if div_id == 0:
            cut_p = cut_ps[0]
        else:
            prev_e2, prev_p5 = prev_div_e2_p5(G, div_id-1) # 'e2', 'p5' of the previous polygon division
            cut_p = other_ho_of_lo(G, prev_e2, prev_p5, dividing_cell)
        es_of_cut_p = [e for e in list(G.predecessors(cut_p)) if e[0] != "Div"] # all edges of cut polygon
        cut_e1e2 = [] # cut edges of a cut polygon
        for cut_e in cut_es: # loop over all cut edges
            if cut_e in es_of_cut_p: # if the cut edge belongs to this polygon
                cut_e1e2.append(cut_e) # cut edge of this cut polygon
        cut_e_first = cut_e1e2[0]
        cut_e_second = cut_e1e2[1]
        num_side_edges, side_edges, side_vertices, v1, v2, v3, v4, e1, e2, e3, e4, e5, p1, p2, p3, p4, p5, c1, c2, c3 = patternMatching_D(G, ("Div", div_id), dividing_cell, cut_e_first, cut_e_second, num_cut_ps)
        if num_side_edges is None:
            return False
        graph_transformations_D(G, ("Div", div_id), num_cut_ps, num_side_edges, side_edges, side_vertices, v1, v2, v3, v4, e1, e2, e3, e4, e5, p1, p2, p3, p4, p5, c1, c2, c3)
    return True


def lower_polygons(G, dividing_cell, num_cut_ps):
    p4_list = nodes_of_divs(G, num_cut_ps, "p4") # list of 'p4' polygons
    p5_list = nodes_of_divs(G, num_cut_ps, "p5") # list of 'p5' polygons
    lower_ps = [] # list of all lower polygons
    lower_ps_full = False # list of all lower polygons is not yet full
    p4_of_d0 = exist_node(G, ("Div", 0), "p4") # 'p4' of the 0th (this is arbitrary) polygon division node
    p3_of_d0 = exist_node(G, ("Div", 0), "p3") # 'p3' of the 0th (this is arbitrary) polygon division node
    lower_ps.append(p4_of_d0) # append 'p4' to the list of all lower polygons
    overlap_p4 = []
    overlap_p4.append(p4_of_d0)
    i = 0
    while not lower_ps_full: # while the list of lower polygons is not full
        next_lower_ps = next_lower(G, lower_ps[i], dividing_cell, p5_list) # neighboring lower polygons of chosen lower polygon
        for next_lower_p in next_lower_ps: # for all neighboring lower polygons of chosen lower polygon
            if (next_lower_p not in lower_ps) and (next_lower_p != p3_of_d0): # if the polygon was not yet recognized as lower
                lower_ps.append(next_lower_p) # append the polygon to the lower polygons (mark it as lower)
            if next_lower_p in p4_list: # if the polygon is also 'p1' or 'p2'
                if next_lower_p not in overlap_p4: # if the polygon was not yet recognized as 'p1' or 'p2'
                    overlap_p4.append(next_lower_p) # append the polygon to the overlap
        if len(overlap_p4) == len(p4_list) and i >= len(lower_ps)-1:
            lower_ps_full = True # list of all lower polygons is full (stop iterating)
        i += 1
    num_extra_ps = len(lower_ps) - len(p4_list)
    return lower_ps, num_extra_ps


def lower_ps_transformations(G, dividing_cell, lower_ps, num_cut_ps):
    qFREEd(G)
    num_extra_done = 0
    done_ps = []
    for lower_p in lower_ps:
        d_list_p = d_of_p(G, num_cut_ps, lower_p) # polygon 'div' node of the given lower polygon
        if d_list_p: # if the polygon already has a 'div' node (in the case of 'p4' polygons)
            div = d_list_p[0] # name it 'div' as always
            p4 = exist_node(G, div, "p4")
            c1 = exist_node(G, div, "c1")
            c2 = exist_node(G, div, "c2")
            e5 = exist_node(G, div, "e5")
            p3 = exist_node(G, div, "p3")
            G.add_edge(p4, c2, ("IS_PART_OF", qSIGNe(G, s=p4, t=c2, v=e5, e=p3)))
            G.remove_edge(p4, c1)
        else: # if the polygon division node doesn't exist
            div = ("Div", num_cut_ps + num_extra_done)
            G.add_edge(div, lower_p, "p") # connects the lower polygon to the new division node
            
            c1 = dividing_cell
            G.add_edge(div, c1, "c1") # connects 'c1' to the new division node

            c2 = exist_node(G, ("Div", 0), "c2") # use the existing 'c2' node
            G.add_edge(div, c2, "c2") # connects 'c2' to the new division node

            G.add_edge(lower_p, c2)
            G.remove_edge(lower_p, c1)
            neigh_ps = neighbor_ps(G, lower_p, c2)
            for neigh_p in neigh_ps:
                if neigh_p in done_ps:
                    edg = lo_of_ho12(G, lower_p, neigh_p)
                    G.remove_edge(lower_p, c2)
                    G.add_edge(lower_p, c2, ("IS_PART_OF", qSIGNe(G, s=lower_p, t=c2, v=edg, e=neigh_p)))
                    break
            num_extra_done += 1
        done_ps.append(lower_p)

