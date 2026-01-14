import math
import random
import os
import networkx as nx
import queries as Q
import queriesT1 as T1
import queriesD as D

class database:
    def __init__(self, G, t):
        # essential variables
        self.file=''
        self.G=G
        # functions
        self.setup_DB(G, t)
    
    # This function sets up the knowledge-graph database
    def setup_DB(self, G, t):
        t.generate_DB_setup_script(G)
    
    # This function creates a new free 'Vertex' node
    def add_vertex(self, G, t):
        t.Nv += 1
        G.add_node(("Vertex", t.Nv-1))
        G.add_edge("Free", ("Vertex", t.Nv-1), "CONTAINS")
        vertex(t.Nv-1, [0,0,0], t)
    
    # This function creates a new free 'Edge' node
    def add_edg(self, G, t):
        t.Ne += 1
        G.add_node(("Edge", t.Ne-1))
        G.add_edge("Free", ("Edge", t.Ne-1), "CONTAINS")
        edge(t.Ne-1, [], t)
    
    # This function creates a new free 'Polygon' node
    def add_polygon(self, G, t):
        t.Np += 1
        G.add_node(("Polygon", t.Np-1))
        G.add_edge("Free", ("Polygon", t.Np-1), "CONTAINS")
        polygon(t.Np-1, [], [], t)
    
    # This function creates a new free 'Cell' node
    def add_cell(self, G, t):
        t.Nc += 1
        G.add_node(("Cell", t.Nc-1))
        G.add_edge("Free", ("Cell", t.Nc-1), "CONTAINS")
        cell(t.Nc-1, [], [], t)
    
    # This function creates a new free 'Div' node
    def add_div(self, G, t):
        t.Nd += 1
        G.add_node(("Div", t.Nd-1))
        G.add_edge("Free", ("Div", t.Nd-1), "CONTAINS")
    
    # This function creates free nodes in case there are not enough for an ET transition
    def add_nodes_ET(self, G, t):
        free_nodes = list(G.successors("Free"))
        free_vertex_nodes = 0
        free_edge_nodes = 0
        free_polygon_nodes = 0
        for free_node in free_nodes:
            if free_node[0] == "Vertex":
                free_vertex_nodes += 1
            if free_node[0] == "Edge":
                free_edge_nodes += 1
            if free_node[0] == "Polygon":
                free_polygon_nodes += 1
        if free_vertex_nodes<1:
            self.add_vertex(G, t)
        # edge nodes
        if free_edge_nodes<2:
            for i in range(2-free_edge_nodes):
                self.add_edg(G, t)
        # polygon nodes
        if free_polygon_nodes<1:
            self.add_polygon(G, t)
    
    # This function creates free nodes in case there are not enough for a division
    def add_nodes_D(self, G, t, num_cut_ps):
        num_new_v = num_cut_ps
        num_new_e = 2*num_cut_ps
        num_new_p = num_cut_ps + 1
        num_new_d = num_cut_ps
        free_nodes = list(G.successors("Free"))
        free_vertex_nodes = 0
        free_edge_nodes = 0
        free_polygon_nodes = 0
        free_cell_nodes = 0
        free_div_nodes = 0
        for free_node in free_nodes:
            if free_node[0] == "Vertex":
                free_vertex_nodes += 1
            if free_node[0] == "Edge":
                free_edge_nodes += 1
            if free_node[0] == "Polygon":
                free_polygon_nodes += 1
            if free_node[0] == "Cell":
                free_cell_nodes += 1
            if free_node[0] == "Div":
                free_div_nodes += 1
        # vertex nodes
        nfv = free_vertex_nodes
        if nfv < num_new_v:
            for i in range(num_new_v - nfv):
                self.add_vertex(G, t)
        # edge nodes
        nfe = free_edge_nodes
        if nfe < num_new_e:
            for i in range(num_new_e - nfe):
                self.add_edg(G, t)
        # polygon nodes
        nfp = free_polygon_nodes
        if nfp < num_new_p:
            for i in range(num_new_p - nfp):
                self.add_polygon(G, t)
        # cell node
        if free_cell_nodes<1:
            self.add_cell(G, t)
        # div nodes
        nfd = free_div_nodes
        if nfd < num_new_d:
            for i in range(num_new_d - nfd):
                self.add_div(G, t)
    
    # Add additional graph nodes for division
    def add_extra_nodes_D(self, G, t, num_extra_ps):
        free_nodes = list(G.successors("Free"))
        # div nodes
        free_div_nodes = 0
        for free_node in free_nodes:
            if free_node[0] == "Div":
                free_div_nodes += 1
        nfd = free_div_nodes
        if nfd < num_extra_ps:
            for i in range(num_extra_ps - nfd):
                self.add_div(G, t)


class vertex:
    def __init__(self, id, r, t):
        # essential variables
        self.id=id
        self.r=r
        self.exist=1
        # reset vertex properties
        self.reset_prop(t)
        # add it to a list of vertices
        t.vertices.append(self)
    
    # This function resets vertex properties
    def reset_prop(self,t):
        self.F=[0,0,0]


class edge:
    def __init__(self, id, e, t):
        # essential variables
        self.id=id
        self.e=e
        self.exist=1
        self.sigma=t.sig
        # reset edge properties
        self.reset_prop(t)
        # add it to a list of edges
        t.edges.append(self)
    
    # This function resets edge properties
    def reset_prop(self,t):
        self.length=0
        self.e_dl=0
        self.center=[0,0,0]
        self.clock=0
        self.lineTension=t.dg0+self.sigma*random.normalvariate(0,1)
        if self.lineTension<0:
            self.lineTension=0
    
    # This function gives an edge's vertices in a periodic boundary system
    def edge_vertices_pbc(self,t):
        # vertex ids
        v1_id=self.e[0]
        v2_id=self.e[1]
        # pbc correct
        v1=t.torus_corrected_position(t.vertices[v1_id].r,t.vertices[v2_id].r) #v1
        v2=[t.vertices[v2_id].r[0],t.vertices[v2_id].r[1],t.vertices[v2_id].r[2]] #v2
        # return
        return [v1_id,v2_id,v1,v2]
    
    # This function calculates edge length
    def edge_length(self,t):
        # vertices
        [v1_id,v2_id,v1,v2]=self.edge_vertices_pbc(t)
        # length calculation
        dr=[0,0,0]
        for i in range(3):
            dr[i]=v2[i]-v1[i]
        return math.sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])   
    
    # This function calculates forces due to line tension
    def force_lineTension(self,t):
        # vertices
        [v1_id,v2_id,v1,v2]=self.edge_vertices_pbc(t)
        # force calculation
        dr=[0,0,0]
        for i in range(3):
            dr[i]=v2[i]-v1[i]
        c0=-self.lineTension/(self.length+1e-8)
        for i in range(3):
            t.vertices[v1_id].F[i] += -c0*dr[i]
            t.vertices[v2_id].F[i] +=  c0*dr[i]
        return self.lineTension*self.length 
    
    # This function calculates edge midpoint
    def edge_center(self,t):
        # vertices
        [v1_id,v2_id,v1,v2]=self.edge_vertices_pbc(t)
        # midpoint calculation
        v_mid=[0,0,0]
        for i in range(3):
            v_mid[i]=0.5*(v1[i]+v2[i])
        # return pbc-corrected position
        return t.torus_on_position(v_mid)
    
    # This function corrects vertex positions upon an ET transition
    def correct_vertex_positions_ET(self, t, v1, v2, v3, p1, p2, p3):
        v1ID, v2ID, v3ID, p1ID, p2ID = v1[1], v2[1], v3[1], p1[1], p2[1]
        try:
            p3ID = p3[1]
        except TypeError:
            p3ID = None
        v1v2v3=[v1ID, v2ID, v3ID]
        p1p2p3=[p1ID, p2ID, p3ID]
        # correct vertex positions
        v_corr=[]
        if p3ID is not None:
            for i in range(3):
                # storing new coordinates in a 2D array
                v1_id=v1v2v3[0]
                v2_id=v1v2v3[1]
                # move v1 to mid and v2 to center
                t.vertices[v1_id].r=[self.center[0],self.center[1],self.center[2]]
                t.vertices[v2_id].r=[t.polygons[p1p2p3[i]].center[0],t.polygons[p1p2p3[i]].center[1],t.polygons[p1p2p3[i]].center[2]]
                v2corr=t.torus_corrected_position(t.vertices[v2_id].r, t.vertices[v1_id].r)
                dr=[0,0,0]
                for j in range(3):
                    # a vector from mid point (v1 now) to the center of the neighbour polygon (v2 now)
                    dr[j]=v2corr[j]-t.vertices[v1_id].r[j]
                len_temp=math.sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])
                # New coordinates of v2
                for j in range(3):
                    t.vertices[v2_id].r[j]=t.vertices[v1_id].r[j]+t.new_len*(dr[j]/len_temp)
                # taking care of periodic boundary condition
                v2=t.torus_on_position(t.vertices[v2_id].r)
                # return
                v_corr.append(v2)
        else:
            vector_center_mid = []
            for i in range(2):
                # storing new coordinates in a 2D array
                v1_id=v1v2v3[0]
                v2_id=v1v2v3[1]
                # move v1 to mid and v2 to center
                t.vertices[v1_id].r=[self.center[0],self.center[1],self.center[2]]
                t.vertices[v2_id].r=[t.polygons[p1p2p3[i]].center[0],t.polygons[p1p2p3[i]].center[1],t.polygons[p1p2p3[i]].center[2]]
                v2corr=t.torus_corrected_position(t.vertices[v2_id].r, t.vertices[v1_id].r)
                dr=[0,0,0]
                for j in range(3):
                    # a vector from mid point (v1 now) to the center of the neighbour polygon (v2 now)
                    dr[j]=v2corr[j]-t.vertices[v1_id].r[j]
                len_temp=math.sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])
                # a vector from the center of the neighbour polygon (v2 now) to mid point (v1 now)
                vector_center_mid.append([-dr[0],-dr[1],-dr[2]])
                # New coordinates of v2
                for j in range(3):
                    t.vertices[v2_id].r[j]=t.vertices[v1_id].r[j]+t.new_len*(dr[j]/len_temp)
                # taking care of periodic boundary condition
                v2=t.torus_on_position(t.vertices[v2_id].r)
                # return
                v_corr.append(v2)
            v1_id=v1v2v3[0]
            v2_id=v1v2v3[1]
            t.vertices[v1_id].r=[self.center[0],self.center[1],self.center[2]]
            vector_mid_v3 = [vector_center_mid[0][0]+vector_center_mid[1][0], vector_center_mid[0][1]+vector_center_mid[1][1], vector_center_mid[0][2]+vector_center_mid[1][2]]
            # a vector from mid point (v1 now) to v3 (v2 now)
            dr=vector_mid_v3
            len_temp=math.sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])
            # New coordinates of v2
            for j in range(3):
                t.vertices[v2_id].r[j]=t.vertices[v1_id].r[j]+t.new_len*(dr[j]/len_temp)
            # taking care of periodic boundary condition
            v2=t.torus_on_position(t.vertices[v2_id].r)
            # return
            v_corr.append(v2)
        # assigning new cordinates of v1, v2, v3
        for i in range(3):
            t.vertices[v1v2v3[i]].r=v_corr[i]
    
    # This function performs ET transition
    def ET_transition(self, G, t, DB):
        # add new nodes if needed
        DB.add_nodes_ET(G, t)
        # pattern matching
        v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5 = T1.patternMatching_ET(G, self.id)
        impossible_ET = Q.impossible_ET(G, p4, p5, p6, p7, p8, p9)
        if impossible_ET[0] == 0 and impossible_ET[1] == 0 and impossible_ET[2] == 0:
            # graph transformations
            T1.graphTransformation_ET(G, v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5)
            # correct vertex positions
            self.correct_vertex_positions_ET(t, v1, v2, v3, p1, p2, p3)
            # updating tissue
            t.update_ET(G, v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5)
            # remove nonexistent entities
            t.remove_nonexistent(G)
            # counter
            t.ET_count+=1
        else:
            print('impossible ET')
            # remove nonexistent entities
            t.remove_nonexistent(G)

class polygon:
    def __init__(self, id:int, loe:int, loo:int, t):
        # essential variables
        self.id=id
        self.loe=loe #list of edges
        self.loo=loo #list of edge orientations
        self.exist=1
        self.g0=1 # surface tension
        self.boundary_polygon = False
        # reset polygon properties
        self.reset_lov(t)
        self.reset_prop(t)
        # add it to list of polygons
        t.polygons.append(self)  
    
    # This function resets polygon properties
    def reset_prop(self,t):
        self.center=[0,0,0]
        self.area=0
        self.p_dA=0
    
    # This function resets list of vertices
    def reset_lov(self,t):
        self.lov=[]
        # oriented list of vertices
        if len(self.loe)>0:
            if self.loo[0]==1:
                vi=t.edges[self.loe[0]].e[0]
                vf=t.edges[self.loe[0]].e[1]
            else:
                vi=t.edges[self.loe[0]].e[1]
                vf=t.edges[self.loe[0]].e[0]
            self.lov.append(vi)
            vprev=vf
            while vprev!=vi:
                for i in range(1,len(self.loe)):
                    if self.loo[i]==1:
                        v1=t.edges[self.loe[i]].e[0]
                        v2=t.edges[self.loe[i]].e[1]
                    else:
                        v1=t.edges[self.loe[i]].e[1]
                        v2=t.edges[self.loe[i]].e[0]
                    if v1==vprev:
                        self.lov.append(v1)
                        vprev=v2
                        break
    
    # This function calculates polygon center
    def polygon_center(self,t):
        self.reset_lov(t)
        # reference vertex
        vref=[t.vertices[self.lov[0]].r[0],t.vertices[self.lov[0]].r[1],t.vertices[self.lov[0]].r[2]]
        # calculate polygon center
        cent=[vref[0],vref[1],vref[2]]
        for i in range(1,len(self.lov)):
            v=t.torus_corrected_position(t.vertices[self.lov[i]].r,vref)
            for j in range(3):
                cent[j]+=v[j]
        return t.torus_on_position([cent[0]/(len(self.lov)*1.0),cent[1]/(len(self.lov)*1.0),cent[2]/(len(self.lov)*1.0)])
    
    # This function gives a triangle's vertices for a system with periodic boundaries
    def triangle_vertices_pbc(self,j,vref,t):
        # vertex ids
        if (j==(len(self.lov)-1)):
            v1_id=self.lov[j]
            v2_id=self.lov[0]
        else:
            v1_id=self.lov[j]
            v2_id=self.lov[j+1]        
        # pbc correct
        v1=t.torus_corrected_position(t.vertices[v1_id].r, vref) #v1
        v2=t.torus_corrected_position(t.vertices[v2_id].r, vref) #v2
        #return
        return [v1_id,v2_id,v1,v2]
    
    # This function returns area of one triangular element
    def dA(self, j, vref, t):
        # vertices
        [v1_id,v2_id,v1,v2]=self.triangle_vertices_pbc(j,vref,t)
        # calculate area
        ax=v2[1]*v1[2]-v1[1]*v2[2]+v1[1]*vref[2]-vref[1]*v1[2]-v2[1]*vref[2]+vref[1]*v2[2]
        ay=v1[0]*v2[2]-v2[0]*v1[2]-v1[0]*vref[2]+vref[0]*v1[2]+v2[0]*vref[2]-vref[0]*v2[2]
        az=v2[0]*v1[1]-v1[0]*v2[1]+v1[0]*vref[1]-vref[0]*v1[1]-v2[0]*vref[1]+vref[0]*v2[1]
        # return
        return 0.5*math.sqrt(ax*ax + ay*ay + az*az)
    
    # This function returns polygon area
    def polygon_area(self,t):
        self.reset_lov(t)
        sum=0.0
        vref=[self.center[0],self.center[1],self.center[2]]
        for j in range(len(self.lov)):
            sum+=self.dA(j,vref,t)
        return sum
    
    # This function calculates force contribution of one triangular element due to surface tension
    def force_surfaceTensionP(self, j, vref, t):
        # vertices
        [v1_id,v2_id,v1,v2]=self.triangle_vertices_pbc(j,vref,t)
        # force calculation
        ax=v2[1]*v1[2]-v1[1]*v2[2]+v1[1]*vref[2]-vref[1]*v1[2]-v2[1]*vref[2]+vref[1]*v2[2]
        ay=v1[0]*v2[2]-v2[0]*v1[2]-v1[0]*vref[2]+vref[0]*v1[2]+v2[0]*vref[2]-vref[0]*v2[2]
        az=v2[0]*v1[1]-v1[0]*v2[1]+v1[0]*vref[1]-vref[0]*v1[1]-v2[0]*vref[1]+vref[0]*v2[1]
        fac_area=0.5*math.sqrt(ax*ax + ay*ay + az*az)
        # gradient
        deriv = 0.25*self.g0/(fac_area+1.0e-8)
        ms=1/(1.*len(self.lov))
        a=[0,0,0]
        a[0] = az*(v2[1] - vref[1]) - ay*(v2[2] - vref[2])
        a[1] = ax*(v2[2] - vref[2]) - az*(v2[0] - vref[0])
        a[2] = ay*(v2[0] - vref[0]) - ax*(v2[1] - vref[1])
        b=[0,0,0]
        b[0] = az*(v1[1] - vref[1]) - ay*(v1[2] - vref[2])
        b[1] = ax*(v1[2] - vref[2]) - az*(v1[0] - vref[0])
        b[2] = ay*(v1[0] - vref[0]) - ax*(v1[1] - vref[1])
        c=[0,0,0]
        for k in range(3):
            c[k]=deriv*ms*(-a[k]+b[k])
        # vertices update
        #v1
        for k in range(3):
            t.vertices[v1_id].F[k]+=deriv*a[k]
        #v2    
        for k in range(3):
            t.vertices[v2_id].F[k]+=-deriv*b[k]
        #p_vertices
        for vr in self.lov:
            for k in range(3):
                t.vertices[vr].F[k]+=c[k]
    
    # This function calculates forces due to surface tension of one polygon
    def force_surfaceTension(self, t):
        self.reset_lov(t)
        vref=[self.center[0],self.center[1],self.center[2]]
        for j in range(len(self.lov)):
            self.force_surfaceTensionP(j,vref,t)
        return self.g0*self.area
    
    # This function returns volume element corresponding to one triangular element
    def dV(self, j, vref, t, chk):
        # vertices
        [v1_id,v2_id,v1,v2]=self.triangle_vertices_pbc(j,vref,t)
        if chk==1:
            vert3=self.polygon_center(t)
        else:
            vert3=[self.center[0],self.center[1],self.center[2]]
        v3=t.torus_corrected_position(vert3, vref)
        # return
        return (-v1[2]*v2[1]*v3[0] + v1[1]*v2[2]*v3[0] + v1[2]*v2[0]*v3[1] - v1[0]*v2[2]*v3[1] - v1[1]*v2[0]*v3[2] + v1[0]*v2[1]*v3[2])/6.
    
    # This function calculates forces contribution of one triangular element due to cell volume preservation
    def force_volumeSpringP(self,j,vref,ori,cell_id,cell_vol,t):
        # vertices
        [v1_id,v2_id,v1,v2]=self.triangle_vertices_pbc(j,vref,t)     
        cent=[self.center[0],self.center[1],self.center[2]]
        v3=t.torus_corrected_position(cent, vref) #v3
        # gradient
        deriv = -2*ori*t.kV*(cell_vol-t.cells[cell_id].V0)
        ms=1/(6.*len(self.lov))
        #v1
        t.vertices[v1_id].F[0]+=deriv*(-v2[2]*v3[1] + v2[1]*v3[2])/6.
        t.vertices[v1_id].F[1]+=deriv*( v2[2]*v3[0] - v2[0]*v3[2])/6.
        t.vertices[v1_id].F[2]+=deriv*(-v2[1]*v3[0] + v2[0]*v3[1])/6.
        #v2
        t.vertices[v2_id].F[0]+=deriv*( v1[2]*v3[1] - v1[1]*v3[2])/6.
        t.vertices[v2_id].F[1]+=deriv*(-v1[2]*v3[0] + v1[0]*v3[2])/6.
        t.vertices[v2_id].F[2]+=deriv*( v1[1]*v3[0] - v1[0]*v3[1])/6.
        #p_vertices
        a=[0,0,0]
        a[0] = deriv * ms * (v1[1]*v2[2] - v2[1]*v1[2])
        a[1] = deriv * ms * (v2[0]*v1[2] - v1[0]*v2[2])
        a[2] = deriv * ms * (v1[0]*v2[1] - v2[0]*v1[1])
        for vr in self.lov:
            for k in range(3):
                t.vertices[vr].F[k]+=a[k]
    
    # This function corrects vertex positions upon a TE transition
    def correct_vertex_positions_TE(self, t, v1ID, v2ID, c4ID):
        # v1 & v2
        v1_id = v1ID
        v2_id = v2ID
        # Mid point of the target triangle is assigned to the v1
        t.vertices[v1_id].r=[self.center[0],self.center[1],self.center[2]]
        # Central point of a neighbor cell is assigned as v2 vertex (later modified)
        c4 = c4ID
        t.vertices[v2_id].r=[t.cells[c4].center[0],t.cells[c4].center[1],t.cells[c4].center[2]]
        # a vector from v1 to the center of the neighbour cell (v2 now)
        v2corr=t.torus_corrected_position(t.vertices[v2_id].r, t.vertices[v1_id].r)
        dr=[0,0,0]
        for j in range(3):
            dr[j]=v2corr[j]-t.vertices[v1_id].r[j]
        len_temp=math.sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])
        # v1temp
        v1temp=[0,0,0]
        for j in range(3):
            v1temp[j]=t.vertices[v1_id].r[j]
        # New coordinates of v2
        for j in range(3):
            t.vertices[v1_id].r[j]=v1temp[j]-t.new_len*(dr[j]/len_temp)
            t.vertices[v2_id].r[j]=v1temp[j]+t.new_len*(dr[j]/len_temp)
        # taking care of periodic boundary condition
        v1=t.torus_on_position(t.vertices[v1_id].r)
        v2=t.torus_on_position(t.vertices[v2_id].r)
        for j in range(3):
            t.vertices[v1_id].r[j]=v1[j]
            t.vertices[v2_id].r[j]=v2[j]
    
    # This function performs TE transition
    def TE_transition(self, G, t):
        # pattern matching
        v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5 = T1.patternMatching_TE(G, self.id)
        impossible_TE = Q.impossible_TE(G, p1, p2, p3)
        if impossible_TE[0] == 0 and impossible_TE[1] == 0 and impossible_TE[2] == 0:
            # graph transformations
            T1.graphTransformation_TE(G, v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5)
            # correct new vertices position
            self.correct_vertex_positions_TE(t, v1[1], v2[1], c4[1])
            # updating tissue
            t.update_TE(G, v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5)
            # remove nonexistent entities
            t.remove_nonexistent(G)
            # counter
            t.TE_count+=1
        else:
            print("impossible TE")


class cell:
    def __init__(self, id:int, lop:int, loo:int, t):
        # essential variables
        self.id = id
        self.lop = lop #list of polygons
        self.loo = loo #list of polygon orientations
        self.exist = 1
        self.growing = 0
        self.V0 = 1 # preferred cell volume
        self.growth_rate = 1
        self.boundary_cell = 0
        # reset cell properties
        self.reset_lov(t)
        self.reset_prop()
        # add it to list of cells
        t.cells.append(self) 
    
    # This function resets cell properties
    def reset_prop(self):
        self.center = [0, 0, 0]
        self.volume = 0
    
    # This function resets list of vertices
    def reset_lov(self,t):
        self.lov = []
        if len(self.lop) > 0:
            cell_vertices_all = []
            for i in range(len(self.lop)):
                t.polygons[self.lop[i]].reset_lov(t)
                for j in range(len(t.polygons[self.lop[i]].lov)):
                    cell_vertices_all.append(t.polygons[self.lop[i]].lov[j])
            # keep only unique vertices
            self.lov = list(set(cell_vertices_all))
    
    # This function calculates cell volume
    def cell_volume(self,t,chk):
        volsum = 0.0
        if chk == 1:
            vref = self.cell_center(t)
        else:
            vref = [self.center[0],self.center[1],self.center[2]]
        for i in range(len(self.lop)):
            t.polygons[self.lop[i]].reset_lov(t)
            for j in range(len(t.polygons[self.lop[i]].lov)):
                dv = t.polygons[self.lop[i]].dV(j,vref,t,chk)
                volsum += self.loo[i]*dv
        return volsum
    
    # This function calculates forces due to volume preservation
    def force_volumeSpring(self,t):
        vref = [self.center[0],self.center[1],self.center[2]]
        for i in range(len(self.lop)):
            t.polygons[self.lop[i]].reset_lov(t)
            for j in range(len(t.polygons[self.lop[i]].lov)):
                t.polygons[self.lop[i]].force_volumeSpringP(j,vref,self.loo[i],self.id,self.volume,t)
        return t.kV*(self.volume-self.V0)*(self.volume-self.V0)
    
    # This function calculates cell center
    def cell_center(self, t):
        self.reset_lov(t)
        vref = [t.vertices[self.lov[0]].r[0], t.vertices[self.lov[0]].r[1], t.vertices[self.lov[0]].r[2]]
        # calculate cell center
        cent = [vref[0], vref[1], vref[2]]
        for i in range(1, len(self.lov)):
            v = t.torus_corrected_position(t.vertices[self.lov[i]].r, vref)
            for j in range(3):
                cent[j] += v[j]
        return t.torus_on_position([cent[0]/(len(self.lov)*1.0), cent[1]/(len(self.lov)*1.0), cent[2]/(len(self.lov)*1.0)])
    
    # This function calculates the cell center manually from the graph
    def cell_center_manual(self, G, t):
        ps_of_c = list(G.predecessors(("Cell", self.id)))
        vs_of_c = []
        for p_of_c in ps_of_c:
            es_of_p = list(G.predecessors(p_of_c))
            for e_of_p in es_of_p:
                vs_of_c += list(G.predecessors(e_of_p))
        vs = list(set(vs_of_c))
        vref=[t.vertices[vs[0][1]].r[0], t.vertices[vs[0][1]].r[1], t.vertices[vs[0][1]].r[2]]
        cent=[vref[0], vref[1], vref[2]]
        for i in range(1, len(vs)):
            v = t.torus_corrected_position(t.vertices[vs[i][1]].r, vref)
            for j in range(3):
                cent[j] += v[j]
        return t.torus_on_position([cent[0]/(len(vs)*1.0),cent[1]/(len(vs)*1.0),cent[2]/(len(vs)*1.0)])
    
    # Correct vertex positions for the cell division
    def correct_vertex_positions_D(self, G, t, vs_new, cut_es):
        for i in range(len(cut_es)):
            preds_e1 = list(G.predecessors(cut_es[i]))
            for pred_e1 in preds_e1:
                if pred_e1[0] == "Div":
                    nodes_div = list(G.successors(pred_e1))
                    for node_div in nodes_div:
                        name = list(G.get_edge_data(pred_e1, node_div))[0]
                        if name == "e1" and cut_es[i] == node_div:
                                div = pred_e1
                                break
            nodes_div = list(G.successors(div))
            for node_div in nodes_div:
                name = list(G.get_edge_data(div, node_div))[0]
                if name == "v3":
                    v3 = node_div
                    break
            v3_coord = vs_new[i]
            t.vertices[v3[1]].r = v3_coord
    
    # Calculate the division geometry
    def division_geometry(self, G, t):
        ps_of_div_c = list(G.predecessors(("Cell", self.id)))
        es_of_div_c = []
        for p_of_div_c in ps_of_div_c:
            es_of_div_c += list(G.predecessors(p_of_div_c))
        es_div_cell = list(set(es_of_div_c))
        # center of mass for the dividing cell
        com = self.cell_center_manual(G, t)
        
        tc = t.tumor_center()
        u = [com[i] - tc[i] for i in range(3)]
        
        bad_cut = True
        while bad_cut:
            # define the division plane
            w0 = random.uniform(-1, 1)
            w1 = random.uniform(-1, 1)
            w2 = random.uniform(-1, 1)
            w = [w0, w1, w2]
            normal_div_plane = [u[1]*w[2] - u[2]*w[1], u[2]*w[0] - u[0]*w[2], u[0]*w[1] - u[1]*w[0]]
            
            vs_new = [] # new vertices 
            cut_es = [] # edges cut by division plane 
            cut_ps = [] # polygons cut by division plane 
            ps_cut_e = [] # dividing cell's polygons of cut edges 
            e_pps = []
            # loop over all edges of the dividing cell
            for e in es_div_cell:  
                vv_id = list(G.predecessors(e))
                v1_id = vv_id[0][1]
                v2_id = vv_id[1][1]
                v1 = t.torus_corrected_position([t.vertices[v1_id].r[0], t.vertices[v1_id].r[1], t.vertices[v1_id].r[2]], com)
                v2 = t.torus_corrected_position([t.vertices[v2_id].r[0], t.vertices[v2_id].r[1], t.vertices[v2_id].r[2]], com)
                v1_com = [v1[0]-com[0], v1[1]-com[1], v1[2]-com[2]]
                v2_com = [v2[0]-com[0], v2[1]-com[1], v2[2]-com[2]]
                n_dot_v1 = normal_div_plane[0]*v1_com[0] + normal_div_plane[1]*v1_com[1] + normal_div_plane[2]*v1_com[2] # n*(v1-com)
                n_dot_v2 = normal_div_plane[0]*v2_com[0] + normal_div_plane[1]*v2_com[1] + normal_div_plane[2]*v2_com[2] # n*(v2-com)
                if n_dot_v1*n_dot_v2 < 0: # one vertex above and one vertex below the division plane -> cut edge
                    cut_e = e
                    cut_es.append(cut_e)
                    v_up = v1
                    v_lo = v2
                    com_v_up = [com[0]-v_up[0], com[1]-v_up[1], com[2]-v_up[2]]
                    v_lo_v_up = [v_lo[0]-v_up[0], v_lo[1]-v_up[1], v_lo[2]-v_up[2]]
                    numerator = normal_div_plane[0]*com_v_up[0] + normal_div_plane[1]*com_v_up[1] + normal_div_plane[2]*com_v_up[2]
                    denominator = normal_div_plane[0]*v_lo_v_up[0] + normal_div_plane[1]*v_lo_v_up[1] + normal_div_plane[2]*v_lo_v_up[2]
                    t_param = numerator/denominator # definition of a line in 3D
                    # new vertex coordinates considering periodic boundaries
                    v_new = t.torus_on_position([v_up[0] + t_param*(v_lo[0] - v_up[0]), v_up[1] + t_param*(v_lo[1] - v_up[1]), v_up[2] + t_param*(v_lo[2] - v_up[2])])
                    vs_new.append(v_new)
                    pp_cut_e = Q.hos_of_lo_hho(G, cut_e, ("Cell", self.id)) # dividing cell's polygons of the cut edge
                    ps_cut_e.append(pp_cut_e[0])
                    ps_cut_e.append(pp_cut_e[1])
                    e_pps.append([cut_e, [pp_cut_e[0], pp_cut_e[1]]])
            # unique cut polygons
            cut_ps = list(set(ps_cut_e))
            
            # check if there are any problems with the cut
            e_this = e_pps[0][0]
            p_next = e_pps[0][1][1]
            es_done = [e_this]
            while p_next != e_pps[0][1][0]:
                for e_pp in e_pps:
                    if e_pp[0] not in es_done:
                        for p12 in e_pp[1]:
                            if p12 == p_next:
                                e_this = e_pp[0]
                                es_done.append(e_this)
                                if p12 == e_pp[1][0]:
                                    p_next = e_pp[1][1]
                                    break
                                elif p12 == e_pp[1][1]:
                                    p_next = e_pp[1][0]
                                    break
            bad_cut = False
            if len(es_done) != len(cut_es) or len(cut_es) != len(cut_ps):
                bad_cut = True
            if bad_cut == False:
                for cut_e1 in cut_es:
                    if bad_cut == True:
                        break
                    for cut_e2 in cut_es:
                        if cut_e1 != cut_e2:
                            ps1 = list(G.successors(cut_e1))
                            ps2 = list(G.successors(cut_e2))
                            ps_common = list(set(ps1) & set(ps2))
                            if ps_common:
                                if ps_common[0] not in cut_ps:
                                    bad_cut = True
                                    break
        # return the needed information about the division geometry
        return vs_new, cut_es, cut_ps
    
    # This function performs cell division
    def Division(self, G, t, DB):
        # division geometry
        vs_new, cut_es, cut_ps = self.division_geometry(G, t)
        num_cut_ps = len(cut_ps)
        # add new nodes if needed
        DB.add_nodes_D(G, t, num_cut_ps)
        # pattern matching & graph transformations
        success = D.PM_GT_D(G, ("Cell", self.id), cut_ps, num_cut_ps, cut_es)
        if success:
            # find all lower polygons and the number of polygons that are not 'p4' polygons
            lower_ps, num_extra_ps = D.lower_polygons(G, ("Cell", self.id), num_cut_ps)
            # add extra nodes if needed
            DB.add_extra_nodes_D(G, t, num_extra_ps)
            # graph transformations for lower polygons
            D.lower_ps_transformations(G, ("Cell", self.id), lower_ps, num_cut_ps)
            # correct vertex positions
            self.correct_vertex_positions_D(G, t, vs_new, cut_es)
            # updating tissue
            t.update_D(G, num_cut_ps, num_extra_ps)
            # counter
            t.D_count+=1
            # release the helping nodes
            for div_id in range(num_cut_ps + num_extra_ps):
                nodes_div = list(G.successors(("Div", div_id)))
                for node_div in nodes_div:
                    G.remove_edge(("Div", div_id), node_div)
                G.add_edge("Free", ("Div", div_id), "CONTAINS") # free all Div nodes
        return success


class tissue:
    def __init__(self, vt3d:str):
        # essential variables
        self.vt3d = vt3d
        self.Lxy = []
        self.vertices = []
        self.edges = []
        self.polygons = []
        self.cells = []
        # number of removed entities
        self.num_removed_v = 0
        self.num_removed_e = 0
        self.num_removed_p = 0
        self.num_removed_c = 0
        # model parameters
        self.h = 0.001 # default time step
        self.Nsim = 0 # simulation number
        self.Tmax = 100 # total simulation time
        self.lengthTH = 0.01 # length threshold
        self.clockTH = 0.02 # clock threshold
        self.new_len = 0.001 # length of new edge
        self.dg0 = 1 # baseline line tension
        self.sig = 0.25 # magnitude of fluctuations
        self.kM = 1 # myosin turnover rate
        self.kV = 100 # volume compressibility modulus
        self.outFreq = 1.0 # output frequency
        self.division_rate = 0.01 # division rate 'k_div'
        self.g0_layer = 1 # tension between live and necrotic layers
        self.live_layers = 1 # number of live layers
        self.Nmax = 2000 # max number of cells
        self.g0_boundary = 1.0 # boundary tension
        self.topological_transition_happened = True
        # properties
        self.reset_prop()
        # initialize
        self.setup_from_vt3d()
    
    # This function resets tissue properties
    def reset_prop(self):
        self.Time = 0
        self.timeCount = 0
        self.outFolder = ""
        self.fileCount = 0
        self.wT = 0
        self.Nv = 0
        self.Ne = 0
        self.Np = 0
        self.Nc = 0
        self.Nd = 0
        self.VT_count = 0
        self.ET_count = 0
        self.TE_count = 0
        self.D_count = 0
        self.nET = 0
        self.nTE = 0
        self.nD = 0
    
    # This function sets up tissue from the input .vt3d file
    def setup_from_vt3d(self):
        self.Lxy=[0,0,0]
        with open(self.vt3d, "r") as f:
            # nr vertices, edges, polygons
            (self.Nv, self.Ne, self.Np, self.Nc) = [int(x) for x in f.readline().split()]
            f.readline()
            # box size
            (self.Lxy[0],self.Lxy[1],self.Lxy[2])=[float(x) for x in f.readline().split()]
            f.readline()
            # vertices
            for count in range(self.Nv):
                (xx, yy, zz)= [float(x) for x in f.readline().split()]
                vertex(count,[xx,yy,zz],self)
            f.readline()
            # edges
            for count in range(self.Ne):
                (vv1, vv2)= [x for x in f.readline().split()]
                edge(count, [int(vv1)-1, int(vv2)-1], self)
            f.readline()
            # polygons
            for count in range(self.Np):
                b = f.readline().split()
                edges=[]
                orientations=[]
                for i in range(len(b)):
                    if int(b[i])!=0 and i>0:
                        if int(b[i])>0:
                            edges.append(int(b[i])-1)
                            orientations.append(1)
                        else:
                            edges.append(abs(int(b[i])+1))
                            orientations.append(-1)
                polygon(count, edges, orientations, self)
            f.readline()
            # cells
            for count in range(self.Nc):            
                b = f.readline().split()
                polygons=[]
                orientations=[]
                for i in range(1,len(b)):
                    if int(b[i])!=0:
                        if int(b[i])>0:
                            polygons.append(int(b[i])-1)
                            orientations.append(1) # Ordinary: +1, -1
                        else:
                            polygons.append(abs(int(b[i])+1))
                            orientations.append(-1) # Ordinary: -1, +1
                cell(count, polygons, orientations, self)
        f.close()
    
    # This function generates a script to setup database
    def generate_DB_setup_script(self, G):
        G.add_node("T1")
        G.add_node("Free")
        n_count=0
        for c in self.cells:
            if (c.exist==1):
                G.add_node(("Cell", n_count))
                n_count +=1
        n_count=0
        for p in self.polygons:
            if (p.exist==1):
                G.add_node(("Polygon", n_count))
                n_count +=1
        n_count=0
        for e in self.edges:
            if (e.exist==1):
                G.add_node(("Edge", n_count))
                n_count +=1
        n_count=0
        for v in self.vertices:
            if (v.exist==1):
                G.add_node(("Vertex", n_count))
                n_count +=1        
        # creating relationships
        # for Vertex -> Edge
        for count in range(len(self.edges)):
            if(len(self.edges[count].e)>0):
                G.add_edge(("Vertex", self.edges[count].e[0]), ("Edge", count), ("IS_PART_OF", 1))
                G.add_edge(("Vertex", self.edges[count].e[1]), ("Edge", count), ("IS_PART_OF", -1))
        # for Edge -> Polygon
        for count in range(len(self.polygons)):
            edges=self.polygons[count].loe
            orientations=self.polygons[count].loo
            if(len(edges)>0):
                for i in range(len(edges)):
                    edge=edges[i]
                    orientation=orientations[i]
                    if orientation>0:
                        G.add_edge(("Edge", abs(edge)), ("Polygon", count), ("IS_PART_OF", 1))
                    else:
                        G.add_edge(("Edge", abs(edge)), ("Polygon", count), ("IS_PART_OF", -1))
        # for Polygon -> Cell
        for count in range(len(self.cells)):
            polygons=self.cells[count].lop
            orientations=self.cells[count].loo
            if(len(polygons)>0):
                for i in range(len(polygons)):
                    polygon=polygons[i]
                    orientation=orientations[i]
                    if orientation>0:
                        G.add_edge(("Polygon", abs(polygon)), ("Cell", count), ("IS_PART_OF", 1))
                    else:
                        G.add_edge(("Polygon", abs(polygon)), ("Cell", count), ("IS_PART_OF", -1))
        # for Cell -> Tissue
        for count in range(len(self.cells)):
            if(count<self.Nc):
                G.add_edge(("Cell", count), ("Tissue", 0), "IS_PART_OF")
    
    # This function corrects position due to periodic boundary conditions
    def torus_on_position(self,v):
        dxdydz=[0,0,0]
        for i in range(3):
            if v[i]<0:
                dxdydz[i]=self.Lxy[i]
            elif v[i]>self.Lxy[i]:
                dxdydz[i]=-self.Lxy[i]
        return [v[0]+dxdydz[0],v[1]+dxdydz[1],v[2]+dxdydz[2]]
    
    # This function calculates displacement of one vertex vs another due to periodic boundary conditions
    def torus_corrected_position(self, v, vref):
        dxdydz=[0,0,0]
        for i in range(3):
            if abs(v[i]-vref[i])>0.5*self.Lxy[i]:
                if v[i]<vref[i]:
                    dxdydz[i]=self.Lxy[i]
                elif v[i]>vref[i]:
                    dxdydz[i]=-self.Lxy[i]
        return [v[0]+dxdydz[0],v[1]+dxdydz[1],v[2]+dxdydz[2]]
    
    # This function calculates forces on vertices
    def forces(self):
        # line energy
        wActive=0
        for e in self.edges:
            if e.exist==1 and e.lineTension!=0:
                wActive+=e.force_lineTension(self)
        # surface energy
        wA=0
        for p in self.polygons:
            if p.exist==1 and p.g0!=0:
                wA+=p.force_surfaceTension(self)
        # volume conservation
        wV=0
        for c in self.cells:
            if c.exist==1 and self.kV!=0:
                wV+=c.force_volumeSpring(self)
        # total energy
        self.wT=wA+wV
        return [wActive,wA,wV]
    
    # This function updates objects' calculable properties
    def update_props(self, when):
        # edge lengths
        for e in self.edges:
            if e.exist==1:
                e.center=e.edge_center(self)
                elnew=e.edge_length(self)
                if when=='after':
                    e.e_dl=elnew-e.length
                    e.clock+=self.h
                e.length=elnew
        # polygon areas and centers
        for p in self.polygons:
            if p.exist==1:
                p.center=p.polygon_center(self)
                pAnew=p.polygon_area(self)
                if when=='after':
                    p.p_dA=pAnew-p.area
                p.area=pAnew
        # cell volumes
        vMIN=100
        for c in self.cells:
            if c.exist==1:
                c.center=c.cell_center(self)
                cVnew=c.cell_volume(self, 0)
                if cVnew<vMIN:
                    vMIN=cVnew
                c.volume=cVnew
        # return
        return vMIN
    
    # This function propagates tensions
    def propagate_tensions(self):
        for e in self.edges:
            if e.exist==1:
                e.lineTension += -self.h*self.kM*(e.lineTension-self.dg0) + math.sqrt(2*self.h*e.sigma*e.sigma*self.kM)*random.normalvariate(0,1)
                if e.lineTension<0:
                    e.lineTension=0
    
    # This function propagates the tissue from t to t+dt
    def solver(self):
        # update props
        vv=self.update_props('before')     
        # calculate forces
        for v in self.vertices:
            if v.exist==1:
                v.F=[0,0,0]
        [wActive,wA,wV]=self.forces()
        # propagate vertices
        drMAX=0
        for v in self.vertices:
            if v.exist==1:
                # displace vertices
                for i in range(3):
                    v.r[i]+=self.h*v.F[i]            
                # check max displacement
                dr=self.h*math.sqrt(v.F[0]*v.F[0]+v.F[1]*v.F[1]+v.F[2]*v.F[2])
                if dr>drMAX:
                    drMAX=dr
        # proapgate tensions
        self.propagate_tensions()
        # update properties
        vMIN=self.update_props('after')  
        # time
        self.Time += self.h
        # output
        if self.Time >= self.timeCount:
            with open(self.outFolder+"/out_{nsim}_{nmax}_{sigma}_{drate}_{tens}_{lays}_{g0b}.txt".format(nsim=self.Nsim, nmax=self.Nmax, sigma=self.sig, drate=self.division_rate, tens=self.g0_layer, lays=self.live_layers, g0b=self.g0_boundary), 'a') as filebv:
                filebv.write(str(round(self.Time,6))+"\t"+str(self.Nc - self.num_removed_c)+"\t"+str(self.wT)+"\t"+str(wA)+"\t"+str(wV)+"\t"+str(self.ET_count)+"\t"+str(self.TE_count)+"\t"+str(self.D_count)+"\t"+str(drMAX)+"\t"+str(vv)+"\t"+str(vMIN)+'\n')
            filebv.close()
            self.fileCount += 1
            self.out_vt3d(self.outFolder,self.fileCount)
            self.timeCount += self.outFreq
            self.nET = self.ET_count
            self.nTE = self.TE_count
            self.nD = self.D_count
            print("[t="+str(self.Time)+"] Intermediate result nr. "+str(self.fileCount)+" stored in "+self.outFolder)
    
    # This function outputs a .vt3d file
    def out_vt3d(self, path:str, filenum:int):
        with open(path+'/out_{nsim}_{nmax}_{sigma}_{drate}_{tens}_{lays}_{g0b}_{fnum}.vt3d'.format(nsim=self.Nsim, nmax=self.Nmax, sigma=self.sig, drate=self.division_rate, tens=self.g0_layer, lays=self.live_layers, g0b=self.g0_boundary, fnum=filenum), 'w') as f3:
            # nr of elements
            list=[len(self.vertices),len(self.edges),len(self.polygons),len(self.cells)]             
            for i in range(len(list)):
                f3.write(str(list[i])+"\t")              
            f3.write("\n\n")
            # simulation box
            for i in range(len(self.Lxy)):
                f3.write(str(self.Lxy[i])+"\t")
            f3.write("\n\n")
            # vertices
            for v in self.vertices:
                if v.exist==1:
                    f3.write(str(v.exist)+"\t")
                    list=v.r
                    for i in range(len(list)):
                        f3.write(str(list[i])+"\t")
                    f3.write("\n")
                else:
                    f3.write(str(v.exist)+"\t"+str(0)+"\t"+str(0)+"\t"+str(0)+"\t")
                    f3.write("\n")
                    
            f3.write("\n")
            # edges
            for ed in self.edges:
                if ed.exist==1:
                    f3.write(str(ed.exist)+"\t")
                    list=ed.e
                    for i in range(len(list)):
                        f3.write(str(list[i]+1)+"\t")
                    f3.write("\n")
                    
                else:
                    f3.write(str(ed.exist)+"\t"+str(0)+"\t"+str(0)+"\t")
                    f3.write("\n")
            f3.write("\n")
            # polygons
            for p in self.polygons:
                if p.exist==1:
                    f3.write(str(p.exist)+"\t")
                    listP=p.loe
                    listO=p.loo
                    for i in range(len(listP)):
                        f3.write(str(listO[i]*(listP[i]+1))+"\t")
                    f3.write("\n")
                else:
                    f3.write(str(p.exist)+"\t"+str(0)+"\t"+str(0)+"\t"+str(0)+"\t")
                    f3.write("\n")
            f3.write("\n")
            # cells
            for c in self.cells:
                if c.exist==1:
                    listP=c.lop
                    listO=c.loo
                    for i in range(len(listP)):
                        f3.write(str(listO[i]*(listP[i]+1))+"\t")
                    f3.write("\n")
        f3.close()

    # This function performs topological transformations
    def topological_transitions(self, G, DB):
        self.topological_transition_happened = False
        for e in self.edges:
            if e.exist==1 and e.clock>self.clockTH:
                if e.length<self.lengthTH and e.e_dl<0:
                    # ET transition
                    triangle = False
                    ps_of_e = list(G.successors(("Edge", e.id)))
                    for p_of_e in ps_of_e:
                        es_of_p = list(G.predecessors(p_of_e))
                        if len(es_of_p) == 3:
                            triangle = True
                            p_triangle = p_of_e
                            break
                    if triangle == False:
                        if Q.count_shared_polygons_c4c5(e.id, G) == 0:
                            e.ET_transition(G, self, DB)
                            self.topological_transition_happened = True
                    # TE transition
                    else:
                        p = self.polygons[p_triangle[1]]
                        elist = list(G.predecessors(p_triangle))
                        e1 = self.edges[elist[0][1]]
                        e2 = self.edges[elist[1][1]]
                        e3 = self.edges[elist[2][1]]
                        if e1.clock>self.clockTH and e2.clock>self.clockTH and e3.clock>self.clockTH:
                            if e1.length<self.lengthTH and e2.length<self.lengthTH and e3.length<self.lengthTH:
                                neighbor_triangle = False
                                for e0 in elist:
                                    p0s = list(G.successors(e0))
                                    for p0 in p0s:
                                        if p0 != p_triangle:
                                            if len(list(G.predecessors(p0))) == 3:
                                                neighbor_triangle = True
                                                break
                                    if neighbor_triangle == True:
                                        break
                                if neighbor_triangle == False:
                                    if p.p_dA<0:
                                        p.TE_transition(G, self)
                                        self.topological_transition_happened = True
        self.update_props('before')
        for c in self.cells:
            if c.exist==1:
                if c.cell_volume(self, 1) >= 2:
                    success = c.Division(G, self, DB)
                    self.topological_transition_happened = True
                    if not success:
                        return False
                    c.growing = 0
                    c.V0 = 1
        return True
    
    # This function updates edges, polygons, and cells
    def update_edges_polygons_cells(self, G, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5):
        # Fetch data
        resE = []
        for e in [e1, e2, e3, e4, e5, e6, e7, e8, e9]:
            try:
                vs_of_e = list(G.predecessors(e))
                vs_info = []
                for v in vs_of_e:
                    r_ve = list(G.get_edge_data(v, e))[0][1]
                    vs_info.append([v[1], r_ve])
                resE.append([e[1], vs_info])
            except nx.exception.NetworkXError:
                continue
        resP = []
        for p in [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]:
            try:
                es_of_p = list(G.predecessors(p))
                es_info = []
                for e in es_of_p:
                    r_ep = list(G.get_edge_data(e, p))[0][1]
                    es_info.append([e[1], r_ep])
                resP.append([p[1], es_info])
            except nx.exception.NetworkXError:
                continue
        resC = []
        for c in [c1, c2, c3, c4, c5]:
            try:
                ps_of_c = list(G.predecessors(c))
                ps_info = []
                for p in ps_of_c:
                    r_pc = list(G.get_edge_data(p, c))[0][1]
                    ps_info.append([p[1], r_pc])
                resC.append([c[1], ps_info])
            except nx.exception.NetworkXError:
                continue
        # edges
        for edg in resE:
            for i in range(len(edg[1])):
                if edg[1][i][1]==-1:
                    tail=edg[1][i][0]
                else:
                    head=edg[1][i][0]
            self.edges[edg[0]].e=[head,tail]
        # polygons
        for poly in resP:
            loe=[]
            loo=[]
            for i in range(len(poly[1])):
                loe.append(poly[1][i][0])
                loo.append(poly[1][i][1])
            self.polygons[poly[0]].loe=loe
            self.polygons[poly[0]].loo=loo
        # cells
        for cqs in resC:
            lop=[]
            loo=[]
            for i in range(len(cqs[1])):
                lop.append(cqs[1][i][0])
                loo.append(cqs[1][i][1])
            self.cells[cqs[0]].lop=lop
            self.cells[cqs[0]].loo=loo
    
    # This function updates edges, polygons, and cells after division
    def update_edges_polygons_cells_D(self, G, num_cut_ps, num_extra_ps):
        for div_id in range(num_cut_ps + num_extra_ps):
            es = []
            ps = []
            cs = []
            nodes_div = list(G.successors(("Div", div_id)))
            for node_div in nodes_div:
                if node_div[0] == "Edge":
                    es.append(node_div)
                if node_div[0] == "Polygon":
                    ps.append(node_div)
                if node_div[0] == "Cell":
                    cs.append(node_div)
            # Fetch data
            resE = []
            for e in es:
                vs_of_e = list(G.predecessors(e))
                vs_info = []
                for v in vs_of_e:
                    if v[0] != "Div":
                        r_ve = list(G.get_edge_data(v, e))[0][1]
                        vs_info.append([v[1], r_ve])
                resE.append([e[1], vs_info])
            resP = []
            for p in ps:
                es_of_p = list(G.predecessors(p))
                es_info = []
                for e in es_of_p:
                    if e[0] != "Div":
                        r_ep = list(G.get_edge_data(e, p))[0][1]
                        es_info.append([e[1], r_ep])
                resP.append([p[1], es_info])
            resC = []
            for c in cs:
                ps_of_c = list(G.predecessors(c))
                ps_info = []
                for p in ps_of_c:
                    if p[0] != "Div":
                        r_pc = list(G.get_edge_data(p, c))[0][1]
                        ps_info.append([p[1], r_pc])
                resC.append([c[1], ps_info])
            # edges
            for edg in resE:
                for i in range(len(edg[1])):
                    if edg[1][i][1]==-1:
                        tail=edg[1][i][0]
                    else:
                        head=edg[1][i][0]
                self.edges[edg[0]].e=[head,tail]
            # polygons
            for poly in resP:
                loe=[]
                loo=[]
                for i in range(len(poly[1])):
                    loe.append(poly[1][i][0])
                    loo.append(poly[1][i][1])
                self.polygons[poly[0]].loe=loe
                self.polygons[poly[0]].loo=loo
            # cells
            for cqs in resC:
                lop=[]
                loo=[]
                for i in range(len(cqs[1])):
                    lop.append(cqs[1][i][0])
                    loo.append(cqs[1][i][1])
                self.cells[cqs[0]].lop=lop
                self.cells[cqs[0]].loo=loo
    
    # This function updates the tissue upon a ET transition
    def update_ET(self, G, v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5):
        self.update_edges_polygons_cells(G, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5)
        # find nodes
        v3=self.vertices[v3[1]]
        e7=self.edges[e7[1]]
        e8=self.edges[e8[1]]
        e9=self.edges[e9[1]]
        p10=self.polygons[p10[1]]
        # .exist=1
        v3.exist=1
        e7.exist=1
        e8.exist=1
        e9.exist=1
        p10.exist=1
        # reset properties
        v3.reset_prop(self)
        e7.reset_prop(self)
        e8.reset_prop(self)
        e9.reset_prop(self)
        p10.reset_prop(self)
    
    # This function updates the tissue upon a TE transition
    def update_TE(self, G, v1, v2, v3, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5):
        self.update_edges_polygons_cells(G, e1, e2, e3, e4, e5, e6, e7, e8, e9, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, c1, c2, c3, c4, c5)
        # find nodes
        v3=self.vertices[v3[1]]
        e7=self.edges[e7[1]]
        e8=self.edges[e8[1]]
        e9=self.edges[e9[1]]
        p10=self.polygons[p10[1]]
        # .exist=0
        v3.exist=0
        e8.exist=0
        e9.exist=0
        p10.exist=0
        # reset new_edge
        e7.reset_prop(self)
    
    # Update the database after cell division
    def update_D(self, G, num_cut_ps, num_extra_ps):
        # update edges, polygons, cells
        self.update_edges_polygons_cells_D(G, num_cut_ps, num_extra_ps)   
        # find nodes
        nodes_div = list(G.successors(("Div", 0)))
        for node_div in nodes_div:
            name = list(G.get_edge_data(("Div", 0), node_div))[0]
            if name == "p3":
                p3 = self.polygons[node_div[1]]
            if name == "c2":
                c2 = self.cells[node_div[1]]
        # .exist=1
        p3.exist=1
        c2.exist=1
        c2.growing=0
        c2.V0=1
        c2.growth_rate=1
        # reset properties
        p3.reset_prop(self)
        c2.reset_prop()
        for div_id in range(num_cut_ps):
            # find nodes
            nodes_div = list(G.successors(("Div", div_id)))
            for node_div in nodes_div:
                name = list(G.get_edge_data(("Div", div_id), node_div))[0]
                if name == "v3":
                    v3 = self.vertices[node_div[1]]
                if name == "e3":
                    e3 = self.edges[node_div[1]]
                if name == "e5":
                    e5 = self.edges[node_div[1]]
                if name == "p4":
                    p4 = self.polygons[node_div[1]]
            # .exist=1
            v3.exist=1
            e3.exist=1
            e5.exist=1
            p4.exist=1
            # reset properties
            v3.reset_prop(self)
            e3.reset_prop(self)
            e5.reset_prop(self)
            p4.reset_prop(self)

    # Remove nonexistent entities from the graph
    def remove_nonexistent(self, G):
        for p in self.polygons:
            if p.exist == 1:
                if len(list(G.successors(("Polygon", p.id)))) == 0:
                    G.remove_node(("Polygon", p.id))
                    p.exist = 0
                    self.num_removed_p += 1
        for e in self.edges:
            if e.exist == 1:
                if len(list(G.successors(("Edge", e.id)))) == 0:
                    G.remove_node(("Edge", e.id))
                    e.exist = 0
                    self.num_removed_e += 1
        for v in self.vertices:
            if v.exist == 1:
                if len(list(G.successors(("Vertex", v.id)))) == 0:
                    G.remove_node(("Vertex", v.id))
                    v.exist = 0
                    self.num_removed_v += 1

    # Create a tissue ellipsoid with the semi-axes 'ax', 'by' and 'cz'
    def tissue_elipsoid(self, G, ax, by, cz):
        remove_cells = []
        for c in self.cells:
            vs_of_c = []
            es_of_c = []
            ps_of_c = list(G.predecessors(("Cell", c.id)))
            for p in ps_of_c:
                es_of_c += list(G.predecessors(p))
            es = list(set(es_of_c))
            for e in es:
                vs_of_c += list(G.predecessors(e))
            vs = list(set(vs_of_c))
            for v in vs:
                argument = (self.vertices[v[1]].r[0] - self.Lxy[0]/2)**2/ax**2 + (self.vertices[v[1]].r[1] - self.Lxy[1]/2)**2/by**2 + (self.vertices[v[1]].r[2] - self.Lxy[2]/2)**2/cz**2
                if argument > 1:
                    remove_cells.append(("Cell", c.id))
                    G.remove_node(("Cell", c.id))
                    c.exist=0
                    break
        self.num_removed_c = len(remove_cells)
        self.remove_nonexistent(G)

    # Create a smooth surface of the tumor 
    def tumor_surface(self, G, DB):
        for c in self.cells:
            if c.exist == 1:
                shared_p = 0
                ps_of_c = list(G.predecessors(("Cell", c.id)))
                for p_of_c in ps_of_c:
                    cs_of_p = list(G.successors(p_of_c))
                    if len(cs_of_p) == 2:
                        shared_p += 1
                if shared_p < 4:
                    G.remove_node(("Cell", c.id))
                    c.exist=0
        for p in self.polygons:
            if p.exist == 1:
                rcount_pc = len(list(G.successors(("Polygon", p.id))))
                if rcount_pc < 2:
                    G.remove_node(("Polygon", p.id))
                    p.exist=0
        self.remove_nonexistent(G)
        
        for v in self.vertices:
            if v.exist == 1:
                es_of_v = list(G.successors(("Vertex", v.id)))
                if len(es_of_v) == 2:
                    # pattern matching
                    eL0eR0 = es_of_v
                    eL0 = eL0eR0[0]
                    eR0 = eL0eR0[1]
                    ps1 = list(G.successors(eL0))
                    ps2 = list(G.successors(eR0))
                    adjacent_p = list(set(ps1) & set(ps2))[0]
                    
                    vL = ("Vertex", v.id)
                    eL = eR0
                    esL = []
                    vsL = [vL]
                    for i in range(100):
                        es_of_vL = list(G.successors(vL))
                        for e_of_vL in es_of_vL:
                            if e_of_vL != eL:
                                e_next = e_of_vL
                                break
                        vs_e_next = list(G.predecessors(e_next))
                        for v_e_next in vs_e_next:
                            if v_e_next != vL:
                                v_next = v_e_next
                                break
                        eL = e_next # next edge
                        esL.append(eL)
                        vL = v_next # next vertex
                        if len(list(G.successors(vL))) > 2:
                            vI = vL
                            break
                        else:
                            vsL.append(vL)
                    
                    vR = ("Vertex", v.id)
                    eR = eL0
                    esR = []
                    vsR = []
                    for i in range(100):
                        es_of_vR = list(G.successors(vR))
                        for e_of_vR in es_of_vR:
                            if e_of_vR != eR:
                                e_next = e_of_vR
                                break
                        vs_e_next = list(G.predecessors(e_next))
                        for v_e_next in vs_e_next:
                            if v_e_next != vR:
                                v_next = v_e_next
                                break
                        eR = e_next # next edge
                        esR.append(eR)
                        vR = v_next # next vertex
                        if len(list(G.successors(vR))) > 2:
                            vF = vR
                            break
                        else:
                            vsR.append(vR)
                    
                    remove_edges = esL + esR
                    for re in remove_edges:
                        G.remove_node(re)
                        self.edges[re[1]].exist=0
                    remove_verts = vsL + vsR
                    for rv in remove_verts:
                        G.remove_node(rv)
                        self.vertices[rv[1]].exist=0
                    
                    free_nodes = list(G.successors("Free"))
                    free_edge_nodes = 0
                    for free_node in free_nodes:
                        if free_node[0] == "Edge":
                            free_edge_nodes += 1
                    if free_edge_nodes < 1:
                        DB.add_edg(G, self)
                    free_nodes = list(G.successors("Free"))
                    for free_node in free_nodes:
                        if free_node[0] == "Edge":
                            surf_e = free_node
                            break
                    
                    # graph transformations
                    G.remove_edge("Free", surf_e)
                    G.add_edge(vI, surf_e, ("IS_PART_OF", 1))
                    G.add_edge(vF, surf_e, ("IS_PART_OF", -1))
                    adjacent_e = Q.hos_of_lo_hho(G, vI, adjacent_p)[0]
                    G.add_edge(surf_e, adjacent_p, ("IS_PART_OF", Q.qSIGNe(G, s=surf_e, t=adjacent_p, v=vI, e=adjacent_e)))
                    
                    # update edges and polygons
                    vs_surf_e = list(G.predecessors(surf_e))
                    vs_info = []
                    for vse in vs_surf_e:
                        r_ve = list(G.get_edge_data(vse, surf_e))[0][1]
                        vs_info.append([vse[1], r_ve])
                    edg = [surf_e[1], vs_info]
                    for i in range(len(edg[1])):
                        if edg[1][i][1] == -1:
                            tail=edg[1][i][0]
                        else:
                            head=edg[1][i][0]
                    self.edges[surf_e[1]].e = [head, tail]
                    es_adjacent_p = list(G.predecessors(adjacent_p))
                    es_info = []
                    for e in es_adjacent_p:
                        r_ep = list(G.get_edge_data(e, adjacent_p))[0][1]
                        es_info.append([e[1], r_ep])
                    poly = [adjacent_p[1], es_info]
                    loe = []
                    loo = []
                    for i in range(len(poly[1])):
                        loe.append(poly[1][i][0])
                        loo.append(poly[1][i][1])
                    self.polygons[poly[0]].loe = loe
                    self.polygons[poly[0]].loo = loo
                    # update
                    self.edges[surf_e[1]].exist = 1
                    self.edges[surf_e[1]].reset_prop(self)
        
        for c in self.cells:
            if c.exist == 1:
                es_of_c = []
                ps_of_c = list(G.predecessors(("Cell", c.id)))
                for p_of_c in ps_of_c:
                    es_of_c += list(G.predecessors(p_of_c))
                es_cell = list(set(es_of_c))
                surface_cell = False
                surface_es = []
                surface_vvs = []
                for e in es_cell:
                    rcount_ep = len(list(G.successors(e)))
                    if rcount_ep < 3:
                        surface_cell = True
                        surface_es.append(e)
                        surface_vvs.append(list(G.predecessors(e)))
                if surface_cell==True:
                    free_nodes = list(G.successors("Free"))
                    free_polygon_nodes = 0
                    for free_node in free_nodes:
                        if free_node[0] == "Polygon":
                            free_polygon_nodes += 1
                    if free_polygon_nodes < 1:
                        DB.add_polygon(G, self)
                    # pattern matching
                    free_nodes = list(G.successors("Free"))
                    for free_node in free_nodes:
                        if free_node[0] == "Polygon":
                            surf_p = free_node
                            break
                    
                    ordered_es = [surface_es[0]]
                    first_p = Q.hos_of_lo_hho(G, ordered_es[0], ("Cell", c.id))[0]
                    ordered_vs = []
                    next_v = surface_vvs[0][1]
                    for i in range(len(surface_es)-1):
                        ordered_vs.append(next_v)
                        for j in range(len(surface_es)):
                            if surface_es[j] not in ordered_es:
                                if surface_vvs[j][0] == next_v:
                                    ordered_es.append(surface_es[j])
                                    next_v = surface_vvs[j][1]
                                    break
                                elif surface_vvs[j][1] == next_v:
                                    ordered_es.append(surface_es[j])
                                    next_v = surface_vvs[j][0]
                                    break
                    
                    # graph transformations
                    G.remove_edge("Free", surf_p)
                    for i in range(len(ordered_es)):
                        if i == 0:
                            G.add_edge(ordered_es[0], surf_p, ("IS_PART_OF", 1))
                        else:
                            G.add_edge(ordered_es[i], surf_p, ("IS_PART_OF", Q.qSIGNe(G, s=ordered_es[i], t=surf_p, v=ordered_vs[i-1], e=ordered_es[i-1])))
                    G.add_edge(surf_p, ("Cell", c.id), ("IS_PART_OF", Q.qSIGNe(G, s=surf_p, t=("Cell", c.id), v=ordered_es[0], e=first_p)))
                    
                    # update polygons and cells
                    es_surf_p = list(G.predecessors(surf_p))
                    es_info = []
                    for e in es_surf_p:
                        r_ep = list(G.get_edge_data(e, surf_p))[0][1]
                        es_info.append([e[1], r_ep])
                    poly = [surf_p[1], es_info]
                    loe = []
                    loo = []
                    for i in range(len(poly[1])):
                        loe.append(poly[1][i][0])
                        loo.append(poly[1][i][1])
                    self.polygons[surf_p[1]].loe = loe
                    self.polygons[surf_p[1]].loo = loo
                    ps_of_c = list(G.predecessors(("Cell", c.id)))
                    ps_info = []
                    for p in ps_of_c:
                        r_pc = list(G.get_edge_data(p, ("Cell", c.id)))[0][1]
                        ps_info.append([p[1], r_pc])
                    cqs = [c.id, ps_info]
                    lop=[]
                    loo=[]
                    for i in range(len(cqs[1])):
                        lop.append(cqs[1][i][0])
                        loo.append(cqs[1][i][1])
                    c.lop=lop
                    c.loo=loo
                    # update
                    self.polygons[surf_p[1]].exist = 1
                    self.polygons[surf_p[1]].reset_prop(self)
    
    # Growth of cells in a system with periodic boundaries
    def growth_PBC(self):
        box_volume = sum([c0.V0 for c0 in self.cells]) + self.num_removed_c
        L_prev = box_volume**(1.0/3.0)
        for c in self.cells:
            if c.growing == 1:
                c.V0 += c.growth_rate*self.h
        box_volume = sum([c0.V0 for c0 in self.cells]) + self.num_removed_c
        L_curr = box_volume**(1.0/3.0)
        self.Lxy=[L_curr, L_curr, L_curr]
        for v in self.vertices:
            v.r[0] *= L_curr/L_prev
            v.r[1] *= L_curr/L_prev
            v.r[2] *= L_curr/L_prev
    
    # Growth of cells in a system with free boundaries
    def growth_FREE(self):
        for c in self.cells:
            if c.growing == 1:
                c.V0 += c.growth_rate*self.h
    
    # Calculate the tumor's center
    def tumor_center(self):
        num_vs=0
        t_center=[0,0,0]
        for v in self.vertices:
            if v.exist == 1:
                num_vs += 1
                for j in range(3):
                    t_center[j] += v.r[j]
        return [t_center[0]/num_vs, t_center[1]/num_vs, t_center[2]/num_vs]
    
    # Create a cubic simulation box with side size 'l'
    def fix_box(self, l):
        tumor_center = self.tumor_center()
        L_new = l
        self.Lxy = [L_new, L_new, L_new]
        box_center = [L_new/2, L_new/2, L_new/2]
        tumor_center_move = [box_center[0]-tumor_center[0], box_center[1]-tumor_center[1], box_center[2]-tumor_center[2]]
        for v in self.vertices:
            if v.exist == 1:
                for j in range(3):
                    v.r[j] += tumor_center_move[j]
    
    # Define the layers of the tumor
    def layer(self, G, n):
        for c in self.cells:
            if c.exist == 1 and c.boundary_cell == 0:
                ps_of_c = list(G.predecessors(("Cell", c.id)))
                for p_of_c in ps_of_c:
                    cs_of_p = list(G.successors(p_of_c))
                    for c_of_p in cs_of_p:
                        if self.cells[c_of_p[1]].boundary_cell == n-1:
                            c.boundary_cell = n
                            break
                    if c.boundary_cell == n:
                        break

    # Set the tension between the live layers and the necrotic layers
    def adjust_tension(self, G, m):
        for p in self.polygons:
            p.g0 = 1
        for c in self.cells:
            if c.boundary_cell == m:
                ps_of_c = list(G.predecessors(("Cell", c.id)))
                for p_of_c in ps_of_c:
                    cs_of_p = list(G.successors(p_of_c))
                    for c_of_p in cs_of_p:
                        if self.cells[c_of_p[1]].boundary_cell == m+1:
                            self.polygons[p_of_c[1]].g0 = self.g0_layer
                            break
    
    # Define the boundary of the tumor
    def boundary(self, G):
        for c in self.cells:
            if c.exist == 1:
                c.boundary_cell = 0
                ps_of_c = list(G.predecessors(("Cell", c.id)))
                for p_of_c in ps_of_c:
                    self.polygons[p_of_c[1]].boundary_polygon = False
                    cs_of_p = list(G.successors(p_of_c))
                    if len(cs_of_p) == 1:
                        c.boundary_cell = 1
                        self.polygons[p_of_c[1]].boundary_polygon = True

    # Set the tension on the tumor's boundary polygons
    def boundary_polygon_tension(self):
        for p in self.polygons:
            if p.exist == 1 and p.boundary_polygon == True:
                p.g0 = self.g0_boundary

    # Initialize cell growth for live cells
    def initialize_cell_growth(self):
        for c in self.cells:
            if c.exist == 1:
                if 0 < c.boundary_cell <= self.live_layers:
                    rand = random.random() # 0-1
                    division_prob = self.division_rate*self.h/(2**(c.boundary_cell-1))
                    if rand < division_prob:
                        c.growing = 1

    # Fluctuations on the live cells' edges
    def live_noise(self, G):
        for e in self.edges:
            e.sigma = 0
        for c in self.cells:
            if c.exist == 1:
                if 0 < c.boundary_cell <= self.live_layers:
                    ps_of_c = list(G.predecessors(("Cell", c.id)))
                    for p_of_c in ps_of_c:
                        es_of_p = list(G.predecessors(p_of_c))
                        for e_of_p in es_of_p:
                            self.edges[e_of_p[1]].sigma = self.sig

    # This function simulates the vertex model
    def simulate(self, G, DB, outFolder):
        # set output folder
        self.outFolder=outFolder
        # set initial conditions for edge tensions
        for e in self.edges:
            if e.exist==1:
                e.reset_prop(self)
        
        # Initialize the simulation box
        self.fix_box(100)

        # Main simulation
        while (self.Time<self.Tmax):
            if self.topological_transition_happened:
                if self.Nc - self.num_removed_c >= self.Nmax:
                    break
                # define the tumor's boundary
                self.boundary(G)
                # define the layers of the tumor
                for l in range(self.live_layers):
                    self.layer(G, l+2)
                # adjust tension between live and necrotic layers
                self.adjust_tension(G, self.live_layers)
                # fluctuations on the live cells' edges
                self.live_noise(G)
                # adjust tension on tumor's boundary polygons
                self.boundary_polygon_tension()
            
            # initialize cell growth for live cells
            self.initialize_cell_growth()
            
            # growth of cells in a system with free boundaries
            self.growth_FREE()
            # equation of motion
            self.solver()
            # T1 transitions and divisions
            success = self.topological_transitions(G, DB)
            if not success:
                return False, 'out_{nsim}_{nmax}_{sigma}_{drate}_{tens}_{lays}_{g0b}_{fnum}.vt3d'.format(nsim=self.Nsim, nmax=self.Nmax, sigma=self.sig, drate=self.division_rate, tens=self.g0_layer, lays=self.live_layers, g0b=self.g0_boundary, fnum=self.fileCount), 'out_{nsim}_{nmax}_{sigma}_{drate}_{tens}_{lays}_{g0b}_0.vt3d'.format(nsim=self.Nsim, nmax=self.Nmax, sigma=self.sig, drate=self.division_rate, tens=self.g0_layer, lays=self.live_layers, g0b=self.g0_boundary), self.fileCount, self.Time, self.timeCount, self.nET, self.nTE, self.nD

        # Save the final output
        with open(self.outFolder+"/out_{nsim}_{nmax}_{sigma}_{drate}_{tens}_{lays}_{g0b}.txt".format(nsim=self.Nsim, nmax=self.Nmax, sigma=self.sig, drate=self.division_rate, tens=self.g0_layer, lays=self.live_layers, g0b=self.g0_boundary), 'a') as filebv:
            filebv.write(str(round(self.Time, 6))+"\t"+str(self.Nc - self.num_removed_c)+"\t"+str(self.ET_count)+"\t"+str(self.TE_count)+"\t"+str(self.D_count)+'\n')
        filebv.close()
        self.fileCount+=1
        self.out_vt3d(self.outFolder, self.fileCount)

        return True, '', '', self.fileCount, self.Time, self.timeCount, self.nET, self.nTE, self.nD
