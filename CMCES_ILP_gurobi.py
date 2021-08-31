# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 20:45:12 2021

@author: klau
"""
import gurobipy as grb
import networkx as nx
import matplotlib.pyplot as plt


# Callback - use lazy constraints for connectivity
def connectivity(model, where):
    if where == grb.GRB.Callback.MIPSOL:
        E1vals = model.cbGetSolution(model._E1)
#        yvals = model.cbGetSolution(model._y)
#        #current = { v for v in model._z1.keys() if z1vals[v] > 0.9 } # current integer solution
#        current = { ep[0] for ep in model._y.keys() if yvals[ep] > 0.9 } # chosen edges in G1
#        print("current via matches:", current)
        current = { e1 for e1 in model._E1.keys() if E1vals[e1] > 0.9 } # chosen edges in G1
#        print("current via vars:   ", current)
#        fig1 = plt.figure(1)
#        pos = nx.kamada_kawai_layout(model._G1)
#        nx.draw(model._G1, pos)
#        nx.draw_networkx_edges(model._G1, pos, edgelist = list(current), edge_color='r')
#        nx.draw_networkx_labels(model._G1, pos, font_size=8, font_color='w')
#        #plt.show()


        #G1_tilde = model._G1.subgraph( { e[0] for e in current } | { e[1] for e in current } )
        G1_tilde = model._G1.copy()
        for e1 in model._G1.edges():
            if e1 not in current: G1_tilde.remove_edge(*e1)

        #print("nodes in G'': ", G1_tilde.nodes())
        #print("edges in G'': ", G1_tilde.edges())

        L1 = nx.line_graph(model._G1)
        L1_tilde = nx.line_graph(G1_tilde)
#        fig2 = plt.figure(2)
#        pos = nx.kamada_kawai_layout(L1_tilde)
#        nx.draw(L1_tilde, pos)
#        nx.draw_networkx_nodes(L1_tilde, pos, nodelist = list(current), node_color='r')
#        nx.draw_networkx_labels(L1_tilde, pos, font_size=8, font_color='w')
#        #plt.show()



        #if L1_tilde.nodes(): print("connected: ", nx.is_connected(L1_tilde))
        if (L1_tilde.number_of_nodes() > 0 and not nx.is_connected(L1_tilde)):
            C = [ c for c in nx.algorithms.components.connected_components(L1_tilde) ]
            for C_l in C:
                 for C_h in C:
                     if C_l != C_h:
                         AC_h =  { v for c in C_h for v in L1.neighbors(c) }
                         L1_prime = L1.copy()
                         L1_prime.remove_edges_from([ e for e in L1.subgraph(C_h | AC_h).edges() ])
                         #print(len(model._Ga.edges()), ":", len(G_prime.edges()))
                         l = next(iter(C_l))
                         R_l = set(nx.dfs_postorder_nodes(L1_prime, l))
                         min_node_sep = AC_h & R_l
                         for l in C_l:
                             for h in C_h:
                                model.cbLazy(model._E1[l] + model._E1[h] - grb.quicksum(model._E1[c] for c in min_node_sep) <= 1)

        E2vals = model.cbGetSolution(model._E2)
#        yvals = model.cbGetSolution(model._y)
#        #current = { v for v in model._z1.keys() if z1vals[v] > 0.9 } # current integer solution
#        current = { ep[0] for ep in model._y.keys() if yvals[ep] > 0.9 } # chosen edges in G2
#        print("current via matches:", current)
        current = { e2 for e2 in model._E2.keys() if E2vals[e2] > 0.9 } # chosen edges in G2
#        print("current via vars:   ", current)
#        fig1 = plt.figure(1)
#        pos = nx.kamada_kawai_layout(model._G1)
#        nx.draw(model._G1, pos)
#        nx.draw_networkx_edges(model._G1, pos, edgelist = list(current), edge_color='r')
#        nx.draw_networkx_labels(model._G1, pos, font_size=8, font_color='w')
#        #plt.show()


        #G1_tilde = model._G1.subgraph( { e[0] for e in current } | { e[1] for e in current } )
        G2_tilde = model._G2.copy()
        for e2 in model._G2.edges():
            if e2 not in current: G2_tilde.remove_edge(*e2)

        #print("nodes in G'': ", G1_tilde.nodes())
        #print("edges in G'': ", G1_tilde.edges())

        L2 = nx.line_graph(model._G2)
        L2_tilde = nx.line_graph(G2_tilde)
#        fig2 = plt.figure(2)
#        pos = nx.kamada_kawai_layout(L1_tilde)
#        nx.draw(L1_tilde, pos)
#        nx.draw_networkx_nodes(L1_tilde, pos, nodelist = list(current), node_color='r')
#        nx.draw_networkx_labels(L1_tilde, pos, font_size=8, font_color='w')
#        #plt.show()



        #if L1_tilde.nodes(): print("connected: ", nx.is_connected(L1_tilde))
        if (L2_tilde.number_of_nodes() > 0 and not nx.is_connected(L2_tilde)):
            C = [ c for c in nx.algorithms.components.connected_components(L2_tilde) ]
            for C_l in C:
                 for C_h in C:
                     if C_l != C_h:
                         AC_h =  { v for c in C_h for v in L2.neighbors(c) }
                         L2_prime = L2.copy()
                         L2_prime.remove_edges_from([ e for e in L2.subgraph(C_h | AC_h).edges() ])
                         #print(len(model._Ga.edges()), ":", len(G_prime.edges()))
                         l = next(iter(C_l))
                         R_l = set(nx.dfs_postorder_nodes(L2_prime, l))
                         min_node_sep = AC_h & R_l
                         for l in C_l:
                             for h in C_h:
                                model.cbLazy(model._E2[l] + model._E2[h] - grb.quicksum(model._E2[c] for c in min_node_sep) <= 1)


def CMCES_ILP(G1, l1, c1, G2, l2, c2, coordination):
    #ILP=pulp.LpProblem("MCES", pulp.LpMinimize)
    ILP = grb.Model(name = "CMCES")
    ILP.Params.lazyConstraints = 1

    # all optimal solutions
    ILP.Params.PoolSearchMode = 2 #
    ILP.Params.PoolGap = 0.0 # nothing suboptimal
    ILP.Params.PoolSolutions = 100000 # 100000 are enough

    x = { (i,j):ILP.addVar(vtype=grb.GRB.BINARY, name="x_{0}_{1}".format(i,j)) for i in G1.nodes for j in G2.nodes }

    y = { (e1, e2):ILP.addVar(vtype=grb.GRB.BINARY, name = "y_{0}_{1}".format(str(e1).replace(" ", ""), str(e2).replace(" ", ""))) for e1 in list(G1.edges) + [-1] for e2 in list(G2.edges) + [-1] }

    E1 = { e1 : ILP.addVar(vtype=grb.GRB.BINARY, name = "e1_{0}".format(str(e1).replace(" ", ""))) for e1 in G1.edges }
    E2 = { e2 : ILP.addVar(vtype=grb.GRB.BINARY, name = "e2_{0}".format(str(e2).replace(" ", ""))) for e2 in G2.edges }

    ILP.update()

    # let the model know some data, for the callback
    ILP._G1 = G1
    ILP._G2 = G2
    ILP._y = y
    ILP._E1 = E1
    ILP._E2 = E2

    ILP.ModelSense = grb.GRB.MINIMIZE
    ILP.setObjective(grb.quicksum(y[e1,-1] for e1 in list(G1.edges)) + grb.quicksum(y[-1,e2] for e2 in list(G2.edges)))

    # for each node: only match it to at most one other one
    for i in G1.nodes:
        ILP.addConstr(grb.quicksum(x[i, j] for j in G2.nodes) <= 1)
    for j in G2.nodes:
        ILP.addConstr(grb.quicksum(x[i, j] for i in G1.nodes) <= 1)

    # for each edge: match it to another edge or to -1
    for e1 in G1.edges:
        ILP.addConstr(grb.quicksum(y[e1, e2] for e2 in list(G2.edges) + [-1]) == 1)
    for e2 in G2.edges:
        ILP.addConstr(grb.quicksum(y[e1, e2] for e1 in list(G1.edges) + [-1]) == 1)

    # normalize the order of an edge (from lower index to higher index)
    def normalize(e):
        if e[0] > e[1]:
            return (e[1], e[0])
        return e

    # link matched edges and matched nodes
    for i in G1.nodes:
        for e2 in G2.edges:
            ILP.addConstr(grb.quicksum(y[normalize(e1), e2] for e1 in G1.edges(i)) <= x[i,e2[0]] + x[i,e2[1]])
    for k in G2.nodes:
        for e1 in G1.edges:
            ILP.addConstr(grb.quicksum(y[e1, normalize(e2)] for e2 in G2.edges(k)) <= x[e1[0],k] + x[e1[1], k])

    # only match what can be matched (could also be solved by not generating the variables)
    for i in G1.nodes:
        for j in G2.nodes:
            if l1[i]!=l2[j]:
                ILP.addConstr(x[i,j] == 0)
            if coordination and (c1[i] != c2[j]):
                ILP.addConstr(x[i,j] == 0)

    # model matched edges in each of the graphs
    for e1 in G1.edges:
        ILP.addConstr(E1[e1] == 1 - y[e1, -1])
    for e2 in G2.edges:
        ILP.addConstr(E2[e2] == 1 - y[-1, e2])
    #ILP.write("check.lp")
    ILP.optimize(connectivity)
    #ILP.optimize()

    if ILP.status != grb.GRB.OPTIMAL:
        print('Optimization ended with status %d' % model.status)
        sys.exit(0)

    print(str(ILP.SolCount) + " optimal solutions exist, here's one")

    S = [ (e1,e2) for e1 in list(G1.edges) + [-1] for e2 in list(G2.edges) + [-1] if y[e1,e2].x > 0.9 ]
    return float(ILP.objVal), 1, S
