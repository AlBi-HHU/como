from optparse import OptionParser
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-p",
                  help="protein cofactor [sf4_sg]",
                  default="sf4_sg")
parser.add_option("-m",
                  help="mineral quotient graph [greigite]",
                  default="greigite")
parser.add_option("-c", "--coordination", help = "consider coordination [false]",
                  action="store_true", dest="coordination", default=False)

(options, args) = parser.parse_args()

G1 = nx.Graph()
G2 = nx.Graph()
l1, l2, c1, c2, id1, id2, e1, e2 = {}, {}, {}, {}, {}, {}, {}, {}

V1 = pd.read_csv(options.p + "_vert.csv", index_col = False, header = 0)
for i, s in V1.iterrows():
    G1.add_node(s.idx)
    l1[s.idx] = s.elem
    c1[s.idx] = s.CN
    id1[s.idx] = s.id

E1 = pd.read_csv(options.p + "_edges.csv", index_col = False, header = 0)
for i, s in E1.iterrows():
    G1.add_edge(s.Source, s.Target)
    e1[(s.Source, s.Target)] = 1.0

V2 = pd.read_csv(options.m + "_vert.csv", index_col = False, header = 0)
for i, s in V2.iterrows():
    G2.add_node(s.idx)
    l2[s.idx] = s.elem
    c2[s.idx] = s.CN
    id2[s.idx] = s.id


E2 = pd.read_csv(options.m + "_edges.csv", index_col = False, header = 0)
for i, s in E2.iterrows():
    G2.add_edge(s.Source, s.Target)
    e2[(s.Source, s.Target)] = 1.0


#print(G1.nodes)
#print(G1.edges)
#print(l1, "\n", c1, "\n", e1)

#print(G2.nodes)
#print(G2.edges)
#print(l2, c1, e2)

from CMCES_ILP_gurobi import CMCES_ILP

#val, code, S = MCES_ILP(G1, l1, c1, e1, G2, l2, c2, e2, None, "default", options.coordination)
#val, code, S = MCES_ILP(G1, l1, c1, e1, G2, l2, c2, e2, None, "GUROBI", options.coordination)
val, code, S = CMCES_ILP(G1, l1, c1, G2, l2, c2, options.coordination)

print(S)

def is_edge(p):
    return p[0] != -1 and p[1] != -1

VS1 = set()
for v in [ p[0][0] for p in S if is_edge(p)]:
    VS1.add(v)
for v in [ p[0][1] for p in S if is_edge(p)]:
    VS1.add(v)
VS2 = set()
for v in [ p[1][0] for p in S if is_edge(p)]:
    VS2.add(v)
for v in [ p[1][1] for p in S if is_edge(p)]:
    VS2.add(v)

inV1 = [ v in VS1 for v in G1.nodes ]
inV2 = [ v in VS2 for v in G2.nodes ]
inE1 = [ e in list(p[0] for p in [ p for p in S if is_edge(p) ]) for e in G1.edges ]
inE2 = [ e in list(p[1] for p in [ p for p in S if is_edge(p) ]) for e in G2.edges ]

#print("inV1: ", inV1)
#print("inV2: ", inV2)
#print("inE1: ", inE1)
#print("inE2: ", inE2)

#print(VS1)
#print(VS2)

def C_str(c):
    if c: return "C"
    return "no_C"

no_shared = len([ p[0] for p in S if is_edge(p) ])
no_nonshared_in_cofactor =len(G1.edges) - no_shared
#print("no shared =", no_shared)
#print("no nonshared in cofactor =", no_nonshared_in_cofactor)
tversky_index = no_shared / (no_shared + no_nonshared_in_cofactor)

fig1 = plt.figure(1)
#fig1.suptitle(options.p + ", consider coordination: " + str(options.coordination))
#fig1.suptitle("Tversky index = " + str(tversky_index))
pos = nx.kamada_kawai_layout(G1)
nx.draw(G1, pos)
nx.draw_networkx_edges(G1, pos, edgelist = [ p[0] for p in S if is_edge(p)], edge_color='r')
nx.draw_networkx_nodes(G1, pos, nodelist = list(VS1), node_color='r')
#nx.draw_networkx_labels(G1, pos, labels=l1, font_color='w')
nx.draw_networkx_labels(G1, pos, labels=id1, font_size=8, font_color='w')
fig1.savefig(options.p + "_" + C_str(options.coordination) + ".pdf")

fig2 = plt.figure(2)
#fig2.suptitle(options.m + ", consider coordination: " + str(options.coordination))
fig2.suptitle("Tversky index = " + str(tversky_index), x = .98, ha = 'right')
pos = nx.kamada_kawai_layout(G2)
nx.draw(G2, pos)
nx.draw_networkx_edges(G2, pos, edgelist = [ p[1] for p in S if is_edge(p)], edge_color='r')
nx.draw_networkx_nodes(G2, pos, nodelist = list(VS2), node_color='r')
#nx.draw_networkx_labels(G2, pos, labels=l2, font_color='w')
nx.draw_networkx_labels(G2, pos, labels=id2, font_size=8, font_color='w')
fig2.savefig(options.m + "_" + C_str(options.coordination) + ".pdf")

# write results
V1['MCES'] = inV1
V2['MCES'] = inV2
E1['MCES'] = inE1
E2['MCES'] = inE2
V1.to_csv(options.p + "_vert_" + C_str(options.coordination) + "_MCES.csv", index = False)
E1.to_csv(options.p + "_edges_" + C_str(options.coordination) + "_MCES.csv", index = False)
V2.to_csv(options.m + "_vert_" + C_str(options.coordination) + "_MCES.csv", index = False)
E2.to_csv(options.m + "_edges_" + C_str(options.coordination) + "_MCES.csv", index = False)
#plt.show()
