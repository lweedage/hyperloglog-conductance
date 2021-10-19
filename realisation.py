import networkx as nx
import hyperloglog as hll

def ball_realisation(G, number_of_iterations, number_of_nodes):
    real_ball_around_node = list(list(set() for i in range(number_of_iterations + 2)) for i in range(number_of_nodes))
    real_edgeball_around_node = list(list(set() for i in range(number_of_iterations + 2)) for i in range(number_of_nodes))
    real_directed_edgeball_around_node = list(list(set() for i in range(number_of_iterations + 2)) for i in range(number_of_nodes))


    for t in range(number_of_iterations + 1):
        for node in range(number_of_nodes):
            if t == 0:
                real_ball_around_node[node][t].add(str(node))
            else:
                real_ball_around_node[node][t].add(str(node))
                for k in real_ball_around_node[node][t-1]:
                    for neighbour in G.neighbors(int(k)):
                        real_ball_around_node[node][t].add(str(neighbour))
                    for edge in G.edges(int(k)):
                        new_edge = edge
                        if edge[0] > edge[1]:
                            new_edge = (edge[1], edge[0])
                        real_edgeball_around_node[node][t-1].add(new_edge)
                        real_directed_edgeball_around_node[node][t-1].add(edge)
    return real_edgeball_around_node, real_directed_edgeball_around_node, real_ball_around_node

def find_real_wedges(G, number_of_iterations, number_of_nodes):
    real_wedges = list(list(set() for i in range(number_of_iterations + 1)) for i in range(number_of_nodes))

    for t in range(number_of_iterations + 1):
        for node in range(number_of_nodes):
            for node1 in G.neighbors(int(node)):
                for node2 in G.neighbors(int(node)):
                    hash = hll.wedge_hash([node, node1, node2])
                    real_wedges[node][t].add(hash)
                    if t > 0:
                        real_wedges[node][t] |= real_wedges[node1][t-1]
                        real_wedges[node][t] |= real_wedges[node2][t-1]
    for t in range(number_of_iterations + 1):
        for node in range(number_of_nodes):
            real_wedges[node][t] = len(real_wedges[node][t])
    return real_wedges

