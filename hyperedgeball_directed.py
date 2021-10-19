import networkx as nx
import matplotlib.pyplot as plt
import hyperloglog as hll
import time
import numpy as np
from scipy import stats
import functions as func

def union(old_counter, new_counter, b):
    union_counter = [max(old_counter[i], new_counter[i]) for i in range(2**b)]
    return union_counter


def Hyper_Ball(number_of_nodes, b, G, Triangles = False, Edges = True, Cycles = False, Wedges = False, kCore = False, k = 8):
    # initializing the hyperloglog counters: length of the hash function, the number of nodes that we will have in the graph,
    # p^b registers for the Hyperloglog counter and the hyperloglog counter itself.
    m = 2**b
    alpha = hll.find_alpha(m)
    degree_list = [nx.degree(G, node) for node in range(number_of_nodes)]
    # ----------------------------initialization--------------------------------------------------
    edge_counter = np.zeros((2, number_of_nodes, 2**b), dtype=np.int8)
    edge_directed_counter = np.zeros((2, number_of_nodes, 2**b), dtype=np.int8)
    edgeball_around_node = [[] for i in range(number_of_nodes)]
    directed_edgeball_around_node = [[] for i in range(number_of_nodes)]

    triangle_counter = np.zeros((2, number_of_nodes, 2**b), dtype=np.int8)
    triangles_around_node = [[] for i in range(number_of_nodes)]

    node_counter = np.zeros((2, number_of_nodes, 2**b), dtype=np.int8)
    ball_around_node = [[] for i in range(number_of_nodes)]

    k_core_node_counter = np.zeros((2, number_of_nodes, 2 ** b), dtype=np.int8)
    k_core_ball_around_node = [[] for i in range(number_of_nodes)]

    wedge_counter = np.zeros((2, number_of_nodes, 2 ** b), dtype=np.int8)
    wedges_around_node = [[] for i in range(number_of_nodes)]

    new, old = 0, 1

    # ----------------------------start of the program--------------------------------------------
    # ----------------------------initialization--------------------------------------------------
    start = time.time()

    for node in range(number_of_nodes):
        if kCore:
            if degree_list[node] >= k:
                hll.add(k_core_node_counter[new, node], hll.node_hash(node), b)
        if Cycles:
            hll.add(node_counter[new, node], hll.node_hash(node), b)
        for edge in nx.edges(G, node):
            if Edges:
                hll.add(edge_counter[new, node], hll.edge_hash_directed(edge=(node, edge[0] + edge[1] - node)), b)
                hll.add(edge_directed_counter[new, edge[0]], hll.edge_hash_directed(edge), b)
                hll.add(edge_directed_counter[new, edge[1]], hll.edge_hash_directed(edge), b)
            i = edge[1]
            if Triangles or Wedges:
                for j in nx.neighbors(G, node):
                    if i < j:
                        hll.add(wedge_counter[new, node], hll.wedge_hash((node, i, j)), b)
                    if j > node:
                         if i in nx.neighbors(G, j) and Triangles:
                            if i > j:
                                assert len({i, j, node}) == 3
                                hll.add(triangle_counter[new, node], hll.triangle_hash((node, i, j)), b)
                                hll.add(triangle_counter[new, i], hll.triangle_hash((node, i, j)), b)
                                hll.add(triangle_counter[new, j], hll.triangle_hash((node, i, j)), b)


        # ------------------------- Calculating sizes ---------------------------------------
        if Edges:
            edgeball_around_node[node].append(hll.size(edge_counter[new, node], b, alpha))
            directed_edgeball_around_node[node].append(hll.size(edge_directed_counter[new, node], b, alpha))
        if Triangles:
            triangles_around_node[node].append(hll.size(triangle_counter[new, node], b, alpha))
        if Cycles:
            ball_around_node[node].append(hll.size(node_counter[new, node], b, alpha))
        if kCore:
            k_core_ball_around_node[node].append(hll.size(k_core_node_counter[new, node], b, alpha))
        if Wedges:
            wedges_around_node[node].append(hll.size(wedge_counter[new, node], b, alpha))
    print("Initialization done in", time.time() - start, "seconds")
    iterations = 0

    # -------------------------------- do 2 iterations ----------------------------------------------
    while iterations < 3:
        iterations += 1
        start1 = time.time()
        new, old = old, new
        if Edges:
            edge_counter[new] = edge_counter[old]
            edge_directed_counter[new] = edge_directed_counter[old]

        if Triangles:
            triangle_counter[new] = triangle_counter[old]

        if Cycles:
            node_counter[new] = node_counter[old]

        if kCore:
            k_core_node_counter[new] = k_core_node_counter[old]

        if Wedges:
            wedge_counter[new] = wedge_counter[old]

        for node in range(number_of_nodes):
            for neighbour in nx.neighbors(G, node):
                if Cycles:
                    node_counter[new, node] = union(node_counter[old, neighbour], node_counter[new, node], b)

                if kCore:
                    if degree_list[node] > k and degree_list[neighbour] > k:
                        k_core_node_counter[new, node] = union(k_core_node_counter[old, neighbour], k_core_node_counter[new, node], b)

                if Edges:
                    edge_counter[new, node] = union(edge_counter[old, neighbour], edge_counter[new, node], b)
                    edge_directed_counter[new, node] = union(edge_directed_counter[old, neighbour], edge_directed_counter[new, node], b)


                if Triangles:
                    triangle_counter[new, node] = union(triangle_counter[old, neighbour], triangle_counter[new, node], b)

                if Wedges:
                    wedge_counter[new, node] = union(wedge_counter[old, neighbour], wedge_counter[new, node], b)


            if Edges:
                edgeball_around_node[node].append(hll.size(edge_counter[new, node], b, alpha))
                directed_edgeball_around_node[node].append(hll.size(edge_directed_counter[new, node], b, alpha))

            if Triangles:
                triangles_around_node[node].append(hll.size(triangle_counter[new, node], b, alpha))

            if Cycles:
                ball_around_node[node].append(hll.size(node_counter[new, node], b, alpha))

            if Wedges:
                wedges_around_node[node].append(hll.size(wedge_counter[new, node], b, alpha))

            if kCore:
                k_core_ball_around_node[node].append(hll.size(k_core_node_counter[new, node], b, alpha))


        end = time.time()

        print("Iteration", iterations, "in", end - start1, "seconds")

    return iterations, degree_list, edgeball_around_node, directed_edgeball_around_node, triangles_around_node, ball_around_node, wedges_around_node, k_core_ball_around_node