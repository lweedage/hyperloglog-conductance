import networkx as nx
import hyperloglog as hll
import time
import numpy as np



def Hyper_Ball(number_of_nodes, b, G, Triangles = False, Edges = True, Wedges = False):
    # initializing the hyperloglog counters: length of the hash function, the number of nodes that we will have in the graph,
    # p^b registers for the Hyperloglog counter and the hyperloglog counter itself.
    size_bepalen = 0
    m = 2**b
    alpha = hll.find_alpha(m)
    degree_list = [nx.degree(G, node) for node in range(number_of_nodes)]
    number_of_iterations = 3
    # ----------------------------initialization--------------------------------------------------
    edge_counter = np.zeros((2, number_of_nodes, 2**b), dtype=np.int8)
    edge_directed_counter = np.zeros((2, number_of_nodes, 2**b), dtype=np.int8)
    edgeball_around_node = np.zeros((number_of_nodes, number_of_iterations + 1), dtype='f')
    directed_edgeball_around_node = np.zeros((number_of_nodes, number_of_iterations + 1), dtype='f')

    triangle_counter = np.zeros((2, number_of_nodes, 2**b), dtype=np.int8)
    triangles_around_node = np.zeros((number_of_nodes, number_of_iterations + 1), dtype='f')

    wedge_counter = np.zeros((2, number_of_nodes, 2 ** b), dtype=np.int8)
    wedges_around_node = np.zeros((number_of_nodes, number_of_iterations + 1), dtype='f')

    new, old = 0, 1


    # ----------------------------start of the program--------------------------------------------
    # ----------------------------initialization--------------------------------------------------

    start = time.time()
    hll.pre_compute_log(b)
    for node in range(number_of_nodes):
        if node % 1000 == 0:
            print(node)
        for edge in nx.edges(G, node):
            if Edges:
                hll.add(edge_counter[new, node], hll.edge_hash_undirected(edge), b)
                hll.add(edge_directed_counter[new, node], hll.edge_hash_directed(edge), b)
            i = edge[1]
            if Triangles or Wedges:
                for j in nx.neighbors(G, node):
                    hll.add(wedge_counter[new, node], hll.wedge_hash((node, i, j)), b)
                    if j > node:
                         if i in nx.neighbors(G, j) and Triangles:
                            if i > j:
                                assert len({i, j, node}) == 3
                                hll.add(triangle_counter[new, node], hll.triangle_hash((node, i, j)), b)
                                hll.add(triangle_counter[new, i], hll.triangle_hash((node, i, j)), b)
                                hll.add(triangle_counter[new, j], hll.triangle_hash((node, i, j)), b)

        start1 = time.time()
        # ------------------------- Calculating sizes ---------------------------------------
        if Edges:
            edgeball_around_node[node, 0] = hll.size(edge_counter[new, node], b, alpha)
            directed_edgeball_around_node[node, 0] = hll.size(edge_directed_counter[new, node], b, alpha)
        if Triangles:
            triangles_around_node[node, 0] = hll.size(triangle_counter[new, node], b, alpha)
        if Wedges:
            wedges_around_node[node, 0] = hll.size(wedge_counter[new, node], b, alpha)
        size_bepalen += time.time() - start1
    print("Initialization done in", time.time() - start, "seconds")
    print("Estimating the size costs", size_bepalen, 'seconds')
    iterations = 0
    # -------------------------------- do 2 iterations ----------------------------------------------
    while iterations < number_of_iterations:
        iterations += 1
        start1 = time.time()
        new, old = old, new
        if Edges:
            edge_counter[new] = edge_counter[old]
            edge_directed_counter[new] = edge_directed_counter[old]

        if Triangles:
            triangle_counter[new] = triangle_counter[old]

        if Wedges:
            wedge_counter[new] = wedge_counter[old]

        for node in range(number_of_nodes):
            for neighbour in G.neighbors(node):
                if Edges:
                    hll.union(edge_counter[old, neighbour], edge_counter[new, node])
                    hll.union(edge_directed_counter[old, neighbour], edge_directed_counter[new, node])

                if Triangles:
                    hll.union(triangle_counter[old, neighbour], triangle_counter[new, node])

                if Wedges:
                    hll.union(wedge_counter[old, neighbour], wedge_counter[new, node])


            if Edges:
                edgeball_around_node[node, iterations] = hll.size(edge_counter[new, node], b, alpha)
                directed_edgeball_around_node[node, iterations] = hll.size(edge_directed_counter[new, node], b, alpha)

            if Triangles:
                triangles_around_node[node, iterations] = hll.size(triangle_counter[new, node], b, alpha)

            if Wedges:
                wedges_around_node[node, iterations] = hll.size(wedge_counter[new, node], b, alpha)


        end = time.time()

        print("Iteration", iterations, "in", end - start1, "seconds")

    return iterations, degree_list, edgeball_around_node, directed_edgeball_around_node, triangles_around_node, wedges_around_node