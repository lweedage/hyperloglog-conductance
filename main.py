import random
import sys
import time
import os
import networkx as nx
import pickle
from networkx import ExceededMaxIterations

import functions as func
import hyperedgeball
import nodes_edges
import prnibble
import triangles
import conductance as cond
import wedges
import matplotlib.pyplot as plt
import hyperedgeball_directed


# ------------------- FILL IN IMPORTANT THINGS --------------------------------
imported = False  # indicates if there is an imported graph or lfr graph
Triangles = True
PRnibble = True
Conductance = True
Realisations = False
Wedges = True
Edges = True
undirected = True
Save = False
format = 'lfr-1'

mu = 0.3
average_degree = 10
minimum_community = 10
maximum_community = 50
maximum_degree = 50

number_of_nodes = 5000
b = 18 # 2**b registers
seed = 0

if number_of_nodes == 5000 and Wedges == True and Realisations == True:
    print('This is not possible - costs too much memory')
    sys.exit()

if format == 'lfr-1':
    average_degree = 10
    minimum_community = 10
    maximum_community = 50
    maximum_degree = 50
    number_of_nodes = 1000
elif format == 'lfr-2':
    average_degree = 20
    minimum_community = 20
    maximum_community = 100
    maximum_degree = 100
    number_of_nodes = 1000
elif format == 'lfr-3':
    average_degree = 10
    minimum_community = 10
    maximum_community = 50
    maximum_degree = 50
    number_of_nodes = 5000

if not imported:
    print("Number of nodes:", number_of_nodes, "Mu:", mu, "b:", b, 'Seed:', seed)

filename = "lfr.edges"
for b in [b]:
    print('b =', b)
    seed = 0

    # -----------------------------------------------------------------------------
    if imported:
        fh = open(str('/home/lotte/Python/Graphs/' + str(filename)), 'rb')
        G = nx.read_edgelist(fh)
        G = nx.convert_node_labels_to_integers(G)
        if undirected:
            G = nx.to_undirected(G)
        number_of_nodes = nx.number_of_nodes(G)
    else:
        try:
            G = func.make_graph(number_of_nodes, mu, average_degree, minimum_community, maximum_community,
                             maximum_degree, seed)
        except ExceededMaxIterations:
            seed += 1
            print('Error in graph building')
            continue

    print("Graph imported. This graph has", number_of_nodes, "nodes and", G.number_of_edges(), "edges")

    func.estimate_time(b, G.number_of_edges(), G.number_of_nodes())

    if imported:
        realisation_name = str(filename[:-6] + 'n=' + str(number_of_nodes) + 'b=' + str(b))
    else:
        realisation_name = str(filename[:-6] + 'n=' + str(number_of_nodes) + '_mu=' + str(
            mu)[0] + str(mu)[2] + '_avgdegree=' + str(average_degree) + '_mincommunity=' + str(
            minimum_community) + '_seed_' + str(seed) + 'b=' + str(b))

    # Edgeball algorithm
    start = time.time()
    if os.path.exists(str('Estimations/' + realisation_name + 'wedges_around_node.p')):
        print('This graph is already stored in memory')
        number_of_iterations = pickle.load(open(str('Estimations/' + realisation_name + 'number_of_iterations.p'), 'rb'))
        degree_list = pickle.load(open(str('Estimations/' + realisation_name + 'degree_list.p'), 'rb'))
        edgeball_around_node = pickle.load(
            open(str('Estimations/' + realisation_name + 'edgeball_around_node.p'), 'rb'))
        directed_edgeball_around_node = pickle.load(
            open(str('Estimations/' + realisation_name + 'directed_edgeball_around_node.p'), 'rb'))
        triangles_around_node = pickle.load(open(str('Estimations/' + realisation_name + 'triangles_around_node.p'), 'rb'))
        wedges_around_node = pickle.load(
            open(str('Estimations/' + realisation_name + 'wedges_around_node.p'), 'rb'))

    else:
        if undirected:
            number_of_iterations, degree_list, edgeball_around_node, directed_edgeball_around_node, triangles_around_node, nodeball_around_node, wedges_around_node = hyperedgeball.Hyper_Ball(
                number_of_nodes, b, G, Triangles, Edges=Edges, Wedges=Wedges)
        else:
            number_of_iterations, degree_list, edgeball_around_node, directed_edgeball_around_node, triangles_around_node, nodeball_around_node, wedges_around_node = hyperedgeball_directed.Hyper_Ball(
                number_of_nodes, b, G, Triangles, Edges=Edges, Wedges=Wedges)
        if Save:
            pickle.dump(number_of_iterations,
                       open(str('Estimations/' + realisation_name + 'number_of_iterations.p'), 'wb'),
                        protocol=4)
            pickle.dump(degree_list,
                        open(str('Estimations/' + realisation_name + 'degree_list.p'), 'wb'),
                        protocol=4)
            if Edges:
                pickle.dump(edgeball_around_node,
                            open(str('Estimations/' + realisation_name + 'edgeball_around_node.p'), 'wb'),
                            protocol=4)
                pickle.dump(directed_edgeball_around_node,
                            open(str('Estimations/' + realisation_name + 'directed_edgeball_around_node.p'), 'wb'),
                            protocol=4)
            if Triangles:
                pickle.dump(triangles_around_node,
                        open(str('Estimations/' + realisation_name + 'triangles_around_node.p'), 'wb'),
                        protocol=4)
            if Wedges:
                pickle.dump(wedges_around_node,
                        open(str('Estimations/' + realisation_name + 'wedges_around_node.p'), 'wb'),
                        protocol=4)


    number_of_edges = sum(degree_list)
    end = time.time()
    # print("The execution of HyperEdgeBall took", end - start, "seconds")


    # Find the conductance and number of triangles and wedges
    if Triangles:
        real_triangles_around_node = triangles.main(Realisations, realisation_name, number_of_nodes, triangles_around_node,
                                                    b, degree_list, G)

    if Edges:
        real_edges_around_node, real_directed_edges, real_ball_around_node = nodes_edges.main(Realisations,
                                                                                              realisation_name, G,
                                                                                              number_of_nodes)

    if Conductance:
        real_conductance, conductance = cond.main(Realisations, realisation_name, real_edges_around_node,
                                                         real_directed_edges, number_of_nodes, edgeball_around_node,
                                                         directed_edgeball_around_node, b, degree_list, t=1)

    if Wedges:
        real_transitivity, transitivity, real_wedges = wedges.main(Realisations, G, triangles_around_node,
                                                                   real_triangles_around_node, wedges_around_node,
                                                                   number_of_nodes, realisation_name, b)



    size = 100


    if PRnibble:
        if os.path.exists(str('Estimations/' + realisation_name + 'prconductanceconductance1.p')):
            continue

        else:
            print('Seed set is conductance distance 1')
            seed_set = sorted(range(number_of_nodes), key=lambda index: conductance[index][0])[:size]
            prconductanceconductance1 = prnibble.main(number_of_nodes, conductance, realisation_name, G, seed_set,
                                                      maximum_community, True)
            conductance1conductance1 = []
            conductance2conductance1 = []
            for node in seed_set:
                conductance1conductance1.append(conductance[node][0])
                conductance2conductance1.append(conductance[node][1])

            print('Seed set is conductance distance 2')
            seed_set = sorted(range(number_of_nodes), key=lambda index: conductance[index][1])[:size]
            prconductanceconductance2 = prnibble.main(number_of_nodes, conductance, realisation_name, G, seed_set,
                                                      maximum_community, True)
            conductance1conductance2 = []
            conductance2conductance2 = []
            for node in seed_set:
                conductance1conductance2.append(conductance[node][0])
                conductance2conductance2.append(conductance[node][1])

            print('Seed set is triangles distance 0')
            seed_set = sorted(range(number_of_nodes), reverse=True,
                              key=lambda index: triangles_around_node[index][0])[:size]

            prconductancetriangles0 = prnibble.main(number_of_nodes, conductance, realisation_name, G, seed_set,
                                                    maximum_community, True)
            conductance1triangles0 = []
            conductance2triangles0 = []
            for node in seed_set:
                conductance1triangles0.append(conductance[node][0])
                conductance2triangles0.append(conductance[node][1])

            print('Seed set is triangles distance 1')
            seed_set = sorted(range(number_of_nodes), reverse=True,
                              key=lambda index: triangles_around_node[index][1])[:size]
            prconductancetriangles1 = prnibble.main(number_of_nodes, conductance, realisation_name, G, seed_set,
                                                    maximum_community, True)
            conductance1triangles1 = []
            conductance2triangles1 = []
            for node in seed_set:
                conductance1triangles1.append(conductance[node][0])
                conductance2triangles1.append(conductance[node][1])

            print('Seed set is degree')
            seed_set = sorted(range(number_of_nodes), reverse=True, key=lambda index: degree_list[index])[:size]
            prconductancedegree = prnibble.main(number_of_nodes, conductance, realisation_name, G, seed_set,
                                                maximum_community, True)
            conductance1degree = []
            conductance2degree = []
            for node in seed_set:
                conductance1degree.append(conductance[node][0])
                conductance2degree.append(conductance[node][1])

            print('Seed set is random')
            seed_set = [random.randint(0, number_of_nodes - 1) for i in range(size)]
            prconductancerandom = prnibble.main(number_of_nodes, conductance, realisation_name, G, seed_set,
                                                maximum_community, True)
            conductance1random = []
            conductance2random = []
            for node in seed_set:
                conductance1random.append(conductance[node][0])
                conductance2random.append(conductance[node][1])

            print('Seed set is transitivity distance 1')
            seed_set = sorted(range(number_of_nodes), reverse=True, key=lambda index: transitivity[index][0])[:size]
            prconductancetransitivity1 = prnibble.main(number_of_nodes, conductance, realisation_name, G, seed_set,
                                                       maximum_community, True)
            conductance1transitivity1 = []
            for node in seed_set:
                conductance1transitivity1.append(conductance[node][0])

            print('Seed set is transitivity distance 2')
            seed_set = sorted(range(number_of_nodes), reverse=True, key=lambda index: transitivity[index][1])[:size]
            prconductancetransitivity2 = prnibble.main(number_of_nodes, conductance, realisation_name, G, seed_set,
                                                       maximum_community, True)
            conductance2transitivity2 = []
            for node in seed_set:
                conductance2transitivity2.append(conductance[node][1])

            pickle.dump(conductance1conductance1,
            open(str('Estimations/' + realisation_name + 'conductance1conductance1.p'), 'wb'), protocol=4)
            pickle.dump(prconductanceconductance1,
            open(str('Estimations/' + realisation_name + 'prconductanceconductance1.p'), 'wb'), protocol=4)
            pickle.dump(conductance2conductance2,
            open(str('Estimations/' + realisation_name + 'conductance2conductance2.p'), 'wb'), protocol=4)
            pickle.dump(prconductanceconductance2,
            open(str('Estimations/' + realisation_name + 'prconductanceconductance2.p'), 'wb'), protocol=4)
            pickle.dump(conductance1triangles0,
            open(str('Estimations/' + realisation_name + 'conductance1triangles0.p'), 'wb'), protocol=4)
            pickle.dump(prconductancetriangles0,
            open(str('Estimations/' + realisation_name + 'prconductancetriangles0.p'), 'wb'), protocol=4)
            pickle.dump(conductance1transitivity1,
            open(str('Estimations/' + realisation_name + 'conductance1transitivity1.p'), 'wb'), protocol=4)
            pickle.dump(prconductancetransitivity1,
            open(str('Estimations/' + realisation_name + 'prconductancetransitivity1.p'), 'wb'), protocol=4)
            pickle.dump(conductance1degree,
            open(str('Estimations/' + realisation_name + 'conductance1degree.p'), 'wb'), protocol=4)
            pickle.dump(prconductancedegree,
            open(str('Estimations/' + realisation_name + 'prconductancedegree.p'), 'wb'), protocol=4)
            pickle.dump(conductance1random,
            open(str('Estimations/' + realisation_name + 'conductance1random.p'), 'wb'), protocol=4)
            pickle.dump(prconductancerandom,
            open(str('Estimations/' + realisation_name + 'prconductancerandom.p'), 'wb'), protocol=4)
            pickle.dump(conductance2triangles1,
            open(str('Estimations/' + realisation_name + 'conductance2triangles1.p'), 'wb'), protocol=4)
            pickle.dump(prconductancetriangles1,
            open(str('Estimations/' + realisation_name + 'prconductancetriangles1.p'), 'wb'), protocol=4)
            pickle.dump(conductance2transitivity2,
            open(str('Estimations/' + realisation_name + 'conductance2transitivity2.p'), 'wb'), protocol=4)
            pickle.dump(prconductancetransitivity2,
            open(str('Estimations/' + realisation_name + 'prconductancetransitivity2.p'), 'wb'), protocol=4)
            pickle.dump(conductance2degree,
            open(str('Estimations/' + realisation_name + 'conductance2degree.p'), 'wb'), protocol=4)
            pickle.dump(conductance2random,
            open(str('Estimations/' + realisation_name + 'conductance2random.p'), 'wb'), protocol=4)

