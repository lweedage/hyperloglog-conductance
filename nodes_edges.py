import pickle
import os
import realisation

def main(Realisations, realisation_name, G, number_of_nodes):
    if Realisations:
        if os.path.exists(str('Realisations/' + realisation_name + 'nodes.p')) and 3 == 2:
            real_edges_around_node = pickle.load(open(str('Realisations/' + realisation_name + 'edges.p'), 'rb'))
            real_directed_edges = pickle.load(
                open(str('Realisations/' + realisation_name + 'edgesdirected.p'), 'rb'))
            real_ball_around_node = pickle.load(open(str('Realisations/' + realisation_name + 'nodes.p'), 'rb'))
        else:
            real_edges_around_node, real_directed_edges, real_ball_around_node = realisation.ball_realisation(G, 3,
                                                                                                              number_of_nodes)
            pickle.dump(real_edges_around_node, open(str('Realisations/' + realisation_name + 'edges.p'), 'wb'),
                        protocol=4)
            pickle.dump(real_directed_edges,
                        open(str('Realisations/' + realisation_name + 'edgesdirected.p'), 'wb'),
                        protocol=4)
            pickle.dump(real_ball_around_node, open(str('Realisations/' + realisation_name + 'nodes.p'), 'wb'),
                        protocol=4)
    else:
        real_edges_around_node = []
        real_directed_edges = []
        real_ball_around_node = []
            # func.plot_edges(number_of_nodes, real_edges_around_node, edgeball_around_node, realisation_name)
            # func.plot_edges_directed(number_of_nodes, real_directed_edges, directed_edgeball_around_node, realisation_name)
    # print("Edges done")
    return real_edges_around_node, real_directed_edges, real_ball_around_node