import os
import pickle
import functions as func

def main(Realisations, realisation_name, real_edges_around_node, real_directed_edges, number_of_nodes, edgeball_around_node, directed_edgeball_around_node, b, degreelist, t = 2):
    if Realisations:
        if os.path.exists(str('Realisations/' + realisation_name + 'conductance_edge_based.p')):
            real_conductance = pickle.load(
                open(str('Realisations/' + realisation_name + 'conductance_edge_based.p'), 'rb'))
        else:
            real_conductance = func.find_real_conductance_edgebase(real_edges_around_node, real_directed_edges,
                                                                   number_of_nodes, t)
            pickle.dump(real_conductance,
                        open(str('Realisations/' + realisation_name + 'conductance_edge_based.p'), 'wb'),
                        protocol=4)
    else:
        real_conductance = []

    conductance = func.find_conductance(edgeball_around_node, directed_edgeball_around_node, number_of_nodes)
    if Realisations:
        for distance in [0, 1]:
            func.plot_conductance_real(number_of_nodes, real_conductance, conductance, real_edges_around_node, real_directed_edges,
                              realisation_name,
                              b, distance=distance, Realisations=True)

    func.plot_conductance(number_of_nodes, conductance, realisation_name)
    func.plot_conductance_of_multiple_distances(number_of_nodes, conductance, realisation_name, degreelist)

    print("Conductance done")
    return real_conductance, conductance