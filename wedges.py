import functions as func
import realisation
import os
import pickle

def main(Realisations, G, triangles_around_node, real_triangles, wedges_around_node, number_of_nodes, realisation_name, b):
    distance = 1
    if Realisations:
        if os.path.exists(str('Realisations/' + realisation_name + 'wedges.p')):
            real_wedges = pickle.load(
                open(str('Realisations/' + realisation_name + 'wedges.p'), 'rb'))
            real_transitivity = pickle.load(
                open(str('Realisations/' + realisation_name + 'transitivity.p'), 'rb'))
        else:
            real_wedges = realisation.find_real_wedges(G, 3, number_of_nodes)
            real_transitivity = func.find_real_transitivity(real_wedges, real_triangles, number_of_nodes)
            pickle.dump(real_wedges,
                        open(str('Realisations/' + realisation_name + 'wedges.p'), 'wb'),
                        protocol=4)
            pickle.dump(real_transitivity,
                        open(str('Realisations/' + realisation_name + 'transitivity.p'), 'wb'),
                        protocol=4)
    else:
        real_transitivity = []
        real_wedges = []
    transitivity = func.find_transitivity(triangles_around_node, wedges_around_node, number_of_nodes)
    # func.plot_wedges(number_of_nodes, realisation_name, wedges_around_node, real_wedges, Realisations, b, distance)
    if Realisations:
        func.plot_transitivity(number_of_nodes, transitivity, real_transitivity, real_wedges, real_triangles, 2,
                          realisation_name, Realisations, b)
        func.plot_transitivity(number_of_nodes, transitivity, real_transitivity, real_wedges, real_triangles, 1,
                          realisation_name, Realisations, b)
        func.plot_transitivity(number_of_nodes, transitivity, real_transitivity, real_wedges, real_triangles, 0,
                          realisation_name, Realisations, b)
    func.plot_transitivity_Estimate(number_of_nodes, transitivity, realisation_name)

    func.plot_wedges_triangles(wedges_around_node, triangles_around_node, distance, number_of_nodes, realisation_name)
    print("Wedges done")
    return real_transitivity, transitivity, real_wedges