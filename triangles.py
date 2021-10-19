import time
import pickle
import functions as func
import os


def main(Realisations, realisation_name, number_of_nodes, triangles_around_node, b, degree_list, G):
    if Realisations:
        if os.path.exists(str('Realisations/' + realisation_name + 'triangles.p')):
            real_triangles_around_node = pickle.load(
                open(str('Realisations/' + realisation_name + 'triangles.p'), 'rb'))
        else:
            real_triangles_around_node = func.triangles(G)
            pickle.dump(real_triangles_around_node,
                        open(str('Realisations/' + realisation_name + 'triangles.p'), 'wb'),
                        protocol=4)
    else:
        real_triangles_around_node = []
    func.plot_triangles_distance_n(number_of_nodes, triangles_around_node, real_triangles_around_node,
                                   realisation_name, b, Realisations)


    func.plot_triangles_of_multiple_distances(number_of_nodes, triangles_around_node, realisation_name, degree_list)

    # print("Triangles done")
    func.plot_vp_chebyshev_triangles(b, realisation_name)
    return real_triangles_around_node