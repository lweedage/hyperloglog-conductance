import math
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import unidip.dip as dip
from cycler import cycler

fig_width = 2.809 * 2
fig_height = fig_width/4*3

matplotlib.use('PDF')

matplotlib.rcParams['axes.prop_cycle'] = cycler('color', ['DeepSkyBlue', 'DarkMagenta', 'LightPink', 'Orange', 'LimeGreen', 'OrangeRed'])
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['savefig.format'] = 'pdf'
matplotlib.rcParams['figure.figsize'] = fig_width, fig_height
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['lines.markersize'] = 3
matplotlib.rcParams['figure.autolayout'] = True


def make_graph(number_of_nodes, mu, average_degree, minimum_community, max_community, maximum_degree, seed):
    tau1 = 3
    tau2 = 2

    G = nx.LFR_benchmark_graph(number_of_nodes, tau1, tau2, mu, average_degree=average_degree,
                               max_degree=maximum_degree, min_community=minimum_community, max_community=max_community,
                               seed=seed)

    # Remove the selfloops
    G.remove_edges_from(nx.selfloop_edges(G))
    nx.set_node_attributes(G, {n: ','.join(map(str, G.nodes[n]['community'])) for n in G.nodes()}, 'community')

    return G


def find_beta(b):
    m = 2 ** b
    if m < 17:
        return 1.106
    elif m < 33:
        return 1.070
    elif m < 65:
        return 1.054
    elif m < 129:
        return 1.046
    else:
        return 1.03896


def find_conductance(edgeball_around_node, directed_edgeball_around_node, number_of_nodes):
    conductance = list([] for i in range(number_of_nodes))
    for t in [0, 1, 2]:
        for node in range(number_of_nodes):
            if edgeball_around_node[node][t + 1] == 0:
                conductance[node].append(1)
            else:
                conductance[node].append(2 * edgeball_around_node[node][t + 1] / directed_edgeball_around_node[node][t + 1] - 1)
    return conductance


def find_real_conductance_edgebase(real_edges, real_directed_edges, number_of_nodes, t):
    conductance = list([] for i in range(number_of_nodes))
    for t in [0, 1, 2]:
        for node in range(number_of_nodes):
            if len(real_edges[node][t + 1]) == 0:
                conductance[node].append(1)
            else:
                conductance[node].append(2 * len(real_edges[node][t + 1]) / len(real_directed_edges[node][t + 1]) - 1)
    return conductance


def triangles(G):
    number_of_nodes = G.number_of_nodes()
    triangles_around_node = [[set(), set(), set(), set()] for i in range(number_of_nodes)]
    for node in range(number_of_nodes):
        for neighbor1 in nx.neighbors(G, node):
            for neighbor2 in nx.neighbors(G, node):
                if neighbor1 in set(nx.neighbors(G, neighbor2)):
                    a = min(node, neighbor1, neighbor2)
                    b = max(node, neighbor1, neighbor2)
                    c = node + neighbor1 + neighbor2 - a - b
                    triangles_around_node[node][0].add((a, c, b))
    for iteration in [1, 2, 3]:
        for node in range(number_of_nodes):
            triangles_around_node[node][iteration] = triangles_around_node[node][iteration - 1].copy()
            for neighbor in nx.neighbors(G, node):
                triangles_around_node[node][iteration] |= triangles_around_node[neighbor][iteration - 1]
    return triangles_around_node

def plot_triangles_distance_n(number_of_nodes, triangles_around_node, real_triangles_around_node, realisation_name,
                              b, Realisations=True):
    fig, ax = plt.subplots()
    lijst2 = [[] for i in range(3)]
    lowerbound = [[] for i in range(3)]
    upperbound = [[] for i in range(3)]
    lowerboundvp = [[] for i in range(3)]
    upperboundvp = [[] for i in range(3)]

    for n in [0, 1, 2]:
        triangles_around_node_of_distance_n = []
        real_triangles_around_node_of_distance_n = []

        for i in range(number_of_nodes):
            triangles_around_node_of_distance_n.append(triangles_around_node[i][n])
            if Realisations:
                real_triangles_around_node_of_distance_n.append(len(real_triangles_around_node[i][n]))

        if Realisations:
            sortedindex_real_triangles_n = sorted(range(number_of_nodes),
                                                  key=lambda index: real_triangles_around_node_of_distance_n[index])

            # Plotting
            lijst1 = []
            delta2 = 0.0005
            etam = find_beta(b) / math.sqrt(2 ** b)
            eta = etam + delta2

            for j in sortedindex_real_triangles_n:
                lijst1.append(triangles_around_node_of_distance_n[j])
                p4 = real_triangles_around_node_of_distance_n[j] * eta / math.sqrt(0.05)
                p5 = 2 / 3 * real_triangles_around_node_of_distance_n[j] * eta / math.sqrt(0.05)
                lowerbound[n].append(real_triangles_around_node_of_distance_n[j] - p4)
                upperbound[n].append(real_triangles_around_node_of_distance_n[j] + p4)
                lowerboundvp[n].append(real_triangles_around_node_of_distance_n[j] - p5)
                upperboundvp[n].append(real_triangles_around_node_of_distance_n[j] + p5)

            lijst2[n] = sorted(real_triangles_around_node_of_distance_n)

        else:
            lijst1 = sorted(triangles_around_node_of_distance_n)

        # ------------------------------------ plot triangles --------------------------------------
        ax.plot(range(number_of_nodes), lijst1, 'o', label=f'$S_{n + 1}(v)$')
    First = True
    if Realisations:
        for i in range(3):
            if First:
                ax.plot(range(number_of_nodes), lowerbound[i], '-', color = 'LimeGreen',
                        label=str("Chebyshev"))
                ax.plot(range(number_of_nodes), upperbound[i], '-', color = 'LimeGreen')
                ax.plot(range(number_of_nodes), lowerboundvp[i], '-', color = 'Orange',
                        label=str("VP"))
                ax.plot(range(number_of_nodes), upperboundvp[i], '-', color = 'Orange')
                ax.plot(range(number_of_nodes), lijst2[i], '-', color = 'OrangeRed', label='Realisation')

                First = False
            else:
                ax.plot(range(number_of_nodes), lowerbound[i], '-', color = 'LimeGreen')
                ax.plot(range(number_of_nodes), upperbound[i], '-', color = 'LimeGreen')
                ax.plot(range(number_of_nodes), lowerboundvp[i], '-', color = 'Orange')
                ax.plot(range(number_of_nodes), upperboundvp[i], '-', color = 'Orange')
                ax.plot(range(number_of_nodes), lijst2[i], color = 'OrangeRed')

    ax.set(xlabel='Node $v$', ylabel='$|\Delta_r(v)|$')
    #plt.title('Number of triangles in $S_1(v)$, $S_2(v)$ and $S_3(v)$')#, fontsize = 12)
    name = str('Pictures/' + realisation_name + '_triangles_distance_1_and_2_and_3.pdf')
    ax.set_yscale('log')

    ax.legend()
    plt.savefig(name)
    plt.show()


def plot_triangles_of_multiple_distances(number_of_nodes, triangles_around_node, realisation_name, degreelist):
    sortedindex = sorted(range(number_of_nodes),
                         key=lambda index: degreelist[index])
    fig, ax = plt.subplots()
    First = True
    for node in sortedindex[:10]:
        lijst = []
        for distance in range(4):
            lijst.append(triangles_around_node[node][distance])
            if First and 3 == 2:
                data = np.msort([triangles_around_node[node][distance] for node in range(number_of_nodes)])
                print('DIP for triangles in $S_', distance, '(v)$ gives us', dip.diptst(data))
        if First:
            ax.plot([i for i in range(4)], lijst, ':', color='lightblue')
            ax.plot([i for i in range(4)], lijst, 'o',markersize = 6, color = 'DeepSkyBlue', label=str('Low degree nodes'))

            First = False
        else:
            ax.plot([i for i in range(4)], lijst, ':', color='lightblue')
            ax.plot([i for i in range(4)], lijst, 'o',markersize = 6, color = 'DeepSkyBlue')

    First = True
    for node in sortedindex[-10:]:
        lijst = []
        for distance in range(4):
            lijst.append(triangles_around_node[node][distance])
        if First:
            ax.plot([i for i in range(4)], lijst, ':', color='pink')
            ax.plot([i for i in range(4)], lijst, 'o', markersize = 6, color = 'OrangeRed',  label=str('High degree nodes'))
            First = False
        else:
            ax.plot([i for i in range(4)], lijst, ':', color='pink')
            ax.plot([i for i in range(4)], lijst, 'o', markersize = 6,  color = 'OrangeRed')

    ax.set(xlabel='Radius $r$', ylabel='$|\hat{\Delta}_r(v)|$')#, title=str('com-Amazon: Triangles in different radii $r$'))
    name = str('Pictures/' + realisation_name + '_triangles_of_all_distances.pdf')
    major_ticks = np.arange(0, 4, 1)
    ax.set_xticks(major_ticks)
    ax.legend()
    plt.savefig(name)
    plt.show()

def plot_conductance_real(number_of_nodes, real_conductance, conductance, real_edges, real_directed_edges,
                          realisation_name,
                          b, distance, Realisations=True):
    fig, ax = plt.subplots()
    delta1 = 0.00005
    delta2 = 0.0005
    etam = find_beta(b) / math.sqrt(2 ** b)
    eta3 = etam + delta2

    sortedindex_real_conductance = sorted(range(number_of_nodes), key=lambda index: real_conductance[index][distance])

    sorted_conductance = []
    sorted_real_conductance = []

    for i in sortedindex_real_conductance[1:10]:
        sorted_real_conductance.append(real_conductance[i][distance])
        sorted_conductance.append(conductance[i][distance])

    # Plotting
    lijst1 = []
    lijst2 = []

    chebyshevlowervariance = []
    chebyshevuppervariance = []
    vplowervariance = []
    vpuppervariance = []
    vpbounds = 0
    chebbounds = 0

    for j in sortedindex_real_conductance:
        lijst1.append(conductance[j][distance])
        lijst2.append(real_conductance[j][distance])
        p = eta3 / math.sqrt(0.05) * math.sqrt(
            len(real_edges[j][distance]) ** 2 + len(real_directed_edges[j][distance]) ** 2)
        epsilon = p / len(real_edges[j][distance]) - delta1
        delta = p / len(real_directed_edges[j][distance]) - delta1
        chebyshevlowervariance.append(real_conductance[j][distance] * (1 - epsilon) / (1 + delta))
        chebyshevuppervariance.append(real_conductance[j][distance] * (1 + epsilon) / (1 - delta))
        vp = math.sqrt(4 / 9) * eta3 / math.sqrt(0.05) * math.sqrt(
            len(real_edges[j][distance]) ** 2 + len(real_directed_edges[j][distance]) ** 2)
        epsilonvp = vp / len(real_edges[j][distance]) - delta1
        deltavp = vp / len(real_directed_edges[j][distance]) - delta1
        vplowervariance.append(real_conductance[j][distance] * (1 - epsilonvp) / (1 + deltavp))
        vpuppervariance.append(real_conductance[j][distance] * (1 + epsilonvp) / (1 - deltavp))
        vpbounds += 1 - (1 + epsilonvp) / (1 - deltavp)
        chebbounds += 1 - (1 + epsilon) / (1 - delta)

    print('VP:', vpbounds / number_of_nodes)
    print('Cheb:', chebbounds / number_of_nodes)
    # ------------------------------------ plot conductance --------------------------------------
    #plt.title('Conductance')

    ax.plot(range(number_of_nodes), lijst1, 'o', label='Estimate')
    ax.plot(range(number_of_nodes), lijst2, color = 'OrangeRed', label='Realisation')
    ax.plot(range(number_of_nodes), chebyshevlowervariance, '-', color = 'LimeGreen',
            label=str('Chebyshev'))
    ax.plot(range(number_of_nodes), chebyshevuppervariance, '-', color =
            'LimeGreen')
    ax.plot(range(number_of_nodes), vplowervariance, '-', color = 'Orange',
            label=str('VP'))
    ax.plot(range(number_of_nodes), vpuppervariance, '-', color = 'Orange')
    #plt.title('Conductance in $S_' + str(distance + 1) + '(v)$')
    plt.ylim(0.45, 1.05)

    name = str(
        "Pictures/" + realisation_name + '_conductance_with_Chebyshev_and_VP_distance_' + str(distance + 1) + '.pdf')

    plt.xlabel('Node $v$')
    plt.ylabel('Conductance $\phi(S_' + str(distance + 1) +
               '(v))$')
    ax.legend()
    plt.savefig(name)
    plt.show()

    if Realisations:
        fig, ax = plt.subplots()

        differences = [lijst2[i] - lijst1[i] for i in range(number_of_nodes)]
        average_difference = sum(differences) / number_of_nodes
        kwargs = dict(alpha=0.5, bins=50, density=False, stacked=True)

        # Plot
        plt.hist(differences, **kwargs)
        #plt.gca().set(title='Error in estimated conductance in $S_' + str(distance + 1) + '(v)$')
        plt.xlabel('Realisation - Estimate')
        ax.grid()
        name = str(
            "Pictures/" + realisation_name + '_conductance_difference_histogram_distance_' + str(distance + 1) + '.pdf')
        plt.savefig(name)
        plt.show()


def plot_conductance(number_of_nodes, conductance, realisation_name):
    fig, ax = plt.subplots()

    # ------------------------------------ plot conductance --------------------------------------
    for distance in [1, 2]:
        ax.plot(range(number_of_nodes), sorted([conductance[i][distance] for i in range(number_of_nodes)]), 'o',
            label=str('$S_' + str(distance) + '(v)$'))
    #plt.title('Estimated conductance $\hat{\phi}(S_r(v))$')
    name = str("Pictures/" + realisation_name + '_conductance.pdf')

    plt.xlabel('Node $v$')
    plt.ylabel('Conductance $\hat{\phi}(S_r(v))$')
    ax.legend()
    plt.savefig(name)
    plt.show()


def plot_transitivity_Estimate(number_of_nodes, transitivity, realisation_name):
    fig, ax = plt.subplots()

    #plt.title('Estimated transitivity $\hat{t}(S_r(v))$')
    for distance in [0, 1, 2]:
        ax.plot(range(number_of_nodes), sorted([transitivity[i][distance] for i in range(number_of_nodes)]), 'o',
                label=str('$S_' + str(distance + 1) + '(v)$'))
        data = np.msort([transitivity[node][distance] for node in range(number_of_nodes)])
        print('DIP for transitivity in S_', distance + 1 , '(v) gives us', dip.diptst(data))
    name = str("Pictures/" + realisation_name + '_transitivity.pdf')

    plt.xlabel('Node $v$')
    plt.ylabel(str('$\hat{t}(S_r(v))$'))
    ax.legend()
    plt.savefig(name)
    plt.show()


def find_real_cycles(number_of_nodes, real_edges_around_node, real_directed_edges, real_ball_around_node, r):
    cycles = [0] * number_of_nodes
    for i in range(number_of_nodes):
        cycles[i] = len(real_directed_edges[i][r]) - len(real_edges_around_node[i][r]) - len(
            real_ball_around_node[i][r]) + 1
    return cycles


def find_cycles(number_of_nodes, edges_around_node, directed_edges, ball_around_node, r):
    cycles = [0] * number_of_nodes
    for i in range(number_of_nodes):
        cycles[i] = directed_edges[i][r] - edges_around_node[i][r] - ball_around_node[i][r] + 1
    return cycles


def plot_cycles(number_of_nodes, cycles, real_cycles, real_edges, real_directed_edges, real_ball, realisation_name, b,
                Realisations):
    fig, ax = plt.subplots()
    if Realisations:
        sortedindex_real_cycles = sorted(range(number_of_nodes), key=lambda index: real_cycles[index])
        vpbounds = 0
        chebbounds = 0

        delta1 = 0.00005
        delta2 = 0.0005
        etam = find_beta(b) / math.sqrt(2 ** b)

        eta3 = etam + delta2

        chebyshevlowervariance = []
        chebyshevuppervariance = []
        vplowervariance = []
        vpuppervariance = []

        for j in sortedindex_real_cycles:
            p = math.sqrt(4 / 9) * 1 / math.sqrt(0.05) * eta3 * math.sqrt(
                len(real_edges[j][1]) ** 2 + len(real_ball[j][1]) ** 2 + len(real_directed_edges[j][1]) ** 2)
            epsilon = p / len(real_edges[j][1]) + delta1
            xi = p / len(real_ball[j][1]) + delta1
            gamma = p / len(real_directed_edges[j][1]) + delta1
            K = gamma * len(real_directed_edges[j][1]) + epsilon * len(real_edges[j][1]) + xi * len(real_ball[j][1])
            vpbounds += K/max(1, real_cycles[j])

            vplowervariance.append(real_cycles[j] - K)
            vpuppervariance.append(real_cycles[j] + K)
            p = 1 / math.sqrt(0.05) * eta3 * math.sqrt(
                len(real_edges[j][1]) ** 2 + len(real_ball[j][1]) ** 2 + len(real_directed_edges[j][1]) ** 2)
            epsilon = p / len(real_edges[j][1]) + delta1
            xi = p / len(real_ball[j][1]) + delta1
            gamma = p / len(real_directed_edges[j][1]) + delta1
            K = gamma * len(real_directed_edges[j][1]) + epsilon * len(real_edges[j][1]) + xi * len(real_ball[j][1])

            chebyshevuppervariance.append(real_cycles[j] - K)
            chebyshevlowervariance.append(real_cycles[j] + K)
            chebbounds += K/max(1, real_cycles[j])

        print('VP:', vpbounds / number_of_nodes)
        print('Cheb:', chebbounds / number_of_nodes)

        # ------------------------------------ plot conductance --------------------------------------

        sorted_cycles = []

        for i in sortedindex_real_cycles:
            sorted_cycles.append(cycles[i])
        ax.plot(range(number_of_nodes), chebyshevlowervariance, '-', color = 'LimeGreen',
                label=str('Chebyshev'))
        ax.plot(range(number_of_nodes), chebyshevuppervariance, '-', color = 'LimeGreen')
        ax.plot(range(number_of_nodes), vplowervariance, '-', color = 'Orange',
                label=str('VP'))
        ax.plot(range(number_of_nodes), vpuppervariance, '-', color = 'Orange')
        data = np.msort(cycles)
        print('DIP for cycles gives us', dip.diptst(data))

    else:
        sorted_cycles = sorted(cycles)

    ax.plot(range(number_of_nodes), sorted_cycles, 'o', label="Estimate")
    if Realisations:
        ax.plot(range(number_of_nodes), sorted(real_cycles), '-', color = 'OrangeRed', label="Realisation")

    name = str("Pictures/" + realisation_name + '_cycles_vp_inequality_and_chebyshev.pdf')
    plt.xlabel('Node $v$')
    plt.ylabel("$C(S_2(v))$")
    #plt.title(str('Number of cycles of length 3 or 4'))
    ax.legend()
    plt.savefig(name)
    plt.show()


    if Realisations:
        fig, ax = plt.subplots()

        differences = [real_cycles[i] - cycles[i] for i in range(number_of_nodes)]
        kwargs = dict(alpha=0.5, bins=50, density=False, stacked=True)

        # Plot
        plt.hist(differences, **kwargs)
        #plt.gca().set(title='Error in estimated number of cycles in $S_1(v)$')
        plt.xlabel('Realisation - Estimate')
        ax.grid()
        plt.xlim(-5, 5)
        name = str(
            "Pictures/" + realisation_name + '_cycles_difference_histogram_distance_1.pdf')
        plt.savefig(name)
        plt.show()


def find_transitivity(triangles, wedges, number_of_nodes):
    transitivity = list([] for i in range(number_of_nodes))
    for distance in [0, 1, 2]:
        for i in range(number_of_nodes):
            transitivity[i].append((3 * triangles[i][distance]) / wedges[i][distance])
    return transitivity


def find_real_transitivity(real_wedges, real_triangles, number_of_nodes):
    real_transitivity = list([] for i in range(number_of_nodes))
    for distance in [0, 1, 2]:
        for i in range(number_of_nodes):
            real_transitivity[i].append(3 * len(real_triangles[i][distance]) / real_wedges[i][distance])
    return real_transitivity


def plot_transitivity(number_of_nodes, transitivity, real_transitivity, real_wedges, real_triangles, distance,
                      realisation_name, Realisations, b):
    fig, ax = plt.subplots()
    delta1 = 0.00005
    delta2 = 0.0005
    etam = find_beta(b) / math.sqrt(2 ** b)
    eta3 = etam + delta2
    print(eta3)
    First = True
    blub = True
    if Realisations:
        for distance in [distance]:
            chebyshevlowervariance = []
            chebyshevuppervariance = []
            vplowervariance = []
            vpuppervariance = []
            sortedindex_real_transitivity = sorted(range(number_of_nodes),
                                                   key=lambda index: real_transitivity[index][distance])
            sorted_transitivity = []
            sortedreal_transitivity = []

            for i in sortedindex_real_transitivity:
                sorted_transitivity.append(transitivity[i][distance])
                sortedreal_transitivity.append(real_transitivity[i][distance])
                if transitivity[i][distance] == 0:
                    chebyshevlowervariance.append(real_transitivity[i][distance])
                    chebyshevuppervariance.append(real_transitivity[i][distance])
                    vplowervariance.append(real_transitivity[i][distance])
                    vpuppervariance.append(real_transitivity[i][distance])
                else:
                    p = math.sqrt(eta3**2 *
                        (real_wedges[i][distance] ** 2 + len(real_triangles[i][distance]) ** 2) / 0.05)
                    epsilon = p / real_wedges[i][distance] + delta1
                    delta = p / len(real_triangles[i][distance]) + delta1
                    p1 = 2 / 3 * math.sqrt(eta3**2 * (real_wedges[i][distance] ** 2 + len(
                                                      real_triangles[i][distance]) ** 2) / 0.05)
                    epsilon1 = p1 / real_wedges[i][distance] + delta1
                    delta2 = p1 / len(real_triangles[i][distance]) + delta1
                    chebyshevlowervariance.append(real_transitivity[i][distance] * (1 - delta) / (1 + epsilon))
                    chebyshevuppervariance.append(real_transitivity[i][distance] * (1 + delta) / (1 - epsilon))
                    vplowervariance.append(real_transitivity[i][distance] * (1 - delta2) / (1 + epsilon1))
                    vpuppervariance.append(real_transitivity[i][distance] * (1 + delta2) / (1 - epsilon1))
            if blub:
                ax.plot(range(number_of_nodes), chebyshevlowervariance, '-', color = 'LimeGreen',
                        label=str('Chebyshev'))
                ax.plot(range(number_of_nodes), vplowervariance, '-', color = 'Orange',
                        label=str('VP'))
                blub = False
            else:
                ax.plot(range(number_of_nodes), chebyshevlowervariance, '-', color = 'LimeGreen',
                        )
                ax.plot(range(number_of_nodes), vplowervariance, '-', color = 'Orange'
                        )
            ax.plot(range(number_of_nodes), chebyshevuppervariance, '-', color = 'LimeGreen')
            ax.plot(range(number_of_nodes), vpuppervariance, '-', color = 'Orange')
            ax.plot(range(number_of_nodes), sorted_transitivity, 'o',
                    label=str("Estimate"))
            if Realisations:
                if First:
                    ax.plot(range(number_of_nodes), sortedreal_transitivity, '-', color = 'OrangeRed', label="Realisation")
                    First = False
                else:
                    ax.plot(range(number_of_nodes), sortedreal_transitivity, '-', color = 'OrangeRed')

    name = str("Pictures/" + realisation_name + '_transitivity_distance_' + str(distance + 1) + '.pdf')

    plt.ylabel(str("Transitivity $t(S_" + str(distance +1) + "(v)$"))
    plt.xlabel('Node $v$')
    #plt.title(str('Transitivity in $S_' + str(distance + 1) + '(v)$'))
    ax.legend()
    plt.savefig(name)
    plt.show()


def plot_transitivity_vs_conductance(real_conductance, real_transitivity, realisation_name):
    fig, ax = plt.subplots()

    ax.plot(real_conductance, real_transitivity, 'o')
    name = str("Pictures/" + realisation_name + '_transitivity_vs_conductance.pdf')

    #plt.title(str('Transitivity vs conductance'))
    plt.ylabel('Transitivity')
    plt.xlabel('Conductance')
    plt.savefig(name)
    plt.show()

def plot_wedges(number_of_nodes, realisation_name, wedges_around_node, real_wedges, Realisations, b, distance=1):
    lowerbound = []
    upperbound = []
    delta2 = 0.0005
    etam = find_beta(b) / math.sqrt(2 ** b)
    eta = etam + delta2

    fig, ax = plt.subplots()
    if Realisations:
        sortedindex_real_wedges = sorted(range(number_of_nodes), key=lambda index: real_wedges[index][distance])
        sorted_wedges = []
        sorted_real_wedges = []
        for i in sortedindex_real_wedges:
            sorted_wedges.append(wedges_around_node[i][distance])
            sorted_real_wedges.append(real_wedges[i][distance])
            p4 = real_wedges[i][distance] * eta / math.sqrt(0.05)
            lowerbound.append(real_wedges[i][distance] - p4)
            upperbound.append(real_wedges[i][distance] + p4)

        ax.plot(range(number_of_nodes), lowerbound, 'g-',
                label=str("Chebyshev"))
        ax.plot(range(number_of_nodes), upperbound, 'g-')
    else:
        sorted_wedges = []
        for i in range(number_of_nodes):
            sorted_wedges.append(wedges_around_node[i][distance])
    sorted_wedges = sorted(sorted_wedges)
    ax.plot(range(number_of_nodes), sorted_wedges, 'o', label="Estimate")
    if Realisations:
        ax.plot(range(number_of_nodes), sorted_real_wedges, '-', color = 'OrangeRed', label="Realisation")

    name = str("Pictures/" + realisation_name + '_wedges_distance_' + str(distance) + '.pdf')

    #plt.title(str('Wedges in $S_' + str(distance) + '(v)$'))
    plt.ylabel('#wedges')
    plt.xlabel('Node')
    ax.legend()
    plt.savefig(name)
    plt.show()


def plot_wedges_triangles(wedges_around_node, triangles_around_node, distance, number_of_nodes, realisation_name):

    fig, ax = plt.subplots()
    for distance in [2, 1, 0]:
        wedges = []
        triangles = []
        for i in range(number_of_nodes):
            wedges.append(wedges_around_node[i][distance])
            triangles.append(triangles_around_node[i][distance])
        ax.plot(wedges, triangles, 'o', label = str('$S_' + str(distance + 1) + '(v)$'))
    name = str("Pictures/" + realisation_name + '_wedges_vs_triangles_distance.pdf')

    #plt.title(str('Triangles versus wedges over different radii $r$'))
    plt.ylabel('$|\hat{\Delta}_r(v)|$')
    plt.xlabel('$|\hat{w}(S_r(v))|$')
    plt.legend()
    plt.savefig(name)
    plt.show()


def generic_plot(x, y, realisation_name, specific_name, xlabel, ylabel):
    fig, ax = plt.subplots()
    ax.plot(range(len(x)), range(len(x)), '-')
    ax.plot(x, y, 'o')
    name = str(
        "Pictures/" + realisation_name + specific_name + '.pdf')
    plt.xlabel(str(xlabel))
    plt.ylabel(str(ylabel))
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    #plt.title(specific_name)
    plt.savefig(name)
    plt.show()


def plot_conductance_of_multiple_distances(number_of_nodes, conductance, realisation_name, degreelist):
    sortedindex = sorted(range(number_of_nodes),
                         key=lambda index: degreelist[index])

    fig, ax = plt.subplots()
    First = True
    for node in sortedindex[:5]:
        lijst = [1]
        for distance in range(3):
            lijst.append(conductance[node][distance])
            data = np.msort([conductance[node][distance] for node in range(number_of_nodes)])
            print('DIP for conductance in S_', distance, '(v) gives us', dip.diptst(data))
        if First:
            ax.plot([i for i in range(4)], lijst, 'o', color = 'DeepSkyblue' ,markersize = 6, label=str('Low degree nodes'))
            ax.plot([i for i in range(4)], lijst, ':', color='lightblue')
            First = False
        else:
            ax.plot([i for i in range(4)], lijst, 'o', color = 'DeepSkyblue', markersize = 6,)
            ax.plot([i for i in range(4)], lijst, ':', color='lightblue')

    First = True
    for node in sortedindex[-5:]:
        lijst = [1]
        for distance in range(3):
            lijst.append(conductance[node][distance])
        if First:
            ax.plot([i for i in range(4)], lijst, 'o', markersize = 6,color = 'OrangeRed', label=str('High degree nodes'))
            ax.plot([i for i in range(4)], lijst, ':', color='pink')
            First = False
        else:
            ax.plot([i for i in range(4)], lijst, 'o', markersize = 6,color = 'OrangeRed')
            ax.plot([i for i in range(4)], lijst, ':', color='pink')

    ax.set(xlabel='Radius $r$', ylabel='Conductance $\hat{\phi}(S_r(v))$')#, title=str('com-Amazon: Conductance in different radii $r$'))
    name = str('Pictures/' + realisation_name + '_conductance_of_all_distances.pdf')
    major_ticks = np.arange(0, 4, 1)
    ax.set_xticks(major_ticks)
    ax.legend()
    plt.savefig(name)
    plt.show()


def estimate_time(b, edges, nodes):
    if b == 10:
        min_initialization = min(0.0005 * edges, 0.0026 * nodes)
        min_iteration_time = min(0.0001 * edges, 0.00006 * nodes)
        max_initialization = max(0.0005 * edges, 0.0026 * nodes)
        max_iteration_time = max(0.0001 * edges, 0.00006 * nodes)
    elif b == 12:
        min_initialization = min(0.0006 * edges, 0.0029 * nodes)
        min_iteration_time = min(0.0002 * edges, 0.0001 * nodes)
        max_initialization = max(0.0006 * edges, 0.0029 * nodes)
        max_iteration_time = max(0.0002 * edges, 0.0001 * nodes)
    elif b == 14:
        min_initialization = min(0.0007 * edges, 0.0036 * nodes)
        min_iteration_time = min(0.0004 * edges, 0.0022 * nodes)
        max_initialization = max(0.0007 * edges, 0.0036 * nodes)
        max_iteration_time = max(0.0004 * edges, 0.0022 * nodes)
    elif b == 16:
        min_initialization = min(0.00013 * edges, 0.0067 * nodes)
        min_iteration_time = min(0.00013 * edges, 0.0066 * nodes)
        max_initialization = max(0.00013 * edges, 0.0067 * nodes)
        max_iteration_time = max(0.00013 * edges, 0.0066 * nodes)
    else:
        return

    print(f"The initialization is going to take between {min_initialization:.1f} and {max_initialization:.1f} seconds")
    print(f"One iteration is going to take between {min_iteration_time:.1f} and {max_iteration_time:.1f} seconds")
    print(f"So this program will be finished in {max_initialization + 3 * max_iteration_time:.1f} seconds, which is "
          f"{(max_initialization + 3 * max_iteration_time) / 60:.1f} minutes")


def plot_vp_chebyshev_triangles(b, realisation_name):
    cheb = []
    vp = []
    delta2 = 0.0005
    etam = find_beta(b) / math.sqrt(2 ** b)
    eta = etam + delta2
    for i in range(8, 24):
        etam = find_beta(i) / math.sqrt(2 ** i)
        eta = etam + delta2
        cheb.append(eta / math.sqrt(0.05) * 100)
        vp.append(2 / 3 * eta / math.sqrt(0.05) * 100)
    fig, ax = plt.subplots()
    ax.plot(range(8, 24), cheb, label='Chebyshev')
    ax.plot(range(8, 24), vp, label='VP')

    name = str("Pictures/" + realisation_name + '_chebyshev_and_VP_bounds_in_triangles.pdf')

    #plt.title(str('Size of Chebyshev and VP bounds in triangles'))
    plt.ylabel('Percentage of number of triangles')
    plt.xlabel('$b$ ($p = 2^b$ registers)')
    major_ticks = np.arange(8, 24, 2)
    ax.set_xticks(major_ticks)
    plt.legend()
    plt.savefig(name)
    plt.show()

def difference_edges_directed_edges(edgeball_around_node, real_edges_around_node, directed_edgeball_around_node, real_directed_edges, number_of_nodes, realisation_name):
    difference_edges = [edgeball_around_node[node][1] - len(real_edges_around_node[node][1]) for node in
                        range(number_of_nodes)]
    difference_directed_edges = [directed_edgeball_around_node[node][1] - len(real_directed_edges[node][1]) for node in
                                 range(number_of_nodes)]

    fig, ax = plt.subplots()
    ax.plot(difference_directed_edges, difference_edges, '+')
    #plt.title('Error in estimated number of edges and directed edges in $S_1(v)$')
    plt.xlabel('Error in number of directed edges')
    plt.ylabel('Error in number of edges')
    ax.grid()
    plt.show()
    name = str('Pictures/' + realisation_name + '_difference_edges_S1(v).pdf')
    plt.savefig(name)

    fig, ax = plt.subplots()
    ax.plot(difference_directed_edges, difference_edges, '+')
    #plt.title('Error in estimated number of edges and directed edges in $S_1(v)$')
    plt.xlabel('Error in number of directed edges')
    plt.ylabel('Error in number of edges')
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    ax.grid()
    plt.show()
    name = str('Pictures/' + realisation_name + '_difference_edges_S1(v)_zoomed_in.pdf')
    plt.savefig(name)


    difference_edges = [edgeball_around_node[node][2] - len(real_edges_around_node[node][2]) for node in
                        range(number_of_nodes)]
    difference_directed_edges = [directed_edgeball_around_node[node][2] - len(real_directed_edges[node][2]) for node in
                                 range(number_of_nodes)]

    print(difference_edges)
    print(difference_directed_edges)

    fig, ax = plt.subplots()
    ax.plot(difference_directed_edges, difference_edges, '+')
    #plt.title('Error in estimated number of edges and directed edges in $S_2(v)$')
    plt.xlabel('Error in number of directed edges')
    plt.ylabel('Error in number of edges')
    ax.grid()
    plt.show()
    name = str('Pictures/' + realisation_name + '_difference_edges_S2(v).pdf')
    plt.savefig(name)


    fig, ax = plt.subplots()
    ax.plot(difference_directed_edges, difference_edges, '+')
    #plt.title('Error in estimated number of edges and directed edges in $S_2(v)$')
    plt.xlabel('Error in number of directed edges')
    plt.ylabel('Error in number of edges')
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    ax.grid()
    plt.show()
    name = str('Pictures/' + realisation_name + '_difference_edges_S2(v)_zoomed_in.pdf')
    plt.savefig(name)

