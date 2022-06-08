import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.cm as cm
import matplotlib.colors as colors
import colorsys
from itertools import cycle
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##########################################################################################


def create_adjency_matrix(df_tsv, column):
    """Create adjency matrix from a TSV file .
    Args:
        df_tsv (pandas.DataFrame): output format of the blasttab outfmt 6
        column (string): name of the column to be the edge value
    Returns:
        pandas.DataFrame: adjavency matrix
    """

    all_proteins = list(set(df_tsv["protein1"].tolist() + df_tsv["protein2"].tolist()))

    df_tsv = df_tsv[~(df_tsv["protein1"] == df_tsv["protein2"])]

    df = pd.pivot_table(
        data=df_tsv, index="protein1", columns="protein2", values=column, fill_value=0
    )

    df = df.reindex(columns=all_proteins, fill_value=0)
    df = df.reindex(all_proteins, fill_value=0)

    return df


##########################################################################################


def get_color_cmap(name, n_colors=6):

    """
    Return discrete colors from a matplotlib palette.

    :param name: Name of the palette. This should be a named matplotlib colormap.
    :type: str
    :param n_colors: Number of discrete colors in the palette.
    :type: int
    :return: List-like object of colors as hexadecimal tuples
    :type: list
    """

    brewer_qual_pals = {
        "Accent": 8,
        "Dark2": 8,
        "Paired": 12,
        "Pastel1": 9,
        "Pastel2": 8,
        "Set1": 9,
        "Set2": 8,
        "Set3": 12,
        "tab20": 20,
        "tab20b": 20,
    }

    if name == "tab20" and n_colors > 19:
        second = "tab20b"
        ncolor2 = n_colors - 19
        n_colors = 19
    else:
        second = False

    cmap = getattr(cm, name)

    if name in brewer_qual_pals:
        bins = np.linspace(0, 1, brewer_qual_pals[name])
        if "tab20" == name:
            len_bins = len(bins)
            bins = [bins[i] for i in range(len_bins) if i != 14][:n_colors]
        else:
            bins = bins[:n_colors]
    else:
        bins = np.linspace(0, 1, n_colors + 2)[1:-1]

    palette = list(map(tuple, cmap(bins)[:, :3]))

    if second:
        cmap = getattr(cm, second)
        bins = np.linspace(0, 1, brewer_qual_pals[second])[:ncolor2]
        palette += list(map(tuple, cmap(bins)[:, :3]))

        pal_cycle = cycle(palette)
        palette = [next(pal_cycle) for _ in range(n_colors + ncolor2)]
    else:
        pal_cycle = cycle(palette)
        palette = [next(pal_cycle) for _ in range(n_colors)]

    return [colors.rgb2hex(rgb) for rgb in palette]


##########################################################################


def contrasting_text_color(hex_str):
    """
    Input a string without hash sign of RGB hex digits to compute
    complementary contrasting color such as for fonts
    """

    (r, g, b) = (hex_str[1:3], hex_str[3:5], hex_str[5:])

    luminosity = (
        1 - (int(r, 16) * 0.299 + int(g, 16) * 0.587 + int(b, 16) * 0.114) / 255
    )

    return "#131516" if luminosity < 0.5 else "white"


##########################################################################


def adjust_lightness(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """

    try:
        c = colors.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*colors.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


##########################################################################################


def create_graph(blast_tbl, thresholdName):
    """
    Create the graph from the blast table

    :params adjacency_matrix: adjacency matrix calculate with get_matrix_interaction_system()
    :type: pandas.DataFrame
    :params output: Name of the graphml file
    :type: str
    :return: Nothing

    """

    adj_matrix = create_adjency_matrix(blast_tbl, thresholdName)

    # Create the graph from adjacency matrix
    graph = nx.from_numpy_matrix(adj_matrix.values)
    graph = nx.relabel_nodes(graph, dict(enumerate(adj_matrix.columns)))

    return graph


##########################################################################################


def star_visu_graph(args):
    return visu_graph(*args)


##########################################################################################


def visu_graph(graph, output, threshold, annot_dict, visu_dict, thresholdName, layout, singleton):

    """
    Visualisation of the graph using the weight for the color. If weight >0.5
    put the edge red. And write the graph in graphMl format.

    :params adjacency_matrix: adjacency matrix calculate with get_matrix_interaction_system()
    :type: pandas.DataFrame
    :params output: Name of the graphml file
    :type: str
    :return: the number of red edge in the network

    """

    # Remove singleton
    if not singleton:
        outdeg = graph.degree()
        to_remove = [n[0] for n in outdeg if outdeg[n[0]] == 0]
        graph.remove_nodes_from(to_remove)

    # Get name of all nodes/genes
    all_gene = list(graph)

    # Get cotegory name
    all_category = np.unique(list(visu_dict.values()))
    num_category = all_category.shape[0]
    all_category = sorted(all_category)

    palette = get_color_cmap("tab20", n_colors=num_category)
    tmp_color = {all_category[i]: palette[i] for i in range(num_category)}

    dict_color = {}
    for hit_id in all_gene:
        if hit_id in visu_dict:
            dict_color[hit_id] = tmp_color[visu_dict[hit_id]]
        else:
            dict_color[hit_id] = "grey"

    # Parsing the edges and changing the color edge to red if percentage id >= 50% (or any threshold)
    # Having a list of the edges with to have the width proportionnal to percentage id
    edge2remove = []

    # Index that will allows to count the number of red edges
    num_edges_red = 0

    for n1, n2, edge_dict in graph.edges.data():
        if thresholdName == "evalue":
            test = edge_dict["weight"] <= threshold
        else:
            test = edge_dict["weight"] >= threshold

        if test:
            if (
                n1 in annot_dict
                and n2 in annot_dict
                and annot_dict[n1] == annot_dict[n2]
            ):
                graph.edges[(n1, n2)]["color"] = "#a9a9a9"  # Grey
            else:
                graph.edges[(n1, n2)]["color"] = "#C41E3A"  # Red
                num_edges_red += 1
        else:
            edge2remove.append((n1, n2))

    for e in edge2remove:
        graph.remove_edge(*e)

    graph = graph.to_undirected()

    # Remove singleton
    if not singleton:
        outdeg = graph.degree()
        to_remove = [n[0] for n in outdeg if outdeg[n[0]] == 0]
        graph.remove_nodes_from(to_remove)

    # print("Calculating layout...")

    # Color edges
    edges, edge_colors = zip(*nx.get_edge_attributes(graph, "color").items())

    # Choose between : dot, neato, fdp, sfdp, twopi, circo
    pos = graphviz_layout(graph, prog=layout)

    # Put the color of the node
    nx.set_node_attributes(graph, dict_color, "color")

    # Color nodes
    zip(*nx.get_node_attributes(graph, "color").items())
    nodes, nodes_colors = zip(*nx.get_node_attributes(graph, "color").items())

    # Write the graph
    # nx.write_graphml(graph, output)

    # Get all the color for the border if colored border = True
    nodes_edges_color = []
    for node_color in nodes_colors:
        nodes_edges_color.append(
            "#2F3D44" if node_color == "#FFFFFF" else adjust_lightness(node_color)
        )
        # if you want to choose between colored or not border
        # else :
        #     nodes_edges_color.append('#131516')

    plt.figure(figsize=(12, 10))

    # If you have too much node maybe reduce the node_size to 50
    nx.draw_networkx_nodes(
        graph, pos, node_color=nodes_colors, node_size=50, edgecolors=nodes_edges_color
    )

    nx.draw_networkx_edges(graph, pos, edgelist=edges, edge_color=edge_colors)

    custom_lines = []
    custom_text = []

    for gene, color in tmp_color.items():
        custom_text.append(gene)
        custom_lines.append(
            Line2D(
                range(1),
                range(1),
                color="white",
                marker="o",
                markerfacecolor=color,
                linestyle="none",
            )
        )

    plt.legend(
        custom_lines,
        custom_text,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        prop={"size": "xx-small"},
    )

    # Label drawing as well
    # nx.draw_networkx_labels(graph,pos,font_size=8)

    plt.axis("off")
    plt.tight_layout()

    if thresholdName == "evalue":
        plt.title(f"Blast network from hit of {thresholdName} <= {threshold}")
    else:
        plt.title(f"Blast network from hit of {thresholdName} => {threshold}")

    if output != None:
        plt.savefig(output, dpi=300, bbox_inches="tight")

    plt.close("all")

    return threshold, num_edges_red


##########################################################################
