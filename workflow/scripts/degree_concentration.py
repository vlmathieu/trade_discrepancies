from snakemake.script import snakemake
import polars as pl
import pickle
import numpy as np
import networkx as nx

def unit_degree_concentration(
        unit_edge_list_dict: dict) -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of total trade degree and degree 
    concentration indices for exports and imports based on an edge list 
    describing a trade network for a given year and traded product. The network 
    refers to the trade of one year and one product and is directed and 
    unweighted.

    Parameters
    ----------
    unit_edge_list_dict : dictionnary
        A dictionnary that associates (i) a tuple (cmd, period) of the commodity
        code and the year of trade and (ii) the associated edge list describing 
        the network and on which network total trade degree and degree 
        concentration indices are calculated.

    Returns
    -------
    unit_degree_concentration : polars data frame
        A polars data frame of the total trade degree and degree concentration
        indices for exports and imports for a given year and traded product.
        
    '''

    # Extract keys and edge list from dict
    [[keys, edge_list]] = unit_edge_list_dict.items()

    # Extract commodity code (=cmd) and year (=period)
    cmd, period = keys

    # Build directed network based on edge_list
    net = nx.from_edgelist(edge_list, create_using=nx.DiGraph)
    
    # Build degree lists for exporters=out_degree | importers=in_degree
    degree_weighted_exp = [
        net.out_degree(x) for x in net.nodes() 
        # Consider only non-None edges = degree > 0
        if net.out_degree(x) > 0
    ]
    degree_weighted_imp = [
        net.in_degree(x) for x in net.nodes() 
        # Consider only non-None edges = degree > 0
        if net.in_degree(x) > 0
    ]

    # Compute total trade degree (same for exports and imports)
    trade_degree = sum(degree_weighted_exp)

    # Build dictionnary of tot trade degree and degree concentration indices
    unit_degree_concentration = pl.from_dict(
            {
                "period": period,
                "cmd": cmd,
                # Assign total trade degree for exports and imports
                "trade_degree": trade_degree,
                # Compute Herfindahl-Hirschmann index for exports and imports
                "hhi_exp": sum([(x/trade_degree)**2 
                                for x in degree_weighted_exp]),
                "hhi_imp": sum([(x/trade_degree)**2 
                                for x in degree_weighted_imp]),
                # Compute Shannon index for exports and imports
                'shannon_exp': -sum(
                    [((x / trade_degree) * np.log(x / trade_degree)) 
                     for x in degree_weighted_exp if x > 0]),
                'shannon_imp': -sum(
                    [((x / trade_degree) * np.log(x / trade_degree)) 
                     for x in degree_weighted_imp if x > 0])
            }
    )

    return unit_degree_concentration

def degree_concentration(
        edge_list_dict: dict) -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of total trade degree and degree 
    concentration indices for exports and imports based on a dictionary of edge 
    lists describing a trade network for each year and product of trade 
    considered. The network is directed and unweighted.

    Parameters
    ----------
    edge_list_dict : dictionnary
        A dictionnary that associates, for all years and products of trade 
        covered, (i) a tuple (product, year) of the product code and the year of 
        trade and (ii) the associated edge list describing the network and on 
        which network total trade degree and degree concentration indices are 
        calculated.

    Returns
    -------
    degree_concentration : polars data frame
        A polars data frame of the total trade degree and degree concentration
        indices for exports and imports for each year and product considered.
        
    '''

    # Divide global dictionary into list of unit edge list dictionnaries
    edge_lists = [{k: v} for (k, v) in edge_list_dict.items()]

    # Apply unit_degree_concentration to every edge list dictionnary
    degree_concentration = pl.concat(
        [
            unit_degree_concentration(
                unit_edge_list_dict = unit_edge_list_dict
            )
            for unit_edge_list_dict in edge_lists
        ],
        how = 'vertical_relaxed'
    )

    # Sort by cmd and year
    degree_concentration = degree_concentration.sort(['cmd', 'period'])
    
    return degree_concentration

# Load dictionary of edge lists
with open(snakemake.input[0], 'rb') as f:
    edge_list_dict = pickle.load(f)

# Compute degree concentration stats based on dictionnary of edge lists
degree_concentration = degree_concentration(
    edge_list_dict= edge_list_dict
)

# Save degree concentration stats
degree_concentration.write_csv(
    snakemake.output[0],
    separator=';'
    )