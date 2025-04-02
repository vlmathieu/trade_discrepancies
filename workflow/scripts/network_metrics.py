from snakemake.script import snakemake
import polars as pl
import pickle
import numpy as np
import networkx as nx
from scipy.stats import kurtosis
from scipy.stats import skew
import itertools


def composition_stats(sources: list, targets: list) -> dict:
    '''
    Functions that computes network composition descriptive statistics (number 
    of trading countries, number of pure exporters, number of pure importers,
    number of countries that are both exporter and importer) based on lists of
    sources and targets of a directed network. Network relates to trade of one 
    year and one product. Network must be unweighted.

    Parameters
    ----------
    sources : list
        The list of network nodes that are sources = countries that export.
    targets : list
        The list of network nodes that are targets = countries that import.

    Returns
    -------
    compo_dict : dictionary
        Dictionary of network composition descriptive statistics.

    '''

    # Compute list of total number of trading countries
    tot_nb_nodes = list(set(sources + targets))

    # Compute lists of pure sources = exporters / of pure targets = importers
    pure_sources = list(set(sources) - set(targets))
    pure_targets = list(set(targets) - set(sources))

    # Compute list of countries that are both sources and targets = exp & imp
    mixed_src_tgt = list(set(sources).intersection(targets))

    # Compute network composition descriptive statistics based on src/tgt
    compo_dict = {
        'tot_nb_nodes': len(tot_nb_nodes),  # Number of trading countries
        'nb_pure_exp': len(pure_sources),   # Number of pure exporters
        'nb_pure_imp': len(pure_targets),   # Number of pure importers
        'nb_mixed': len(mixed_src_tgt),     # Number of mixed countries
    }

    return compo_dict

def connectivity_stats(degree_list: list) -> dict:
    '''
    Functions that computes network statistics (number of links = trade flows, 
    mean degree = average connectivity, and variance | skeweness | kurtosis of 
    the number of degree = variance | skeweness | kurtosis of connectivity) 
    based on a network degree list. Network relates to trade of one year and one 
    product. Network must be unweighted.

    Parameters
    ----------
    degree_list : list
        The degree list of the network on which network statistics are 
        calculated. 

    Returns
    -------
    stats_dict : dictionary
        Dictionary of network statistics.

    '''

    # Compute network statistics based on degree_list
    stats_dict = {
        'trade_flows': sum(degree_list),        # Number of edges
        'mean_degree': np.mean(degree_list),    # Average degree
        'var_degree': np.var(degree_list),      # Degrees variance
        'skew_degree': skew(degree_list),       # Degrees skewness
        'kurt_degree': kurtosis(degree_list)    # Degrees kurtosis
    }

    return stats_dict

def market_concentration(degree_list: list) -> dict:
    '''
    Functions that computes total circulating value and market concentration 
    indexes (Herfindahl-Hirschmann index, Shannon index) based on a network 
    degree list. Network relates to trade of one year and one product. Network 
    must be weighted.

    Parameters
    ----------
    degree_list : list
        The degree list of the network on which total traded value and market 
        concentration indexes are calculated.

    Returns
    -------
    market_concentration_dict : dictionary
        Dictionary of traded value and market concentration indexes.
        
    '''

    # Compute total circulating value of the network
    value_tot = sum(degree_list)

    market_concentration_dict = {
        'traded_value': value_tot,  # Assign total circulating value
        'hhi': sum([(x/value_tot)**2 
                    for x in degree_list]),  # Compute Herfindahl-Hirschmann index
        'shannon': -sum([((x / value_tot) * np.log(x / value_tot)) 
                         for x in degree_list if x > 0])  # Compute Shannon index
    }

    return market_concentration_dict

def unit_metrics(unit_edge_list_dict: dict, 
                 to_omit: str = None,
                 weight: str = 'primary_value_deflated') -> tuple:
    '''
    Functions that computes network statistics, total circulating value, and 
    market concentration indexes based on an edge list. Network relates to trade 
    of one year and one product and is directed. Network metrics computation may
    omit a specific node using to_omit argument.

    Parameters
    ----------
    unit_edge_list_dict : dictionnary
        A dictionnary that associates (i) a tuple (product, year) of the product 
        code and the year of trade and (ii) the associated edge list describing 
        the network and on which network statistics, total traded value and 
        market concentration indexes are calculated.
    to_omit : string
        The complete name of the country to omit in metrics computation, if
        wanted (NOT THE ISO CODE, except edge_lists built with iso_code).
    weight : string
        The name of the weight to consider. May be 'net_wgt', 'primary_value' or
        'primary_value_deflated'. Default value is 'primary_value_deflated'.

    Returns
    -------
    result, schema : tuple
        Tuple of a list that stores the calculated network statistics, total 
        traded value and market concentration indexes, and a list of column 
        names.
        
    '''

    # Extract keys and edge list
    [[keys, edge_list]] = unit_edge_list_dict.items()

    # Extract product code and year
    product, year = keys

    # Build directed network based on edge_list
    net = nx.from_edgelist(edge_list, create_using=nx.DiGraph)

    # Omit country from network if needed
    if to_omit is not None:
        # If node present in network, remove node
        if net.has_node(to_omit):
            net.remove_node(to_omit)
        # If node not present in network, no omission
        else:
            to_omit = None
    
    # List sources (origin of link = exporters)
    sources = [x for x in net.nodes() if net.out_degree(x) >= 1]

    # List targets (destination of link = importers)
    targets = [x for x in net.nodes() if net.in_degree(x) >= 1]

    # Build unweighted degree lists for sources=exporters | targets=importers
    degree_exp_unweighted = [net.out_degree(x) for x in sources]
    degree_imp_unweighted = [net.in_degree(x) for x in targets]

    # Replace None weights of edges by 0 to work on weighted network
    for u,v,d in net.edges(data=True):
        if d[f'{weight}_exp'] is None:
            d[f'{weight}_exp'] = 0
        if d[f'{weight}_imp'] is None:
            d[f'{weight}_imp'] = 0

    # Build weighted degree lists for sources=exporters | targets=importers
    degree_exp_weighted = [net.out_degree(x, weight = f'{weight}_exp') 
                           for x in sources 
                           # Consider only non-None edges = > 0 weights
                           if net.out_degree(x, weight = f'{weight}_exp') > 0]
    degree_imp_weighted = [net.in_degree(x, weight = f'{weight}_imp') 
                           for x in targets 
                           # Consider only non-None edges = > 0 weights
                           if net.in_degree(x, weight = f'{weight}_imp') > 0]

    # List parameters to compute network statistics, value and market indexes
    params = [('exp', degree_exp_unweighted, degree_exp_weighted), 
              ('imp', degree_imp_unweighted, degree_imp_weighted)]
    
    # Compute network statistics, value and market concentration indexes
    result = (
        [[product, year] + # Product and year of trade
         [f'{s}orter'] + # Type of network
         [to_omit] + # Omission (if specified)
         list(composition_stats(sources, targets).values()) + # Compo. stats.
         list(connectivity_stats(du).values()) + # Connectivity statistics
         list(market_concentration(dw).values()) # Value and mkt concentration
         for s,du,dw in params]
    )
    
    # List column names
    schema = (
        ['product', 'period', 'trader_type', 'omission'] + 
        list(composition_stats(sources, targets)) +
        list(connectivity_stats(degree_exp_unweighted)) + 
        list(market_concentration(degree_exp_weighted))
    )

    return result, schema

def metrics(edge_list_dict: dict, 
            omission_list: list = [None],
            weight: list = 'primary_value_deflated') -> pl.dataframe.frame.DataFrame:
    '''
    Functions that computes network statistics, total circulating value, and 
    market concentration indexes based on an edge list. Network relates to trade 
    of one year and one product and is directed. Network metrics computation may
    omit specific nodes using ommision_list argument.

    Parameters
    ----------
    edge_list_dict : dictionnary
        A dictionnary that associates, for all years and products of trade 
        covered, (i) a tuple (product, year) of the product code and the year of 
        trade and (ii) the associated edge list describing the network and on 
        which network statistics, total traded value and market concentration 
        indexes are calculated.
    omission_list : list of string
        The list if country complete names to omit in metrics computation, if
        wanted (NOT THE ISO CODE, except edge_lists built with iso_code).
    weight : string
        The name of the weight to consider. May be 'net_wgt', 'primary_value' or
        'primary_value_deflated'. Default value is 'primary_value_deflated'.

    Returns
    -------
    result : polars dataframe
        Dataframe of all network metrics per year, per product, and with 
        omission of country if specified.
        
    '''

    # Divide global dictionary into list of unit dictionnaries
    edge_lists = [{k: v} for (k, v) in edge_list_dict.items()]

    # Combine product, year, and to_omit parameters 
    params = itertools.product(edge_lists, omission_list)

    # Apply unit_metrics to every edge list dictionnary accross parameters
    metric_list = [
        unit_metrics(unit_edge_list_dict=unit_edge_list_dict, 
                     to_omit=to_omit, 
                     weight=weight)[0]
        for unit_edge_list_dict, to_omit
        in params
    ]

    # Collect column names
    schema = unit_metrics(edge_lists[0])[1]

    # Build result dataframe
    result = pl.DataFrame(itertools.chain.from_iterable(metric_list), schema)

    # Remove potential duplicates
    result = result.unique()

    return result

# Load dictionary of edge lists
with open(snakemake.input[0], 'rb') as f:
    edge_list_dict = pickle.load(f)

# Compute network metrics based on dictionnary of edge lists
network_metrics = metrics(
    edge_list_dict= edge_list_dict,
    omission_list= snakemake.params['omission_list'],
    weight= snakemake.params['weight']
)

# Save metrics
network_metrics.write_csv(
    snakemake.output[0],
    separator=';'
    )
