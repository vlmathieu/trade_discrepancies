from snakemake.script import snakemake
import polars as pl
import pickle
import numpy as np
import networkx as nx

def market_concentration_metrics(degree_list: list) -> dict:
    '''
    Functions that compute total circulating value and market concentration 
    indices (Herfindahl-Hirschmann index, Shannon index) based on a 
    network-weighted degree list. The network relates to the trade of one year 
    and one product and is directed and weighted.

    Parameters
    ----------
    degree_list : list
        The degree list of the network on which total traded value and market 
        concentration indices are calculated.

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

def unit_market_concentration(
        unit_edge_list_dict: dict, 
        weight: str = 'primary_value_deflated') -> tuple:
    '''
    Function that returns a polar data frame of total traded value and market 
    concentration indices for exports and imports based on an edge list 
    describing a trade network for a given year and traded product. The weight 
    to be taken into account when calculating the total traded value and market 
    concentration indices is given by the weight parameter. The network relates 
    to the trade of one year and one product and is directed and weighted.

    Parameters
    ----------
    unit_edge_list_dict : dictionnary
        A dictionnary that associates (i) a tuple (product, year) of the product 
        code and the year of trade and (ii) the associated edge list describing 
        the network and on which network statistics, total traded value and 
        market concentration indexes are calculated.
    weight : string
        The weight to which the threshold is applied for a country to be 
        considered as a main contributor to trade. The default value is 
        "primary_value_deflated".

    Returns
    -------
    unit_market_concentration : polars data frame
        A polars data frame of the total traded value and market concentration
        indices for exports and imports.
        
    '''

    # Extract keys and edge list
    [[keys, edge_list]] = unit_edge_list_dict.items()

    # Extract product code and year
    cmd, year = keys

    # Build directed network based on edge_list
    net = nx.from_edgelist(edge_list, create_using=nx.DiGraph)
    
    # Replace None weights by 0
    for _,_,d in net.edges(data=True):
        for key in d:
            if d[key] is None:
                d[key] = 0
    
    # Build weighted degree lists for exporters=out_degree | importers=in_degree
    degree_weighted_exp = [net.out_degree(x, weight = f'{weight}_exp') 
                           for x in net.nodes() 
                           # Consider only non-None edges = > 0 weights
                           if net.out_degree(x, weight = f'{weight}_exp') > 0]
    degree_weighted_imp = [net.in_degree(x, weight = f'{weight}_imp') 
                           for x in net.nodes() 
                           # Consider only non-None edges = > 0 weights
                           if net.in_degree(x, weight = f'{weight}_imp') > 0]

    # Compute total traded value for exports and imports
    traded_value_exp = sum(degree_weighted_exp)
    traded_value_imp = sum(degree_weighted_imp)

    # Build dictionnary of tot traded value and market concentration indices
    unit_market_concentration = pl.from_dicts(
        [
            {
                "period": year,
                "cmd": cmd,
                # Assign total circulating value for exports and imports
                "traded_value_exp": traded_value_exp,
                "traded_value_imp": traded_value_imp,
                # Compute Herfindahl-Hirschmann index for exports and imports
                "hhi_exp": sum([(x/traded_value_exp)**2 
                                for x in degree_weighted_exp]),
                "hhi_imp": sum([(x/traded_value_imp)**2 
                                for x in degree_weighted_imp]),
                # Compute Shannon index for exports and imports
                'shannon_exp': -sum(
                    [((x / traded_value_exp) * np.log(x / traded_value_exp)) 
                     for x in degree_weighted_exp if x > 0]),
                'shannon_imp': -sum(
                    [((x / traded_value_imp) * np.log(x / traded_value_imp)) 
                     for x in degree_weighted_imp if x > 0])
            }
        ]
    )

    return unit_market_concentration

def market_concentration(
        edge_list_dict: dict, 
        weight: list = 'primary_value_deflated') -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of total traded value and market 
    concentration indices for exports and imports based on an edge list 
    describing a trade network for each year and traded product considered. 
    The weight to be taken into account when calculating the total traded value 
    and market concentration indices is given by the weight parameter. The 
    network is directed and weighted.

    Parameters
    ----------
    edge_list_dict : dictionnary
        A dictionnary that associates, for all years and products of trade 
        covered, (i) a tuple (product, year) of the product code and the year of 
        trade and (ii) the associated edge list describing the network and on 
        which network statistics, total traded value and market concentration 
        indexes are calculated.
    weight : string
        The weight to which the threshold is applied for a country to be 
        considered as a main contributor to trade. The default value is 
        "primary_value_deflated".

    Returns
    -------
    market_concentration : polars data frame
        A polars data frame of the total traded value and market concentration
        indices for exports and imports for each year and product considered.
        
    '''

    # Divide global dictionary into list of unit dictionnaries
    edge_lists = [{k: v} for (k, v) in edge_list_dict.items()]

    # Apply unit_metrics to every edge list dictionnary accross parameters
    market_concentration = pl.concat(
        [
            unit_market_concentration(
                unit_edge_list_dict = unit_edge_list_dict, 
                weight = weight
            )
            for unit_edge_list_dict in edge_lists
        ],
        how = 'vertical_relaxed'
    )

    return market_concentration

# Load dictionary of edge lists
with open(snakemake.input[0], 'rb') as f:
    edge_list_dict = pickle.load(f)

# Compute network metrics based on dictionnary of edge lists
market_concentration = market_concentration(
    edge_list_dict= edge_list_dict,
    weight= snakemake.params['weight']
)

# Save metrics
market_concentration.write_csv(
    snakemake.output[0],
    separator=';'
    )