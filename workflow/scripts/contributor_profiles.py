from snakemake.script import snakemake
import polars as pl
import pickle
import networkx as nx

def unit_contributor_profiles(unit_edge_list_dict: dict,
                              threshold: float,
                              weight: str = 'primary_value_deflated') -> dict:
    '''
    Function that returns a polar data frame of the profiles of the main 
    contributors to the trade based on an edge list describing a trade network 
    for a given year and traded product. Main contributors are countries that 
    account for at least a certain percentage of the total export or import 
    value. The percentage is defined by the threshold parameter. Contributor 
    profiles include out- and in-degree, total export and import value (current 
    and deflated), and net export and import weight. The network relates to the 
    trade of one year and one product and is directed.

    Parameters
    ----------
    unit_edge_list_dict : dictionnary
        A dictionnary that associates (i) a tuple (product, year) of the product 
        code and the year of trade and (ii) the associated edge list describing 
        the network and on which network statistics, total traded value and 
        market concentration indexes are calculated.
    threshold : float
        The minimum percentage of a country's contribution to the total value of
        exports and imports for a country to be considered a major contributor 
        to trade.
    weight : string
        The weight to which the threshold is applied for a country to be 
        considered as a main contributor to trade. The default value is 
        "primary_value_deflated".

    Returns
    -------
    unit_contributor_profiles : polars data frame
        Dictionary of the profiles of the main contributors to the trade.
        
    '''

    # Extract keys and edge list
    [[keys, edge_list]] = unit_edge_list_dict.items()

    # Extract product code and year
    product, year = keys

    # Build directed network based on edge_list
    net = nx.from_edgelist(edge_list, create_using=nx.DiGraph)

    # Replace None weights by 0
    for _,_,d in net.edges(data=True):
        for key in d:
            if d[key] is None:
                d[key] = 0
    
    # List sources (origin of link = exporters)
    sources = [x for x in net.nodes() if net.out_degree(x) >= 1]

    # List targets (destination of link = importers)
    targets = [x for x in net.nodes() if net.in_degree(x) >= 1]
    
    # Compute total export and import value
    tot_exp = sum(
        [net.out_degree(x, weight = f'{weight}_exp') 
         for x in sources]
    )
    tot_imp = sum(
        [net.in_degree(x, weight = f'{weight}_imp') 
         for x in targets]
    )
    
    # List main contributors to export value
    main_sources = (
        [x for x in net.nodes() 
         if net.out_degree(x, weight = f'{weight}_exp') >= threshold * tot_exp]
    )
    
    # List main contributors to import value
    main_targets = (
        [x for x in net.nodes()
         if net.in_degree(x, weight = f'{weight}_imp') >= threshold * tot_imp]
    )

    # List main contributors to traded value
    countries = list(set(main_sources + main_targets))
    
    # Build contributor profiles as a dictionnary
    contributors_profile = pl.from_dicts([
        {
            "period": year,
            "cmd": product,
            "country": x, 
            "nb_edge_exp": net.out_degree(x),
            "nb_edge_imp": net.in_degree(x),
            "primary_value_exp": net.out_degree(x, weight = 'primary_value_exp'),
            "primary_value_imp": net.in_degree(x, weight = 'primary_value_imp'),
            "primary_value_deflated_exp": 
                net.out_degree(x, weight = 'primary_value_deflated_exp'),
            "primary_value_deflated_imp": 
                net.in_degree(x, weight = 'primary_value_deflated_imp'),
            "net_wgt_exp": net.out_degree(x, weight = 'net_wgt_exp'),
            "net_wgt_imp": net.in_degree(x, weight = 'net_wgt_imp')
        }
        for x in countries
    ])
    
    return contributors_profile

def contributor_profiles(
        edge_list_dict: dict, 
        threshold: float,
        weight: str = 'primary_value_deflated') -> pl.dataframe.frame.DataFrame:
    '''
    Function that returns a polar data frame of the profiles of the main 
    contributors to trade, based on a dictionary of edge lists describing a 
    trade network for each year and product of trade considered. Main 
    contributors are countries that account for at least a certain percentage of
    the total export or import value. The percentage is defined by the threshold
    parameter. Contributor profiles include out- and in-degree, total export and
    import value (current and deflated), and net export and import weight. The 
    network is directed.

    Parameters
    ----------
    edge_list_dict : dictionnary
        A dictionnary that associates, for all years and products of trade 
        covered, (i) a tuple (product, year) of the product code and the year of 
        trade and (ii) the associated edge list describing the network and on 
        which network statistics, total traded value and market concentration 
        indexes are calculated.
    threshold : float
        The minimum percentage of a country's contribution to the total value of
        exports and imports for a country to be considered a major contributor 
        to trade.
    weight : string
        The weight to which the threshold is applied for a country to be 
        considered as a main contributor to trade. The default value is 
        "primary_value_deflated".

    Returns
    -------
    contributor_profiles : polars dataframe
        Dataframe of the profiles of the main contributors to the trade for each
        year and product considered.
        
    '''

    # Divide global dictionary into list of unit dictionnaries
    edge_lists = [{k: v} for (k, v) in edge_list_dict.items()]

    # Apply unit_contributor_profiles to every edge list dictionnary
    contributor_profiles = pl.concat(
        [
            unit_contributor_profiles(
                unit_edge_list_dict = unit_edge_list_dict,
                threshold = threshold,
                weight = weight
            )
            for unit_edge_list_dict in edge_lists
        ]
    )

    return contributor_profiles

# Load dictionary of edge lists
with open(snakemake.input[0], 'rb') as f:
    edge_list_dict = pickle.load(f)

# Compute contributor profiles based on dictionnary of edge lists
contributor_profiles = contributor_profiles(
    edge_list_dict = edge_list_dict,
    threshold = snakemake.params['threshold'],
    weight = snakemake.params['weight']
)

# Save contributor profiles
contributor_profiles.write_csv(
    snakemake.output[0],
    separator=';'
    )