import polars as pl
import polars.selectors as cs
import pickle
import re
import numpy as np
import networkx as nx
from scipy.stats import kurtosis
from scipy.stats import skew
import pandas as pd
import itertools

# Parameters
HS_version = [1996, 2002, 2007, 2012, 2017, 2022]
year_start = 2020
year_stop = 2024
excluded_iso = ['XX', '_X', '\\d']
flow_to_keep = ['M', 'X']
col_keep = ['period', 
           'reporterISO', 
           'reporterDesc', 
           'flowCode', 
           'partnerISO', 
           'partnerDesc',
           'FAO Code',
           'FAO Product',  
           'netWgt', 
           'primaryValue', 
           'primaryValue_deflated']

# Load data
FAO_HS = pl.read_json('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/correspondence_FAO_HS.json')

comtrade_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/public/uncomtrade_data.parquet.gzip')

wb_countries = pl.read_csv('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/wb_countries_data.csv',
                      separator=';')
wb_data = pl.read_csv('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/wb_series_data.csv',
                      separator=';').select(['economy', 'time', 'TM.UVI.MRCH.XD.WD', 'TX.UVI.MRCH.XD.WD'])

deflate_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/global/deflate_uncomtrade_data.parquet.gzip')

merged_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/global/merged_data.parquet.gzip')

input_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/input/input_data.parquet.gzip')

mirror_flows = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/intermediary/mirror_flows.csv',
    separator=';'
)

path = '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/intermediary/edge_lists.pkl'
with open(path, 'rb') as f:
    net_dict = pickle.load(f)

contributor_profiles = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/contributor_profiles.csv',
    separator=';'
)

market_concentration = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/market_concentration.csv',
    separator=';'
)

sorted(market_concentration.select('period').unique().to_series().to_list())
market_concentration.filter(pl.col('cmd') == 12).sort('period').select('hhi_imp')

network_composition = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/network_composition.csv',
    separator=';'
)

network_connectivity = pl.read_csv(
    '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/network_connectivity.csv',
    separator=';'
)

list(net_dict.keys())
dict(net_dict[(12, 2022)])
edge_lists = [{k: v} for (k, v) in net_dict.items()]
unit_edge_list_dict = edge_lists[-1]
type(net_dict)