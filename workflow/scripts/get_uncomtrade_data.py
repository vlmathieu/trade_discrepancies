from snakemake.script import snakemake
import polars as pl
import comtradeapicall
from tqdm import tqdm 

def get_uncomtrade_annual(apikey, year, cmd, flow):
    '''
    Function that downloads UN Comtrade data for a given year, a given 
    commodity, and a given trade flow. Need an API key.

    Parameters
    ----------
    apikey : sting
        The API subscription key to download data.
    year : string
        The year of trade.
    cmd : string
        The commodity code.
    flow : string
        The trade flow to download (import, export, re-import, re-export...).

    Returns
    -------
    data : polars dataframe
        The UN Comtrade data for a given year, commodity, and trade flow.

    '''
    
    data = comtradeapicall.getFinalData(
        apikey,
        typeCode        = 'C',          # typeCode(str) : Product type. Goods (C) or Services (S)
        freqCode        = 'A',          # freqCode(str) : The time interval at which observations occur. Annual (A) or Monthly (M)
        clCode          = 'HS',         # clCode(str) : Indicates the product classification used and which version (HS, SITC)
        period          = year,         # period(str) : Combination of year and month (for monthly), year for (annual)
        reporterCode    = None,         # reporterCode(str) : The country or geographic area to which the measured statistical phenomenon relates
        cmdCode         = cmd,          # cmdCode(str) : Product code in conjunction with classification code
        flowCode        = flow,         # flowCode(str) : Trade flow or sub-flow (exports, re-exports, imports, re-imports, etc.)
        partnerCode     = None,         # partnerCode(str) : The primary partner country or geographic area for the respective trade flow
        partner2Code    = None,         # partner2Code(str) : A secondary partner country or geographic area for the respective trade flow
        customsCode     = None,         # customsCode(str) : Customs or statistical procedure
        motCode         = None,         # motCode(str) : The mode of transport used when goods enter or leave the economic territory of a country
        format_output   = 'JSON',       # format_output(str) : The output format. CSV or JSON
        breakdownMode   = 'classic',    # breakdownMode(str) : Option to select the classic (trade by partner/product) or plus (extended breakdown) mode
        includeDesc     = True          # includeDesc(bool) : Option to include the description or not
        )
        
    return data

def chunks(lst, n):
    '''
    Yield successive n-sized chunks from lst.

    Parameters
    ----------
    lst : list
        List to divide in chunks.
    n : integer
        Size of the chunks.
    '''
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_uncomtrade(apikey, years, cmdCode, flowCode):
    '''
    Function that downloads UN Comtrade data for a several years, 
    commodities, and trade flows. Need an API key.

    Parameters
    ----------
    apikey : string
        The API subscription key to download data.
    years : list of strings
        The years of trade.
    cmdCode : list of strings
        The commodity codes.
    flowCode : list of strings
        The trade flow to download (import, export, re-import, re-export...).

    Returns
    -------
    data : polars dataframe
        The UN Comtrade data for a several years, commodity, and trade flows.

    '''

    data_years_batch = (
        [pl.from_pandas(
            get_uncomtrade_annual(
                apikey,
                ','.join(years_batch),
                ','.join(cmd_batch),
                ','.join(flowCode)
            ))
        for years_batch in tqdm(chunks(years, 4), total=sum(1 for _ in chunks(years, 4)))
        for cmd_batch in chunks(cmdCode, 10)]
    )

    data = pl.concat(
        [df for df in data_years_batch if df.shape != (0,0)],
        how='vertical_relaxed'
    )

    return data

UN_Comtrade_data = get_uncomtrade(
    snakemake.params['apikey'],
    snakemake.params['years'],
    snakemake.params['cmdCode'],
    snakemake.params['flowCode']
)

print("\nDataframe head: \n\n", UN_Comtrade_data.head(5), "\n")
print("\nDataframe size (rows, columns): ", UN_Comtrade_data.shape, "\n")

# Check of years, commodities, and different flows considered
check_list = [
    sorted(set(UN_Comtrade_data['period'].unique())) == sorted(set(snakemake.params['years'])),
    sorted(set(UN_Comtrade_data['cmdCode'].unique())) == sorted(set(snakemake.params['cmdCode'])),
    sorted(set(UN_Comtrade_data['flowCode'].unique())) == sorted(set(snakemake.params['flowCode']))
    ]

# Save data
if all(check_list):
    print('Data have been checked.\n')   
    UN_Comtrade_data.write_parquet(
        snakemake.output[0],
        compression='gzip'
        )
else:
    print('Issues found in data download.\n')
