import polars as pl
import polars.selectors as cs

HS_version = [1996, 2002, 2007, 2012, 2017, 2022]
year_start = 2020
year_stop = 2024

json = '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/correspondence_FAO_HS.json'
FAO_HS = pl.read_json(json)

# Work on years_chunks
HS_version_keep = HS_version[next(x[0] for x in enumerate(HS_version) if x[1] > year_start)-1:]
# break_years = [_ for _ in HS_version if _ >= year_start]
break_years = sorted(set(
    [_ for _ in HS_version + [year_start, year_stop] if _ >= year_start]
))
years_chunks = []
for i in range(0, len(break_years)-1):
    year_batch = [str(_) for _ in range(break_years[i], break_years[i+1])]
    years_chunks.append(year_batch)

# Work on codes per HS_version
codes = [HS_codes.to_list()
        for HS in HS_version_keep
        for HS_codes in FAO_HS.select(cs.contains(str(HS))).unique()]

# Ziping years_chunks and codes
res = zip(years_chunks, codes)
test = list(res)

test = [str(cmd) for _,cmds in zip(years_chunks, codes) for cmd in cmds]

zip_ = process_param(year_start, year_stop, HS_version, json)

test = [[year, codes] 
        for year, codes in process_param(year_start, year_stop, HS_version, json)]

years = [y 
         for year,_ in process_param(year_start, year_stop, HS_version, json) 
         for y in year]

test = [list(_) for _ in process_param(year_start, year_stop, HS_version, json)]

years, _ = (
    map(list, 
        zip(*[list(_) 
              for _ in process_param(year_start, year_stop, HS_version, json)]))
)
flat_years = sorted(set([str(_) for years_ in years for _ in years_]))
flat_codes = sorted(set([str(_) for codes_ in codes for _ in codes_]))

for years, cmdCode in zip(years_chunks, codes):
    years_ = years
    cmdCode_ = cmdCode

list_sample = [
    {
        'name': 'A',
        'fame': 0,
        'data': {
            'date': ['2021-01-01', '2021-02-01'],
            'credit_score': [800, 890],
            'spend': [1500, 25000],
            'average_spend': 5000
        }
    },
    {
        'name': 'B',
        'fame': 1,
        'data': {
            'date': ['2022-01-01', '2022-02-01', '2022-03-01'],
            'credit_score': [2800, 390, 8900],
            'spend': [15000, 5000, 500],
            'average_spend': 3000
        }
    }
]

(pl.DataFrame(list_sample) 
   .unnest('data')
   .explode('date', 'credit_score', 'spend') 
)

lst1 = ['440110', '440111', '440112', '440121', '440122', '440130', '440131', '440132', '440139', '440140', '440141', '440149', '440200', '440290', '4403', '440310', '440311', '440312', '440320', '440321', '440322', '440323', '440324', '440325', '440326', '440341', '440342', '440349', '440391', '440392', '440393', '440394', '440395', '440396', '440397', '440398', '440399', '4404', '440410', '440420', '4405', '440500', '4406', '440610', '440611', '440612', '440690', '440691', '440692', '4407', '440710', '440711', '440712', '440713', '440714', '440719', '440721', '440722', '440723', '440724', '440725', '440726', '440727', '440728', '440729', '440791', '440792', '440793', '440794', '440795', '440796', '440797', '440799', '4408', '440810', '440831', '440839', '440890', '440910', '440920', '440922', '440929', '4410', '441011', '441012', '441019', '441090', '4411', '441112', '441113', '441114', '441192', '441193', '441194', '4412', '441231', '441232', '441233', '441234', '441239', '441241', '441242', '441249', '441251', '441252', '441259', '441291', '441292', '441294', '441299', '4413', '441300', '4414', '441400', '441410', '441490', '4415', '441510', '441520', '4416', '441600', '4417', '441700', '441810', '441811', '441819', '441820', '441821', '441829', '441830', '441840', '441850', '441860', '441871', '441872', '441874', '441875', '441879', '441881', '441882', '441883', '441889', '441890', '441899', '4419', '441900', '441920', '441990', '4420', '442010', '442011', '442019', '442090', '442110', '442120', '442190', '442199', '4501', '450110', '450190', '4502', '450200', '4503', '450310', '450390', '4504', '450410', '450490', '4701', '470100', '4702', '4703', '470311', '470319', '470321', '470329', '4704', '470411', '470419', '470421', '470429', '4705', '470500', '4706', '470610', '470620', '470630', '470691', '470692', '470693', '4707', '470710', '470720', '470730', '470790', '4801', '480100', '4802', '480210', '480220', '480240', '480254', '480255', '480256', '480257', '480258', '480261', '480262', '480269', '4803', '480300', '4804', '480411', '480419', '480421', '480429', '480431', '480439', '480441', '480442', '480449', '480451', '480452', '480459', '4805', '480511', '480512', '480519', '480524', '480525', '480530', '480540', '480550', '480591', '480592', '480593', '4806', '480610', '480620', '480630', '4807', '480700', '480710', '480790', '4808', '480810', '480820', '480840', '480890', '4809', '4810', '481013', '481014', '481019', '481022', '481029', '481031', '481032', '481039', '481092', '481099', '481110', '481131', '481139', '481141', '481149', '481151', '481159', '481160', '481190', '4812', '481200', '4813', '481310', '481320', '4814', '481420', '481490', '4816', '481620', '481690', '4817', '481710', '481720', '481730', '4818', '481810', '481820', '481830', '481840', '481850', '481890', '4819', '481910', '481920', '481930', '481940', '481950', '481960', '4820', '482010', '482020', '482030', '482040', '482050', '482090', '4821', '482110', '482190', '4822', '482210', '482290', '482320', '482340', '482360', '482369', '482370', '482390', '6808', '680800', '940130', '940131', '940140', '940141', '940161', '940169', '940190', '940191', '940330', '940340', '940350', '940360', '940390', '940391', '9406', '940600', '940610', '9619', '961900'] 
lst2 = ['440110', '440111', '440112', '440121', '440122', '440130', '440131', '440132', '440139', '440140', '440141', '440149', '440200', '440290', '4403', '440310', '440311', '440312', '440320', '440321', '440322', '440323', '440324', '440325', '440326', '440341', '440342', '440349', '440391', '440392', '440393', '440394', '440395', '440396', '440397', '440398', '440399', '4404', '440410', '440420', '4405', '440500', '4406', '440610', '440611', '440612', '440690', '440691', '440692', '4407', '440710', '440711', '440712', '440713', '440714', '440719', '440721', '440722', '440723', '440724', '440725', '440726', '440727', '440728', '440729', '440791', '440792', '440793', '440794', '440795', '440796', '440797', '440799', '4408', '440810', '440831', '440839', '440890', '440910', '440920', '440922', '440929', '4410', '441011', '441012', '441019', '441090', '4411', '441112', '441113', '441114', '441192', '441193', '441194', '4412', '441231', '441232', '441233', '441234', '441239', '441241', '441242', '441249', '441251', '441252', '441259', '441291', '441292', '441294', '441299', '4413', '441300', '4414', '441400', '441410', '441490', '4415', '441510', '441520', '4416', '441600', '4417', '441700', '441810', '441811', '441819', '441820', '441821', '441829', '441830', '441840', '441850', '441860', '441871', '441872', '441874', '441875', '441879', '441881', '441882', '441883', '441889', '441890', '441899', '4419', '441900', '441920', '441990', '4420', '442010', '442011', '442019', '442090', '442110', '442120', '442190', '442199', '4501', '450110', '450190', '4502', '450200', '4503', '450310', '450390', '4504', '450410', '450490', '4701', '470100', '4702', '4703', '470311', '470319', '470321', '470329', '4704', '470411', '470419', '470421', '470429', '4705', '470500', '4706', '470610', '470620', '470630', '470691', '470692', '470693', '4707', '470710', '470720', '470730', '470790', '4801', '480100', '4802', '480210', '480220', '480240', '480254', '480255', '480256', '480257', '480258', '480261', '480262', '480269', '4803', '480300', '4804', '480411', '480419', '480421', '480429', '480431', '480439', '480441', '480442', '480449', '480451', '480452', '480459', '4805', '480511', '480512', '480519', '480524', '480525', '480530', '480540', '480550', '480591', '480592', '480593', '4806', '480610', '480620', '480630', '4807', '480700', '480710', '480790', '4808', '480810', '480820', '480840', '480890', '4809', '4810', '481013', '481014', '481019', '481022', '481029', '481031', '481032', '481039', '481092', '481099', '481110', '481131', '481139', '481141', '481149', '481151', '481159', '481160', '481190', '4812', '481200', '4813', '481310', '481320', '4814', '481420', '481490', '4816', '481620', '481690', '4817', '481710', '481720', '481730', '4818', '481810', '481820', '481830', '481840', '481850', '481890', '4819', '481910', '481920', '481930', '481940', '481950', '481960', '4820', '482010', '482020', '482030', '482040', '482050', '482090', '4821', '482110', '482190', '4822', '482210', '482290', '482320', '482340', '482360', '482369', '482370', '482390', '490640', '6808', '680800', '940130', '940131', '940140', '940141', '940161', '940169', '940190', '940191', '940330', '940340', '940350', '940360', '940390', '940391', '9406', '940600', '940610', '948130', '9619', '961900']
list(set(lst2) - set(lst1))

import polars as pl
import polars.selectors as cs
wb_data = (
    pl.read_csv('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/public/wb_series_data.csv',
                separator=';')
                .drop(['TM.UVI.MRCH.XD.WD', 'TX.UVI.MRCH.XD.WD'])
                .with_columns(pl.col('time').str.replace(r'YR', ''))
)
wb_data.with_columns(pl.col('time').str.replace(r'YR', '')).drop(['TM.UVI.MRCH.XD.WD', 'TX.UVI.MRCH.XD.WD'])
wb_data.describe()
wb_data.describe().row(1, named=True)

wb_data.columns
sorted([_ for _ in wb_data.columns if _ not in {'economy', 'time'}])

comtrade_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/uncomtrade_data.parquet.gzip')

wb_countries = pl.read_csv('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/wb_countries_data.csv',
                      separator=';')
wb_countries.columns
wb_countries.select('aggregate')

wb_data = pl.read_csv('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/wb_series_data.csv',
                      separator=';').select(['economy', 'time', 'TM.UVI.MRCH.XD.WD', 'TX.UVI.MRCH.XD.WD'])

json = '/Users/valentinmathieu/Desktop/wd/trade_discrepancies/resources/raw_data/correspondence_FAO_HS.json'
FAO_HS = pl.read_json(json)
FAO_HS.columns
FAO_HS.select('FAO Code')

# Explore null values in deflate data, primaryValue_deflated
deflate_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/global/deflate_uncomtrade_data.parquet.gzip')
deflate_data.select('primaryValue_deflated').describe()
2.466928e6/2.1949847e7*100

# Check classification code -> delay for adoption of new classfication
deflate_data.columns
deflate_data.select('classificationCode')
(deflate_data
 .filter((pl.col('period').cast(pl.Int32) < 2025)
         & (pl.col('period').cast(pl.Int32) >= 2022))
 .select('classificationCode')
 .unique()
 )
deflate_data.select('cmdCode')

# Merge comtrade data with FAO_HS
FAO_HS_str = FAO_HS.select(cs.starts_with("").cast(pl.String)) # string for join
classifCode = sorted(deflate_data.select('classificationCode').unique().to_series().to_list())
HS_vers = [code for code in FAO_HS_str.columns if '1996' in code]+[code for code in FAO_HS_str.columns if 'HS' in code]
[(classCode, HS_vers) for classCode, HS_vers in zip(classifCode,HS_vers)]
FAO_cat = [_ for _ in FAO_HS_str.columns if 'HS' not in _]

test = (
    [
        (deflate_data
         .filter(pl.col('classificationCode') == classCode)
         .join(FAO_HS_str.unique(subset=FAO_cat+[HS_vers]),
               left_on='cmdCode',
               right_on=HS_vers,
               how='left')
         .drop(cs.starts_with('HS'))
        )
        for classCode, HS_vers in zip(classifCode,HS_vers)
    ]
)

res = pl.concat(
        [df for df in test if df.shape != (0,0)],
        how='vertical_relaxed'
    )
res.shape
deflate_data.shape

res.select(['FAO Code Agg', 'FAO 1982', 'FAO Product', 'FAO Code']).describe()
res.select(['cmdCode', 'FAO Code Agg', 'FAO 1982', 'FAO Product', 'FAO Code'])

res_duplicates = res.filter(pl.struct('period', 'reporterDesc', 'flowDesc', 'partnerDesc', 'cmdCode').is_duplicated())
res_duplicates.select(['period', 'reporterDesc', 'flowDesc', 'partnerDesc', 'cmdCode', 'FAO Code Agg', 'FAO Code'])
sorted(res_duplicates.select('FAO Code').to_series().unique().to_list())
sorted(res_duplicates.select('FAO Product Agg').to_series().unique().to_list())
sorted(res_duplicates.select('FAO Product').to_series().unique().to_list())
sorted(res_duplicates.select(pl.col('FAO Code')).to_series().unique().to_list())
zoom = res_duplicates.filter(pl.col('FAO Code').is_in(['011', '012']))
zoom.select(['period', 'reporterDesc', 'flowDesc', 'partnerDesc', 'cmdCode', 'FAO Code Agg', 'FAO Code'])

# Check merged data
deflate_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/global/deflate_uncomtrade_data.parquet.gzip')
merged_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/global/merged_data.parquet.gzip')
deflate_data.shape
merged_data.columns
merged_data.select(pl.col('FAO Code')).unique()
merged_data.filter(pl.col('FAO Code').is_in(['011', '012']))
input_data = pl.read_parquet('/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/input/input_data.parquet.gzip')

input_data.columns

sum_cols = ['netWgt', 'primaryValue', 'primaryValue_deflated']
with pl.Config(tbl_cols=-1):
    print(input_data.select(sum_cols).describe())

input_data.select('primaryValue_deflated').describe()
with pl.Config(tbl_rows=-1):
    print((
        input_data
     .filter(pl.col('primaryValue_deflated').is_null())
     .select(['period', 'reporterISO', 'reporterDesc', 'flowCode', 'partnerISO', 'partnerDesc', 'netWgt'])
    ))

with pl.Config(tbl_rows=-1):
    print(
        input_data
        .unique(subset=['reporterISO', 'reporterDesc'])
        .select(['reporterISO', 'reporterDesc'])
        .filter(pl.col('reporterISO').str.contains('|'.join(excluded_iso)))
    )