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

years, codes = (
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
