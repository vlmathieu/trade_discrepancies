rule merge_data:
    input:
        'correspondence_FAO_HS.json',
        'resources/raw_data/wb_series_data.csv',
        'resources/raw_data/wb_countries_data.csv',
        'results/processed_data/deflate_uncomtrade_data.parquet.gzip'
    output:
        'results/processed_data/merged_data.parquet.gzip'
    params:
        wb_series_drop      = ['TM.UVI.MRCH.XD.WD', 'TX.UVI.MRCH.XD.WD'],
        wb_countries_keep   = ['id', 'longitude', 'latitude','capitalCity'],
    threads: 2
    script:
        '../scripts/merge_data.py'