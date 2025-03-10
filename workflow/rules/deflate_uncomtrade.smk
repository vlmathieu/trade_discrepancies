rule deflate_uncomtrade:
    input:
        'resources/raw_data/wb_series_data.csv',
        'resources/raw_data/uncomtrade_data.parquet.gzip'
    output:
        'results/processed_data/deflate_uncomtrade_data.parquet.gzip'
    threads: 2
    script:
        '../scripts/deflate_uncomtrade.py'