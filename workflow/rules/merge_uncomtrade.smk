rule merge_uncomtrade:
    input:
        data= lambda wildcards: expand('resources/raw_data/uncomtrade_{cmd}.parquet.gzip', cmd=config['cmdCode'])
    output:
        'resources/raw_data/uncomtrade_merged.parquet.gzip'
    threads: 1
    script:
        '../scripts/merge_uncomtrade.py'
