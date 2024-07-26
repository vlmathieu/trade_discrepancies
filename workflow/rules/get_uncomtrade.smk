rule get_uncomtrade:
    output:
        'resources/raw_data/uncomtrade_{cmd}.parquet.gzip'
    params:
        years       = list(map(str,range(config['years']['start'], config['years']['stop']))),
        cmdCode     = "{cmd}",
        flowCode    = config['flowCode'],
        apikey      = os.environ['comtrade_apikey']
    threads: 1
    conda:
        '../envs/comtradeapicall.yaml'
    script:
        '../scripts/get_uncomtrade_data.py'