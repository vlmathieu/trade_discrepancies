rule get_uncomtrade:
    output:
        'resources/raw_data/uncomtrade_data.parquet.gzip'
    params:
        years       = list(map(str,range(config['years']['start'], config['years']['stop']))),
        cmdCode     = list(str(cmd) for cmd in config['cmdCode']),
        flowCode    = list(str(flow) for flow in config['flowCode']),
        apikey      = os.environ['comtrade_apikey']
    threads: 1
    conda:
        '../envs/comtradeapicall.yaml'
    script:
        '../scripts/get_uncomtrade_data.py'