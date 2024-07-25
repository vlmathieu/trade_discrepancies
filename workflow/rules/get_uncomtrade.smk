rule get_uncomtrade:
    output:
        'resources/raw_data/UN_Comtrade_data.parquet.gzip'
    params:
        years       = list(map(str,range(config['years']['start'], config['years']['stop']))),
        cmdCode     = config['cmdCode'],
        flowCode    = config['flowCode'],
        apikey      = os.environ['comtrade_apikey']
    conda:
        '../envs/comtradeapicall.yaml'
    script:
        '../scripts/get_uncomtrade_data.py'