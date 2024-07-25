rule test:
    params:
        years   = config['years']
    log:
        'logs/test.log'
    script:
        'sandbox/toy_script.py'