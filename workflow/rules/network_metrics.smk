rule network_metrics:
    input:
        'results/processed_data/network_analysis/intermediary/edge_lists.pkl'
    output:
        'results/processed_data/network_analysis/output/network_metrics.csv'
    params:
        omission_list   = config['omission_list'],
        weight          = config['weight']
    threads: 4
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/network_metrics.py'