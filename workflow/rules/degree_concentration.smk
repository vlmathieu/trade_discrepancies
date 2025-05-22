rule degree_concentration:
    input:
        'results/processed_data/network_analysis/intermediary/edge_lists.pkl'
    output:
        'results/processed_data/network_analysis/output/degree_concentration.csv'
    threads: 2
    conda:
        '../envs/network_metrics.yaml'
    script: 
        '../scripts/degree_concentration.py'