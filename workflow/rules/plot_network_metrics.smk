rule plot_network_metrics:
    input:
        'results/processed_data/network_analysis/output/network_metrics.csv'
    output:
        'results/processed_data/plot/network_analysis/01_net_size.png'
    params:
        r_packages  = config['r_packages']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script: 
        '../scripts/plot_network_metrics.R'