rule plot_network_metrics:
    input:
        'results/processed_data/network_analysis/output/network_composition.csv'
    output:
        'results/processed_data/network_analysis/plot/network_composition.png'
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script: 
        '../scripts/plot_network_metrics.R'