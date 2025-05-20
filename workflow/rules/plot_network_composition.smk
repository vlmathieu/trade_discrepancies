rule plot_network_composition:
    input:
        'results/processed_data/network_analysis/output/network_composition.csv'
    output:
        expand('results/processed_data/network_analysis/plot/{fao_div}/network_composition.png',
               fao_div = config['fao_divisions'])
    params:
        fao_divisions   = config['fao_divisions']
    threads: 1
    conda:
        '../envs/r_plots.yaml'
    script: 
        '../scripts/plot_network_composition.R'