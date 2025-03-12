from snakemake.script import snakemake
import polars as pl

# Load data
merged_data = pl.read_parquet(snakemake.input[0])

# Filter data for network analysis
input_data = (
    merged_data
        .select(snakemake.params['col_keep'])
        .filter(

            # Keep data in specified time range
            pl.col('period') >= str(snakemake.params['year_start']),
            pl.col('period') <= str(snakemake.params['year_stop']-2),

            # Keep imports and exports only
            pl.col('flowCode').is_in(snakemake.params['flow_to_keep']),

            # Keep FAO product division specified
            pl.col('FAO Code').is_in(snakemake.params['fao_divisions']),

            # Remove non-country reporters and partners
            (~pl.col('partnerISO')
             .str.contains('|'.join(snakemake.params['excluded_iso']))),
            (~pl.col('reporterISO')
             .str.contains('|'.join(snakemake.params['excluded_iso']))),

            # Delete "auto-" imports or exports (reporter = partner)
            pl.col('reporterDesc') != pl.col('partnerDesc'),

            # Remove primaryValue_deflated null values (~8% of the dataset)
            ~pl.col('primaryValue_deflated').is_null()
        )
        # Drop potential duplicates
        .unique(subset=snakemake.params['col_keep'])
)

# Drop outliers = values under fifth percentile for weight (kg) and value (USD)
stats_desc = (
    input_data
    .select(['netWgt', 'primaryValue'])
    .describe(percentiles=[0.05])
)

min_weight, min_value = (
    stats_desc
    .filter(pl.col('statistic') == '5%')
    .select(['netWgt', 'primaryValue'])
)

input_data = (
    input_data.filter(

        # Drop trade flow with net weight (kg) under fifth percentile
        pl.col('netWgt') > min_weight.item(),

        # Drop trade flow with value (USD) under fifth percentile
        pl.col('primaryValue') > min_value.item()
    )
)

# Sum weight and values by FAO division product
sum_cols = ['netWgt', 'primaryValue', 'primaryValue_deflated']

groupby_cols = [_ for _ in input_data.columns if _ not in sum_cols]

input_data = (
    input_data
    .group_by(groupby_cols)
    .agg(pl.sum(sum_cols))
 )

# Save input data
input_data.write_parquet(
    snakemake.output[0],
    compression='gzip'
    )
