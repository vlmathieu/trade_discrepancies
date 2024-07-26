from snakemake.script import snakemake
import polars as pl

data_merged = pl.read_parquet(snakemake.input['data'])

data_merged.write_parquet(
    snakemake.output[0],
    compression='gzip'
)
