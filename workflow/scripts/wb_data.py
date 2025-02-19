from snakemake.script import snakemake
import polars as pl
import wbgapi as wb

wb_series_data = (
    pl.DataFrame(
        wb.data.DataFrame(snakemake.params['wb_series'], 
                          time=range(
                              snakemake.params['year_start'],
                              snakemake.params['year_stop']), 
                          columns='series')
                          .reset_index()
                )
)

wb_countries_data = wb.economy.DataFrame(skipAggs=True).reset_index()

# Saving correspondence dataframe
wb_series_data.write_parquet(
    snakemake.output[0],
    compression='gzip'
)

wb_countries_data.write_csv(
    snakemake.output[1],
    separator=';'
)
