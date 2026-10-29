
from datetime import datetime
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import polars as pl
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import argparse

def read_e3SIM(samples: 'str', infections: 'str', recoveries: 'str'):
    samples_raw = pl.read_csv(samples, 
        has_header=False,
        separator=' ',
        new_columns=['tick', 'host', 'seed'])
    infections_raw = pl.read_csv(infections, 
        separator=" ",
        new_columns=("tick", "infector", "infectee"))
    recoveries_raw = pl.read_csv(recoveries,
        has_header=False,
        separator=" ",
        new_columns=["tick", "host"])
    return [samples_raw, infections_raw, recoveries_raw]

def make_genome_metadata(sample_csv: pl.DataFrame, base_date: pl.datetime):
    return sample_csv.with_columns(
        pl.duration(days=pl.col('tick')).alias("days")
    ).with_columns(
    (pl.col('days') + base_date).alias("date")
    ).with_columns(
        pl.concat_str(
            [pl.col('tick').cast(pl.String), pl.col('host').cast(pl.String)],
            separator="."
        ).alias("sequence_id")
    ).with_columns(
        pl.col("date").dt.to_string(format="%Y-%m-%d")
    ).drop("tick", "days")

def find_recovery_tick(
    infection_tick: int,
    host: int,
    recoveries: pl.DataFrame,
) -> int | None:
    recovery_ticks = recoveries.filter(
        (pl.col("host") == host) & (pl.col("tick") > infection_tick)
    ).select(pl.col("tick").min())
    return recovery_ticks.item()

def make_epi_case_data(
    infections: pl.DataFrame, 
    genome_metadata: pl.DataFrame, 
    recoveries: pl.DataFrame
    ):
    infection_rows = infections.with_row_index("_infection_row")
    recovery_ticks = (
        infection_rows.join(
            recoveries,
            left_on="infectee",
            right_on="host",
            suffix="_recovery",
        )
        .filter(pl.col("tick_recovery") > pl.col("tick"))
        .group_by("_infection_row")
        .agg(pl.col("tick_recovery").min().alias("recovery_tick"))
    )
    epi_data = (
        infection_rows.join(recovery_ticks, on="_infection_row", how="left")
        .sort("_infection_row")
        .drop("_infection_row")
        .with_columns(
            pl.duration(days=pl.col("tick")).alias("days"),
            pl.concat_str(
                [pl.col("tick").cast(pl.String), pl.col("infectee").cast(pl.String)],
                separator=".",
            ).alias("id"),
        )
        .with_columns((pl.col("days") + base_date).alias("date"))
        .with_columns(
            pl.duration(days=pl.col("recovery_tick")).alias("recovery_days")
        )
        .with_columns(
            (pl.col("recovery_days") + base_date).alias("recovery_date")
        )
    )
    sampled_hosts = genome_metadata.select(
        pl.col("sequence_id"),
        pl.col("sequence_id")
        .str.split(".")
        .list.get(1)
        .cast(pl.Int64)
        .alias("sampled_host"),
        pl.col("date").str.to_date("%Y-%m-%d").alias("sample_date"),
    ).unique()
    sampled_infections = (
        epi_data.with_row_index("_infection_row")
        .join(
            sampled_hosts,
            left_on="infectee",
            right_on="sampled_host",
            how="left",
        )
        .filter(
            (pl.col("sample_date") >= pl.col("date").dt.date())
            & (
                pl.col("recovery_date").is_null()
                | (pl.col("sample_date") <= pl.col("recovery_date").dt.date())
            )
        )
        .sort(["_infection_row", "sample_date", "sequence_id"])
        .group_by("_infection_row")
        .agg(
            pl.col("sequence_id").first(),
            pl.col("sample_date").first().alias("sampling_date"),
        )
    )
    return (
        epi_data.with_row_index("_infection_row")
        .join(sampled_infections, on="_infection_row", how="left")
        .with_columns(pl.col("sequence_id").is_not_null().alias("sampled"))
        .drop("_infection_row", "days", "recovery_days")
    )


def plot_epi_curve(
    epi_data: pl.DataFrame,
    output_path: str = "epi_curve.svg",
) -> tuple[Figure, Axes]:
    daily_counts = epi_data.group_by("date", "sampled").len(name="cases")
    counts_by_date: dict[datetime, dict[bool, int]] = {}
    for date, sampled, cases in daily_counts.iter_rows():
        counts_by_date.setdefault(date, {True: 0, False: 0})[bool(sampled)] += cases

    dates = sorted(counts_by_date)
    sampled_counts = [counts_by_date[date][True] for date in dates]
    unsampled_counts = [counts_by_date[date][False] for date in dates]

    figure, axes = plt.subplots()
    axes.bar(dates, sampled_counts, label="Sampled")
    axes.bar(
        dates,
        unsampled_counts,
        bottom=sampled_counts,
        label="Unsampled",
    )
    axes.set_xlabel("Date")
    axes.set_ylabel("Number of cases infected")
    axes.set_title("Epidemiological curve")
    axes.xaxis.set_major_locator(mdates.AutoDateLocator())
    axes.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(axes.xaxis.get_major_locator())
    )
    axes.legend()
    figure.autofmt_xdate()
    figure.savefig(output_path, format="svg", bbox_inches="tight")
    return figure, axes

if __name__ == "__main__":

    base_date=pl.datetime(2026, 1, 1, 00, 00)

    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", type=str, help="directory for target replicate")

    args=parser.parse_args()

    samples, infections, recoveries=read_e3SIM(
        samples=args.directory + "/sample.csv", 
        infections=args.directory + "/infection.csv", 
        recoveries=args.directory + "/recovery.csv"
        )

    metadata=make_genome_metadata(sample_csv=samples, base_date=base_date)
    print("metadata")
    print(metadata)
    epi=make_epi_case_data(
        infections=infections,
        genome_metadata=metadata,
        recoveries=recoveries,
        )
    print("epi_data")
    print(epi)
    epi.write_csv("epi.csv")

    plot_epi_curve(epi_data=epi)
    plt.show()
    #metadata_double_epoch.write_csv("double_epoch_metadata.csv")
    # print(epi_double_epoch)
    # print(metadata_double_epoch)
    # print(epi_double_epoch.select("sampled").unique())
