#!/usr/bin/env python3
"""Build the Data S1 supplement from the analysis inputs."""

from pathlib import Path
from subprocess import check_output
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile


ROOT = Path(__file__).resolve().parents[1]
INPUTS = (
    "data/Peppe_2011_calibration_data_leaf_level_clean.csv",
    "data/Peppe_2011_fossil_data_June_2026_leaf_level_clean.csv",
    "data/extra_calibration_data_for_LMA.csv",
    "data/best_wcvp.tre_dated",
    "data/published_dilp_fossil.csv",
)


def build_archive():
    destination = ROOT / "submission/submission_1/Data_S1.zip"
    ages = check_output(
        [
            "Rscript", "--vanilla", "-e",
            'source("code/fossil_taxonomy.R"); '
            'ages <- sort(fossil_site_ages(), decreasing = TRUE); '
            'write.csv(data.frame(site = names(ages), age_ma = unname(ages)), '
            'stdout(), row.names = FALSE)',
        ],
        cwd=ROOT,
    )
    with TemporaryDirectory(dir=destination.parent) as directory:
        temporary = Path(directory) / destination.name
        with ZipFile(temporary, "w", compression=ZIP_DEFLATED) as archive:
            for path in INPUTS:
                archive.write(ROOT / path, path)
            archive.write(ROOT / "submission/Data_S1_README.md", "Data_S1_README.md")
            archive.writestr("data/fossil_site_ages.csv", ages)
        temporary.replace(destination)
    print(destination)


if __name__ == "__main__":
    build_archive()
