#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BioProject robustness analysis using manuscript-aligned serotype denominators.

This is the manuscript-aligned implementation. Earlier exploratory variants
are intentionally not part of this release.

Input
-----
合并结果_with_BioProject.xlsx

Required columns
----------------
BioProject
year
Serotype or Serotyope
country   (needed for country × BioProject analyses)

Outputs
-------
05_bioproject_robustness_results.xlsx
05_bioproject_robustness_overview.csv

Analyses
--------
1. Original pooled prevalence.
2. Leave-one-Top10-BioProject-out:
   remove ONE major BioProject at a time.
3. Cumulative exclusion of Top 1 / 3 / 5 / 10 projects.
4. BioProject-balanced:
   calculate prevalence within each BioProject-period first,
   then give each eligible BioProject equal weight.
5. Country × BioProject-balanced:
   each country-BioProject-period cell is one equal-weight sampling unit.
6. Hierarchical country/project-balanced:
   first average eligible BioProjects within each country-period,
   then give countries equal weight.

Periods
-------
2005-2009 / 2010-2014 / 2015-2019 / 2020-2024
"""

import os
import sys
import numpy as np
import pandas as pd

INPUT_FILE = "data/PA_7071_BioProject.xlsx"
OUTPUT_XLSX = "results/BioProject_robustness.xlsx"
OUTPUT_CSV = "results/BioProject_overview.csv"

PERIOD_ORDER = ["2005-2009", "2010-2014", "2015-2019", "2020-2024"]
SEROTYPES = ["O11", "O6"]

# Threshold is the minimum number of serotyped isolates required
# within a BioProject-period or country-BioProject-period sampling unit.
MIN_N_THRESHOLDS = [1, 5, 10, 20]

TOP_N = 10


def fail(msg):
    print("\n[ERROR]", msg)
    sys.exit(1)


def period_from_year(y):
    if pd.isna(y):
        return np.nan
    try:
        y = int(float(y))
    except Exception:
        return np.nan

    if 2005 <= y <= 2009:
        return "2005-2009"
    if 2010 <= y <= 2014:
        return "2010-2014"
    if 2015 <= y <= 2019:
        return "2015-2019"
    if 2020 <= y <= 2024:
        return "2020-2024"
    return np.nan


def clean_string(s):
    s = s.astype("string").str.strip()
    bad = {"", "-", "NA", "N/A", "nan", "None", "<NA>"}
    s = s.mask(s.isin(bad), pd.NA)
    return s


def period_prevalence(data, label, extra=None):
    """Return one row per period for isolate-weighted prevalence."""
    rows = []
    for period in PERIOD_ORDER:
        x = data[data["period"] == period]
        n = len(x)

        row = {
            "analysis": label,
            "period": period,
            "n_isolates": n,
            "n_projects": x["BioProject"].nunique(dropna=True),
            "n_countries": x["country"].nunique(dropna=True),
        }
        if extra:
            row.update(extra)

        for s in SEROTYPES:
            k = int((x["serotype"] == s).sum())
            row[f"{s}_n"] = k
            row[f"{s}_prevalence"] = k / n if n else np.nan

        rows.append(row)

    return pd.DataFrame(rows)


def early_late_from_period_table(table, analysis_name=None, extra=None):
    """Build compact early-vs-late result from a 4-period table."""
    a = table[table["period"].astype(str) == "2005-2009"]
    b = table[table["period"].astype(str) == "2020-2024"]

    if len(a) != 1 or len(b) != 1:
        return None

    a = a.iloc[0]
    b = b.iloc[0]

    row = {
        "analysis": analysis_name if analysis_name is not None else str(a.get("analysis", "")),
        "n_isolates_early": a.get("n_isolates", np.nan),
        "n_isolates_late": b.get("n_isolates", np.nan),
        "n_projects_early": a.get("n_projects", np.nan),
        "n_projects_late": b.get("n_projects", np.nan),
        "n_countries_early": a.get("n_countries", np.nan),
        "n_countries_late": b.get("n_countries", np.nan),
        "O11_early_percent": a["O11_prevalence"] * 100,
        "O11_late_percent": b["O11_prevalence"] * 100,
        "delta_O11_pp": (b["O11_prevalence"] - a["O11_prevalence"]) * 100,
        "O6_early_percent": a["O6_prevalence"] * 100,
        "O6_late_percent": b["O6_prevalence"] * 100,
        "delta_O6_pp": (b["O6_prevalence"] - a["O6_prevalence"]) * 100,
    }

    if extra:
        row.update(extra)

    return row


def make_unit_detail(data, group_cols):
    """
    Create prevalence table for sampling units, e.g.
    BioProject-period or country-BioProject-period.
    """
    g = (
        data.groupby(group_cols, observed=True)
        .size()
        .rename("n")
        .reset_index()
    )

    for s in SEROTYPES:
        hit = (
            data.assign(hit=(data["serotype"] == s).astype(int))
            .groupby(group_cols, observed=True)["hit"]
            .sum()
            .rename(f"{s}_n")
            .reset_index()
        )
        g = g.merge(hit, on=group_cols, how="left")
        g[f"{s}_prevalence"] = g[f"{s}_n"] / g["n"]

    return g


def equal_unit_summary(unit_detail, unit_name, min_n):
    """
    Each eligible unit gets equal weight within each period.
    """
    rows = []
    z0 = unit_detail[unit_detail["n"] >= min_n].copy()

    for period in PERIOD_ORDER:
        z = z0[z0["period"].astype(str) == period].copy()

        row = {
            "analysis": unit_name,
            "min_n": min_n,
            "period": period,
            "n_units": len(z),
            "n_isolates": int(z["n"].sum()) if len(z) else 0,
        }

        if "BioProject" in z.columns:
            row["n_projects"] = z["BioProject"].nunique()
        else:
            row["n_projects"] = np.nan

        if "country" in z.columns:
            row["n_countries"] = z["country"].nunique()
        else:
            row["n_countries"] = np.nan

        for s in SEROTYPES:
            row[f"{s}_prevalence"] = (
                z[f"{s}_prevalence"].mean() if len(z) else np.nan
            )

        rows.append(row)

    return pd.DataFrame(rows)


def hierarchical_country_project_summary(cp_detail, min_n):
    """
    Hierarchical equal weighting:
    1) Within each country-period, average eligible BioProject prevalences.
    2) Across countries within each period, average country-level means.

    This avoids giving countries with many BioProjects proportionally more weight.
    """
    eligible = cp_detail[cp_detail["n"] >= min_n].copy()

    country_period = (
        eligible
        .groupby(["country", "period"], observed=True)
        .agg(
            n_projects=("BioProject", "nunique"),
            n_isolates=("n", "sum"),
            O11_prevalence=("O11_prevalence", "mean"),
            O6_prevalence=("O6_prevalence", "mean"),
        )
        .reset_index()
    )

    rows = []
    for period in PERIOD_ORDER:
        z = country_period[country_period["period"].astype(str) == period]

        rows.append({
            "analysis": "Hierarchical country→BioProject-balanced",
            "min_n": min_n,
            "period": period,
            "n_units": len(z),
            "n_countries": z["country"].nunique(),
            "n_projects": int(z["n_projects"].sum()) if len(z) else 0,
            "n_isolates": int(z["n_isolates"].sum()) if len(z) else 0,
            "O11_prevalence": z["O11_prevalence"].mean() if len(z) else np.nan,
            "O6_prevalence": z["O6_prevalence"].mean() if len(z) else np.nan,
        })

    return pd.DataFrame(rows), country_period


def main():
    print("=" * 82)
    print("BioProject robustness analysis")
    print("=" * 82)

    os.makedirs(os.path.dirname(OUTPUT_XLSX) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(OUTPUT_CSV) or ".", exist_ok=True)

    if not os.path.exists(INPUT_FILE):
        fail(
            f"当前目录没有找到 {INPUT_FILE}\n"
            "请把本脚本与该 Excel 放在同一个目录。"
        )

    df = pd.read_excel(INPUT_FILE, sheet_name=0)
    print("读取文件:", INPUT_FILE)
    print("总行数:", len(df))

    # Robust column detection.
    lower_map = {str(c).strip().lower(): c for c in df.columns}

    year_col = lower_map.get("year")
    country_col = lower_map.get("country")
    bp_col = lower_map.get("bioproject")

    serotype_col = None
    for candidate in ["serotype", "serotyope"]:
        if candidate in lower_map:
            serotype_col = lower_map[candidate]
            break

    if year_col is None:
        fail("没有找到 year 列。")
    if bp_col is None:
        fail("没有找到 BioProject 列。")
    if serotype_col is None:
        fail("没有找到 Serotype / Serotyope 列。")
    if country_col is None:
        fail("没有找到 country 列；country×BioProject 分析需要该列。")

    d = df[[year_col, country_col, bp_col, serotype_col]].copy()
    d.columns = ["year", "country", "BioProject", "serotype"]

    d["country"] = clean_string(d["country"])
    d["BioProject"] = clean_string(d["BioProject"])
    d["serotype"] = clean_string(d["serotype"])

    # Match scripts/run_core_analyses.py. The 2026-10-01 script only
    # excluded blank values and therefore included records typed "NT" or
    # "multiple" in the serotyped denominator. That produces 362/962/2174/
    # 2218 instead of the manuscript-aligned 361/955/2157/2196 counts.
    d["serotype"] = d["serotype"].str.upper()
    invalid_serotypes = {"", "-", "NT", "MULTIPLE", "NAN", "NA", "N/A", "NONE", "<NA>"}
    d["serotype"] = d["serotype"].mask(d["serotype"].isin(invalid_serotypes), pd.NA)

    # For geographic balancing only, do not interpret restricted, missing or
    # not-collected geographical metadata as countries.
    invalid_countries = {"", "-", "NAN", "NA", "N/A", "NONE", "<NA>",
                         "UNKNOWN", "MISSING", "NOT COLLECTED", "RESTRICTED ACCESS"}
    d["country"] = d["country"].mask(d["country"].str.upper().isin(invalid_countries), pd.NA)
    d["period"] = d["year"].apply(period_from_year)

    # Top 10 projects are identified from the entire 7071-row dataset,
    # matching the previous "largest BioProjects" definition.
    full_project_counts = (
        d[d["BioProject"].notna()]
        .groupby("BioProject")
        .size()
        .sort_values(ascending=False)
    )

    top10 = full_project_counts.head(TOP_N)
    top10_projects = top10.index.tolist()

    top10_table = pd.DataFrame({
        "rank": range(1, len(top10) + 1),
        "BioProject": top10.index,
        "n_total_dataset": top10.values,
        "percent_of_7071": top10.values / len(d) * 100,
    })

    print("\nTop 10 BioProjects:")
    print(top10_table.to_string(index=False))

    # Analysis dataset: valid period, serotype, and BioProject.
    temporal = d[
        d["period"].notna()
        & d["serotype"].notna()
        & d["BioProject"].notna()
    ].copy()

    print("\n用于 BioProject temporal analysis:", len(temporal))
    print("BioProjects:", temporal["BioProject"].nunique())
    period_sizes = temporal.groupby("period", observed=True).size().reindex(PERIOD_ORDER, fill_value=0)
    print("Serotyped sample sizes by period:", period_sizes.to_dict())
    if len(d) == 7071:
        expected = [361, 955, 2157, 2196]
        if list(period_sizes) != expected:
            fail(f"Period denominator mismatch against the frozen 7,071-genome manuscript source: {list(period_sizes)} != {expected}")

    # ------------------------------------------------------------
    # 1. Original pooled
    # ------------------------------------------------------------
    pooled_period = period_prevalence(temporal, "Pooled")
    pooled_overview = early_late_from_period_table(
        pooled_period,
        analysis_name="Pooled"
    )

    # ------------------------------------------------------------
    # 2. Leave-one-Top10-out
    # ------------------------------------------------------------
    loo_period_tables = []
    loo_overview_rows = []

    for rank, project in enumerate(top10_projects, start=1):
        x = temporal[temporal["BioProject"] != project].copy()

        pt = period_prevalence(
            x,
            "Leave-one-major-BioProject-out",
            extra={
                "removed_rank": rank,
                "removed_BioProject": project,
                "removed_n_total_dataset": int(full_project_counts.loc[project]),
            }
        )
        loo_period_tables.append(pt)

        row = early_late_from_period_table(
            pt,
            analysis_name=f"Remove {project}",
            extra={
                "removed_rank": rank,
                "removed_BioProject": project,
                "removed_n_total_dataset": int(full_project_counts.loc[project]),
            }
        )
        loo_overview_rows.append(row)

    loo_period = pd.concat(loo_period_tables, ignore_index=True)
    loo_overview = pd.DataFrame(loo_overview_rows)

    # Flag whether direction flips relative to original claimed direction:
    # O11 expected positive; O6 expected negative.
    loo_overview["O11_direction_preserved"] = loo_overview["delta_O11_pp"] > 0
    loo_overview["O6_direction_preserved"] = loo_overview["delta_O6_pp"] < 0
    loo_overview["both_directions_preserved"] = (
        loo_overview["O11_direction_preserved"]
        & loo_overview["O6_direction_preserved"]
    )

    # ------------------------------------------------------------
    # 3. Cumulative Top-N exclusion
    # ------------------------------------------------------------
    cumulative_period_tables = []
    cumulative_overview_rows = []

    for n_remove in [1, 3, 5, 10]:
        remove_set = set(top10_projects[:n_remove])
        x = temporal[~temporal["BioProject"].isin(remove_set)].copy()

        pt = period_prevalence(
            x,
            f"Exclude Top {n_remove}",
            extra={
                "n_removed_projects": n_remove,
                "removed_projects": ";".join(top10_projects[:n_remove]),
            }
        )
        cumulative_period_tables.append(pt)

        row = early_late_from_period_table(
            pt,
            analysis_name=f"Exclude Top {n_remove}",
            extra={
                "n_removed_projects": n_remove,
                "removed_projects": ";".join(top10_projects[:n_remove]),
            }
        )
        cumulative_overview_rows.append(row)

    cumulative_period = pd.concat(cumulative_period_tables, ignore_index=True)
    cumulative_overview = pd.DataFrame(cumulative_overview_rows)

    # ------------------------------------------------------------
    # 4. BioProject-balanced
    # ------------------------------------------------------------
    project_period_detail = make_unit_detail(
        temporal,
        ["BioProject", "period"]
    )

    project_balanced_tables = []
    project_balanced_overview_rows = []

    for min_n in MIN_N_THRESHOLDS:
        pt = equal_unit_summary(
            project_period_detail,
            "BioProject-balanced",
            min_n
        )
        project_balanced_tables.append(pt)

        row = early_late_from_period_table(
            pt,
            analysis_name=f"BioProject-balanced ≥{min_n}",
            extra={"min_n": min_n}
        )
        project_balanced_overview_rows.append(row)

    project_balanced = pd.concat(project_balanced_tables, ignore_index=True)
    project_balanced_overview = pd.DataFrame(project_balanced_overview_rows)

    # ------------------------------------------------------------
    # 5. Country × BioProject-balanced
    # ------------------------------------------------------------
    cp_data = temporal[temporal["country"].notna()].copy()

    country_project_detail = make_unit_detail(
        cp_data,
        ["country", "BioProject", "period"]
    )

    cp_balanced_tables = []
    cp_balanced_overview_rows = []

    for min_n in MIN_N_THRESHOLDS:
        pt = equal_unit_summary(
            country_project_detail,
            "Country×BioProject-unit-balanced",
            min_n
        )
        cp_balanced_tables.append(pt)

        row = early_late_from_period_table(
            pt,
            analysis_name=f"Country×BioProject balanced ≥{min_n}",
            extra={"min_n": min_n}
        )
        cp_balanced_overview_rows.append(row)

    cp_balanced = pd.concat(cp_balanced_tables, ignore_index=True)
    cp_balanced_overview = pd.DataFrame(cp_balanced_overview_rows)

    # ------------------------------------------------------------
    # 6. Hierarchical country → BioProject-balanced
    # ------------------------------------------------------------
    hierarchical_tables = []
    hierarchical_country_tables = []
    hierarchical_overview_rows = []

    for min_n in MIN_N_THRESHOLDS:
        pt, country_period = hierarchical_country_project_summary(
            country_project_detail,
            min_n
        )
        hierarchical_tables.append(pt)

        country_period = country_period.copy()
        country_period["min_n"] = min_n
        hierarchical_country_tables.append(country_period)

        row = early_late_from_period_table(
            pt,
            analysis_name=f"Hierarchical country→BioProject ≥{min_n}",
            extra={"min_n": min_n}
        )
        hierarchical_overview_rows.append(row)

    hierarchical = pd.concat(hierarchical_tables, ignore_index=True)
    hierarchical_country_detail = pd.concat(
        hierarchical_country_tables,
        ignore_index=True
    )
    hierarchical_overview = pd.DataFrame(hierarchical_overview_rows)

    # ------------------------------------------------------------
    # Combine overview
    # ------------------------------------------------------------
    overview_parts = [
        pd.DataFrame([pooled_overview]),
        cumulative_overview,
        project_balanced_overview,
        cp_balanced_overview,
        hierarchical_overview,
    ]

    overview = pd.concat(overview_parts, ignore_index=True, sort=False)

    # Add direction interpretation.
    overview["O11_direction"] = np.where(
        overview["delta_O11_pp"] > 0, "increase",
        np.where(overview["delta_O11_pp"] < 0, "decrease", "no change")
    )
    overview["O6_direction"] = np.where(
        overview["delta_O6_pp"] > 0, "increase",
        np.where(overview["delta_O6_pp"] < 0, "decrease", "no change")
    )
    overview["original_direction_preserved"] = (
        (overview["delta_O11_pp"] > 0)
        & (overview["delta_O6_pp"] < 0)
    )

    # ------------------------------------------------------------
    # Project composition by period and serotype for Top 10
    # ------------------------------------------------------------
    top10_composition_rows = []

    for rank, project in enumerate(top10_projects, start=1):
        x = temporal[temporal["BioProject"] == project]

        for period in PERIOD_ORDER:
            z = x[x["period"] == period]
            n = len(z)

            top10_composition_rows.append({
                "rank": rank,
                "BioProject": project,
                "period": period,
                "n_serotyped": n,
                "O11_n": int((z["serotype"] == "O11").sum()),
                "O11_percent": (
                    (z["serotype"] == "O11").mean() * 100 if n else np.nan
                ),
                "O6_n": int((z["serotype"] == "O6").sum()),
                "O6_percent": (
                    (z["serotype"] == "O6").mean() * 100 if n else np.nan
                ),
            })

    top10_composition = pd.DataFrame(top10_composition_rows)

    # ------------------------------------------------------------
    # Print the decisive results
    # ------------------------------------------------------------
    print("\n" + "=" * 82)
    print("A. Leave-one-Top10-BioProject-out")
    print("=" * 82)

    print(
        loo_overview[
            [
                "removed_rank",
                "removed_BioProject",
                "removed_n_total_dataset",
                "O11_early_percent",
                "O11_late_percent",
                "delta_O11_pp",
                "O6_early_percent",
                "O6_late_percent",
                "delta_O6_pp",
                "both_directions_preserved",
            ]
        ]
        .round(2)
        .to_string(index=False)
    )

    print("\n单个项目删除后仍同时保持 O11↑ / O6↓：",
          int(loo_overview["both_directions_preserved"].sum()),
          "/",
          len(loo_overview))

    print("\n" + "=" * 82)
    print("B. Balanced analyses: early vs late")
    print("=" * 82)

    show_cols = [
        "analysis",
        "n_projects_early",
        "n_projects_late",
        "n_countries_early",
        "n_countries_late",
        "O11_early_percent",
        "O11_late_percent",
        "delta_O11_pp",
        "O6_early_percent",
        "O6_late_percent",
        "delta_O6_pp",
        "original_direction_preserved",
    ]

    print(
        overview[show_cols]
        .round(2)
        .to_string(index=False)
    )

    # ------------------------------------------------------------
    # Export
    # ------------------------------------------------------------
    try:
        writer = pd.ExcelWriter(OUTPUT_XLSX, engine="xlsxwriter")
    except ModuleNotFoundError:
        fail(
            "没有安装 xlsxwriter。\n"
            "请运行：conda install -c conda-forge xlsxwriter\n"
            "然后重新运行本脚本。"
        )

    with writer:
        overview.to_excel(writer, sheet_name="overview", index=False)
        loo_overview.to_excel(writer, sheet_name="leave_one_out", index=False)
        loo_period.to_excel(writer, sheet_name="LOO_four_periods", index=False)

        cumulative_overview.to_excel(
            writer, sheet_name="cumulative_exclusion", index=False
        )
        cumulative_period.to_excel(
            writer, sheet_name="cumulative_4periods", index=False
        )

        project_balanced_overview.to_excel(
            writer, sheet_name="project_bal_overview", index=False
        )
        project_balanced.to_excel(
            writer, sheet_name="project_bal_4periods", index=False
        )
        project_period_detail.to_excel(
            writer, sheet_name="project_period_detail", index=False
        )

        cp_balanced_overview.to_excel(
            writer, sheet_name="country_project_overview", index=False
        )
        cp_balanced.to_excel(
            writer, sheet_name="country_project_4periods", index=False
        )
        country_project_detail.to_excel(
            writer, sheet_name="country_project_detail", index=False
        )

        hierarchical_overview.to_excel(
            writer, sheet_name="hierarchical_overview", index=False
        )
        hierarchical.to_excel(
            writer, sheet_name="hierarchical_4periods", index=False
        )
        hierarchical_country_detail.to_excel(
            writer, sheet_name="hierarchical_country", index=False
        )

        top10_table.to_excel(
            writer, sheet_name="Top10_projects", index=False
        )
        top10_composition.to_excel(
            writer, sheet_name="Top10_composition", index=False
        )
        pooled_period.to_excel(
            writer, sheet_name="pooled_4periods", index=False
        )

        workbook = writer.book

        header_fmt = workbook.add_format({
            "bold": True,
            "bg_color": "#D9EAF7",
            "border": 1
        })

        percent_number_fmt = workbook.add_format({"num_format": "0.00"})
        bool_fmt = workbook.add_format({"align": "center"})

        for sheet_name in writer.sheets:
            ws = writer.sheets[sheet_name]
            ws.freeze_panes(1, 0)

            # Retrieve corresponding dataframe approximately for autofilter range.
            # Adding autofilter on a generous used area is harmless, but
            # we only use it on the explicit written data when possible.
            ws.set_default_row(16)

        sheet_df_pairs = {
            "overview": overview,
            "leave_one_out": loo_overview,
            "LOO_four_periods": loo_period,
            "cumulative_exclusion": cumulative_overview,
            "cumulative_4periods": cumulative_period,
            "project_bal_overview": project_balanced_overview,
            "project_bal_4periods": project_balanced,
            "project_period_detail": project_period_detail,
            "country_project_overview": cp_balanced_overview,
            "country_project_4periods": cp_balanced,
            "country_project_detail": country_project_detail,
            "hierarchical_overview": hierarchical_overview,
            "hierarchical_4periods": hierarchical,
            "hierarchical_country": hierarchical_country_detail,
            "Top10_projects": top10_table,
            "Top10_composition": top10_composition,
            "pooled_4periods": pooled_period,
        }

        for sheet_name, data in sheet_df_pairs.items():
            ws = writer.sheets[sheet_name]

            if len(data.columns) > 0:
                ws.autofilter(0, 0, max(len(data), 1), len(data.columns) - 1)

            for col_idx, col in enumerate(data.columns):
                width = max(12, min(34, len(str(col)) + 3))

                # Widen fields that contain project lists or analysis names.
                if col in ["analysis", "removed_projects"]:
                    width = 34
                if col in ["BioProject", "removed_BioProject"]:
                    width = 20

                ws.set_column(col_idx, col_idx, width)

                # Rewrite header with consistent format.
                ws.write(0, col_idx, col, header_fmt)

                if (
                    "_percent" in str(col)
                    or "_pp" in str(col)
                    or str(col) == "percent_of_7071"
                ):
                    ws.set_column(col_idx, col_idx, max(width, 17), percent_number_fmt)

                if "preserved" in str(col):
                    ws.set_column(col_idx, col_idx, max(width, 18), bool_fmt)

    overview.to_csv(OUTPUT_CSV, index=False)

    print("\n" + "=" * 82)
    print("输出完成")
    print("=" * 82)
    print(" ", OUTPUT_XLSX)
    print(" ", OUTPUT_CSV)
    print("\nBioProject sensitivity analysis completed.")
    print("Review the overview CSV and workbook for scenario-specific changes in direction.")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="BioProject sensitivity analyses (aligned to manuscript serotype coding)")
    p.add_argument("--input", default=INPUT_FILE, help="Excel with BioProject/year/country/Serotype columns")
    p.add_argument("--output-xlsx", default=OUTPUT_XLSX)
    p.add_argument("--output-csv", default=OUTPUT_CSV)
    args = p.parse_args()
    INPUT_FILE, OUTPUT_XLSX, OUTPUT_CSV = args.input, args.output_xlsx, args.output_csv
    main()
