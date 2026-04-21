import os
import sys
import polars as pl


def main():
    infile = sys.argv[1]
    output_dir = sys.argv[2]
    # Must include hap in fstring
    fname_hap_fstring = sys.argv[3]

    os.makedirs(output_dir, exist_ok=True)

    df_calls = (
        pl.scan_csv(infile, separator="\t", has_header=True)
        .with_columns(
            hap=pl.col("#chrom").str.extract(r"#(?<hap>.+)#")
        )
        .collect()
    )
    for grp, df_grp in df_calls.partition_by(["hap"], as_dict=True, include_key=False).items():
        hap = grp[0]
        df_grp.write_csv(
            os.path.join(output_dir, f"{fname_hap_fstring.format(hap=hap)}.bed"),
            include_header=True,
            separator="\t"
        )

if __name__ == "__main__":
    raise SystemExit(main())
