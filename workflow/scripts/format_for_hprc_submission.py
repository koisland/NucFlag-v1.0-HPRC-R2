from concurrent.futures._base import Future
from concurrent.futures.process import ProcessPoolExecutor
import sys
import pathlib
import subprocess
import polars as pl

URI = "s3://human-pangenomics/submissions"
SUBM_NAME = "5E9D3123-C9D3-47D4-8BB7-D82FF0DC84EA--HPRC_R2_QC_NUCFLAG"
VERSION = "v1.0.0-a2"

Row = tuple[str, str, int, str, str, str]


def write_md5(file: pathlib.Path, cwd: pathlib.Path) -> None:
    process_checksum = subprocess.run(
        ["md5sum", file.name], check=True, cwd=cwd, capture_output=True, text=True
    )
    with open(f"{file.as_posix()}.md5", "wt") as fh:
        fh.write(process_checksum.stdout)


def link_files(
    input_dir: pathlib.Path,
    output_dir: pathlib.Path,
    sm: str,
    dtype: str,
    hap: int,
    asm_name: str,
) -> list[Row]:
    rows: list[Row] = []
    sm_dir_no_output_dir = pathlib.Path(
        sm, "hprc_r2", "assembly_qc", "nucflag", f"{VERSION}_{dtype}"
    )
    sm_output_dir = output_dir.joinpath(sm_dir_no_output_dir)

    _, hap_name, *_ = asm_name.split("_")
    infile_prefix = f"{sm}_{dtype}_{hap}"
    outfile_prefix = f"{sm}_{hap_name}_NucFlag_{VERSION}_{dtype}"
    print(f"On {asm_name} ({dtype})", file=sys.stderr)
    bw_dir = input_dir.joinpath(f"{infile_prefix}_bw")

    for input_file in bw_dir.glob("*.bw"):
        output_file = sm_output_dir.joinpath(f"{outfile_prefix}_{input_file.stem}.bw")
        assert input_file.exists(), f"Doesn't exist: {input_file.as_posix()}"
        output_file.absolute().symlink_to(input_file.absolute())
        write_md5(output_file, cwd=sm_output_dir)

        rows.append(
            (
                "bw",
                sm,
                hap,
                dtype,
                asm_name,
                f"{URI}/{SUBM_NAME}/{sm_dir_no_output_dir.joinpath(output_file.name).as_posix()}",
            )
        )

    for out_dtype, input_file, output_file in (
        (
            "bed",
            input_dir.joinpath(f"{infile_prefix}_final.bed"),
            sm_output_dir.joinpath(f"{outfile_prefix}.bed"),
        ),
        (
            "plots",
            input_dir.joinpath(f"{infile_prefix}_plots.tar.gz"),
            sm_output_dir.joinpath(f"{outfile_prefix}_plots.tar.gz"),
        ),
        (
            "bed_qv",
            input_dir.joinpath(f"{infile_prefix}.qv.bed"),
            sm_output_dir.joinpath(f"{outfile_prefix}_qv.bed"),
        ),
        (
            "bed_status",
            input_dir.joinpath(f"{infile_prefix}.status.bed"),
            sm_output_dir.joinpath(f"{outfile_prefix}_status.bed"),
        ),
    ):
        assert input_file.exists(), f"Doesn't exist: {input_file.as_posix()}"
        output_file.absolute().symlink_to(input_file.absolute())
        write_md5(output_file, cwd=sm_output_dir)
        rows.append(
            (
                out_dtype,
                sm,
                hap,
                dtype,
                asm_name,
                f"{URI}/{SUBM_NAME}/{sm_dir_no_output_dir.joinpath(output_file.name).as_posix()}",
            )
        )

    return rows


def main():
    # /project/logsdon_shared/projects/HPRC/NucFlag-HPRC-R2/results_hprc/final/{sm}/hprc_r2/assembly_qc/nucflag/{version}_{dtype}
    # /project/logsdon_shared/projects/HPRC/NucFlag-HPRC-R2/results_hprc/final/HG005/hprc_r2/assembly_qc/nucflag/v0.3.3_hifi
    # /project/logsdon_shared/projects/HPRC/NucFlag-HPRC-R2/results_hprc/final/HG005/hprc_r2/assembly_qc/nucflag/v0.3.3_hifi/HG005_mat_hprc_r2_v1.0.1.nucflag.bed
    # /project/logsdon_shared/projects/HPRC/NucFlag-HPRC-R2/results_hprc/final/HG005/hprc_r2/assembly_qc/nucflag/v0.3.3_hifi/HG005_mat_hprc_r2_v1.0.1.nucflag.bed.md5

    input_dir = pathlib.Path(sys.argv[1])
    asm_manifest = sys.argv[2]
    output_dir = pathlib.Path(sys.argv[3])
    output_manifest_dir = pathlib.Path(sys.argv[4])
    processes = int(sys.argv[5])
    df_asm_manifest = pl.read_csv(asm_manifest)

    output_dir.mkdir(exist_ok=True, parents=True)

    # out_dtype,sample_id,haplotype,assembly_name,location
    rows = []

    elems: list[tuple[pathlib.Path, pathlib.Path, str, str, int, str]] = []
    for hifi_file, ont_file in zip(
        input_dir.glob("*_hifi.bed"), input_dir.glob("*_ont.bed"), strict=True
    ):
        sm, _ = hifi_file.stem.split("_")
        for dtype in ("hifi", "ont"):
            sm_dir_no_output_dir = pathlib.Path(
                sm, "hprc_r2", "assembly_qc", "nucflag", f"{VERSION}_{dtype}"
            )
            sm_output_dir = output_dir.joinpath(sm_dir_no_output_dir)
            sm_output_dir.mkdir(exist_ok=True, parents=True)
            for hap in (1, 2):
                asm_name: str = df_asm_manifest.filter(
                    pl.col("sample_id").eq(sm) & pl.col("haplotype").eq(hap)
                )["assembly_name"][0]
                elems.append((input_dir, output_dir, sm, dtype, hap, asm_name))

    with ProcessPoolExecutor(processes) as pool:
        futures: list[Future] = []
        for args in elems:
            future = pool.submit(link_files, *args)
            futures.append(future)
        rows = []
        for future in futures:
            if future.cancelled():
                raise RuntimeError(future.exception())
            else:
                rows.extend(future.result())

    for dtype_out_dtype, df_out_dtype in (
        pl.DataFrame(
            rows,
            schema=[
                "out_dtype",
                "sample_id",
                "haplotype",
                "dtype",
                "assembly_name",
                "location",
            ],
            orient="row",
        )
        .select(
            "sample_id", "haplotype", "assembly_name", "location", "dtype", "out_dtype"
        )
        .partition_by(["dtype", "out_dtype"], include_key=False, as_dict=True)
        .items()
    ):
        dtype, out_dtype = dtype_out_dtype
        filename = output_manifest_dir.joinpath(f"nucflag_{dtype}_{out_dtype}.csv")
        df_out_dtype.sort(by=["haplotype", "sample_id"]).write_csv(filename)


if __name__ == "__main__":
    raise SystemExit(main())
