
HG002_DATA_MANIFEST = "data/HG002_HPRC_R2_QC_comparison_qc_files.csv"
DF_HG002_QC_DATA = pl.read_csv(HG002_DATA_MANIFEST).with_columns(file_id=pl.col("file") + "_" + pl.col("haplotype"))
DOWNLOAD_FILES = [
    "assembly",
    "hifi_aligned_bam",
    "ont_r9_aligned_bam",
    "ont_r10_aligned_bam",
]
FILE_IDS_BAM = [
    "hifi_aligned_bam_diploid",
    "ont_r9_aligned_bam_diploid",
    "ont_r10_aligned_bam_diploid",
]
SM_WC_MAP = dict(DF_HG002_QC_DATA.filter(pl.col("file").is_in(DOWNLOAD_FILES)).select("file_id", "location").iter_rows())
S3_OUTPUT_DIRS = {
    "hifi_aligned_bam_diploid": "HG002/assemblies/freeze_2/assembly_pipeline/ncbi_upload/assembly_qc/nucflag/v1.0.0a5_hifi/",
    "ont_r9_aligned_bam_diploid": "HG002/assemblies/freeze_2/assembly_pipeline/ncbi_upload/assembly_qc/nucflag/v1.0.0a5_ont_r9/",
    "ont_r10_aligned_bam_diploid": "HG002/assemblies/freeze_2/assembly_pipeline/ncbi_upload/assembly_qc/nucflag/v1.0.0a5_ont_r10/",
}
VERSION = "v1.0.0a5"

rule download_hg002_qc_data:
    output:
        outfile=touch(join(OUTPUT_DIR, "HG002", "data", "{file_id}.done")),
    params:
        uri=lambda wc: SM_WC_MAP[wc.file_id],
        is_bam=lambda wc: wc.file_id.endswith("bam_diploid"),
        output_dir=lambda wc, output: dirname(output.outfile)
    conda:
        "../envs/hg002.yaml"
    shell:
        """
        aws s3 cp {params.uri} {params.output_dir}
        if [[ "{params.is_bam}" == "True" ]]; then
            aws s3 cp {params.uri}.bai {params.output_dir}
        fi
        """


rule merge_asm_hg002_haps:
    input:
        chkpts=expand(rules.download_hg002_qc_data.output, file_id=["assembly_hap_1", "assembly_hap_2"])
    output:
        asm=join(OUTPUT_DIR, "HG002", "data", "HG002_asm_merged.fa.gz"),
        asm_fai=join(OUTPUT_DIR, "HG002", "data", "HG002_asm_merged.fa.gz.fai"),
    params:
        infiles=lambda wc: " ".join([
            join(OUTPUT_DIR, "HG002", "data", basename(SM_WC_MAP["assembly_hap_1"])),
            join(OUTPUT_DIR, "HG002", "data", basename(SM_WC_MAP["assembly_hap_2"]))
        ])
    conda:
        "../envs/hg002.yaml"
    shell:
        """
        zcat {params.infiles} | bgzip > {output.asm}
        samtools faidx {output.asm}
        """


rule run_nucflag_hg002:
    input:
        bam_chkpt=rules.download_hg002_qc_data.output,
        asm=rules.merge_asm_hg002_haps.output.asm,
        asm_fai=rules.merge_asm_hg002_haps.output.asm_fai,
    output:
        calls=join(OUTPUT_DIR, "HG002", "{file_id}.bed"),
        plot_dir=directory(join(OUTPUT_DIR, "HG002", "{file_id}_plot")),
    threads: 12
    log:
        join(LOG_DIR, "HG002", "{file_id}.log")
    benchmark:
        join(BMK_DIR, "HG002", "{file_id}.log")
    conda:
        "../envs/hg002.yaml"
    params:
        bam=lambda wc: join(OUTPUT_DIR, "HG002", "data", basename(SM_WC_MAP[wc.file_id])),
        # hifi, ont_r9, ont_r10
        preset=lambda wc: wc.file_id.replace("_aligned_bam_diploid", ""),
    shell:
        """
        nucflag call \
        -p {threads} \
        -i {params.bam} \
        -f {input.asm} \
        -x {params.preset} \
        -o {output.calls} \
        --output_plot_dir {output.plot_dir} \
        --add_builtin_tracks mapq ident 2> {log}
        """


# Split by haplotype and create directory structure matching HMM-Flagger's
rule prep_data_for_upload:
    input:
        calls=rules.run_nucflag_hg002.output.calls
    output:
        chkpt=touch(join(OUTPUT_DIR, "final_HG002_{file_id}.done"))
    params:
        script="workflow/scripts/format_for_hprc_submission_hg002.py",
        output_dir=lambda wc: join(OUTPUT_DIR, "final_HG002", S3_OUTPUT_DIRS[wc.file_id]),
        fname_fstring=lambda wc: f"HG002_hprc_v2_{wc.file_id}_NucFlag_{VERSION}_hap{{hap}}"
    shell:
        """
        python {params.script} {input.calls} {params.output_dir} {params.fname_fstring}
        """


rule hg002_qc_all:
    input:
        expand(rules.download_hg002_qc_data.output, file_id=SM_WC_MAP.keys()),
        rules.merge_asm_hg002_haps.output,
        expand(
            rules.run_nucflag_hg002.output,
            file_id=FILE_IDS_BAM
        ),
        expand(rules.prep_data_for_upload.output, file_id=FILE_IDS_BAM),
