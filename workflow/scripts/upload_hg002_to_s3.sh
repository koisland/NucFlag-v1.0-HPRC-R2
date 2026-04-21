
# Following flagger's example:
# s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG002/assemblies/freeze_2/assembly_pipeline/ncbi_upload/assembly_qc/hmm_flagger/v1.2.0_hifi/
# So something like this?
# s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG002/assemblies/freeze_2/assembly_pipeline/ncbi_upload/assembly_qc/nucflag/v1.0.0a5_hifi/

set -euo pipefail

SUBMISSION_ID="DC27718F-5F38-43B0-9A78-270F395F13E8"
SUBMISSION_NAME="INT_ASM_PRODUCTION"
INPUT_DIR="/project/logsdon_shared/projects/HPRC/NucFlag-HPRC-R2-v1.0/results/final_HG002"

ssds staging upload \
--submission-id "${SUBMISSION_ID}" \
--name "${SUBMISSION_NAME}" \
"${INPUT_DIR}"
