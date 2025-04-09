import glob
import os

# Get all BAMs inside old_bam folders
bam_files = glob.glob("/scratch/tsapalou/renamed_bams_7g_T2T/bam_5/*/old_bam/*.bam")

# Remove trailing .bam for use as a wildcard in the rule
bam_stems = [f[:-4] for f in bam_files]  # removes ".bam"

rule all:
    input:
        expand("{path}.sm_tagged.ok", path=bam_stems)

rule update_sm_tag:
    input:
        bam="{path}.bam"
    output:
        done="{path}.sm_tagged.ok"
    threads: 4
    shell:
        """
        set -euo pipefail
        module load SAMtools/1.14-GCC-11.2.0

        BAM_FILE="{input.bam}"
        SAMPLE=$(basename $(dirname $(dirname "$BAM_FILE")))
        TEMP_HEADER=$(mktemp)

        OLD_SM=$(samtools view -H "$BAM_FILE" | grep "^@RG" | grep -o "SM:[^[:space:]]*" | cut -d':' -f2)

        echo "Updating $BAM_FILE: $OLD_SM → $SAMPLE"

        samtools view -H "$BAM_FILE" | sed "s/SM:$OLD_SM/SM:$SAMPLE/" > "$TEMP_HEADER"
        samtools reheader "$TEMP_HEADER" "$BAM_FILE" | samtools sort -o "${BAM_FILE}.tmp"
        mv "${BAM_FILE}.tmp" "$BAM_FILE"
        samtools index "$BAM_FILE"
        rm "$TEMP_HEADER"

        # Create and remove marker so Snakemake thinks job succeeded
        touch {output.done}
        rm {output.done}
        """
