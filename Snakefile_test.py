import glob

# Get all BAM stems (without .bam)
bam_stems = [
    f[:-4] for f in glob.glob("/scratch/tsapalou/renamed_bams_7g_T2T/bam_5/*/old_bam/*.bam")
    if f.endswith(".bam")
]

rule all:
    input:
        expand("{path}.bam.bai", path=bam_stems)

rule update_sm_tag:
    input:
        bam = "{path}.bam"
    output:
        bai = "{path}.bam.bai"
    log:
        "{path}.sm_tag_update.log"
    threads: 4
    resources:
        mem_mb = 4000
    shell:
        r"""
        set -euo pipefail
        set -x

        SAMPLE=$(basename $(dirname $(dirname "{input.bam}")))
        BAM="{input.bam}"
        TMP_HEADER=$(mktemp)
        TMP_BAM="${BAM}.tmp"

        OLD_SM=$(samtools view -H "$BAM" | grep "^@RG" | grep -o "SM:[^[:space:]]*" | cut -d':' -f2 | head -n 1 || true)

        echo "Updating $BAM: $OLD_SM → $SAMPLE" > "{log}"

        if [[ -z "$OLD_SM" ]]; then
            echo "No SM tag found in $BAM" >> "{log}"
            exit 1
        fi

        samtools view -H "$BAM" | sed "s/SM:$OLD_SM/SM:$SAMPLE/" > "$TMP_HEADER"
        samtools reheader "$TMP_HEADER" "$BAM" > "$TMP_BAM"
        samtools sort -@ {threads} -m {resources.mem_mb}M -O bam -o "$BAM" "$TMP_BAM"
        samtools index -@ {threads} --output "{output.bai}" "$BAM"

        rm "$TMP_HEADER" "$TMP_BAM"
        """

