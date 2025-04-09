import os
import pandas as pd
import shutil
import re
import sys

def process_bam_files(bam_folders, predictions_files, output_dir):
    # Combine predictions files
    all_predictions = []
    for pred_file in predictions_files.split():
        df = pd.read_excel(pred_file)
        all_predictions.append(df)
    predictions_df = pd.concat(all_predictions, ignore_index=True)

    # Make map of -rXX-cXX → sample
    sample_dict = {}
    for _, row in predictions_df.iterrows():
        match = re.search(r'(-r\d+-c\d+)', row['cell'])
        if match:
            sample_dict[match.group(1)] = row['1KG_identified_sample']

    os.makedirs(output_dir, exist_ok=True)

    for bam_folder in bam_folders.split():
        for filename in os.listdir(bam_folder):
            if filename.endswith(".bam") and "bai" not in filename and "raw" not in filename:
                match = re.search(r'(-r\d+-c\d+)', filename)
                if match:
                    short_id = match.group(1)
                    sample = sample_dict.get(short_id)
                    if sample:
                        sample_folder = os.path.join(output_dir, sample)
                        bam_subfolder = os.path.join(sample_folder, "bam")
                        os.makedirs(bam_subfolder, exist_ok=True)

                        # Build full paths
                        src_bam = os.path.join(bam_folder, filename)
                        dest_bam = os.path.join(bam_subfolder, f"{sample}.{filename}")

                        # Skip if already copied
                        if not os.path.exists(dest_bam):
                            shutil.copy2(src_bam, dest_bam)
                            print(f"✅ Copied BAM: {src_bam} → {dest_bam}")
                        else:
                            print(f"⏩ Skipped existing BAM: {dest_bam}")

                        # Copy BAI if exists
                        bai_src = src_bam + ".bai"
                        bai_dest = dest_bam + ".bai"
                        if os.path.exists(bai_src):
                            if not os.path.exists(bai_dest):
                                shutil.copy2(bai_src, bai_dest)
                                print(f"✅ Copied BAI: {bai_src} → {bai_dest}")
                            else:
                                print(f"⏩ Skipped existing BAI: {bai_dest}")
                        else:
                            print(f"⚠️  Missing BAI: {bai_src}")
                    else:
                        print(f"⚠️  No match in predictions for: {short_id}")
                else:
                    print(f"⚠️  No -rXX-cXX pattern in filename: {filename}")

    with open(os.path.join(output_dir, "done.txt"), "w") as f:
        f.write("Done")
    print(f"🏁 Finished! Marker file written at: {os.path.join(output_dir, 'done.txt')}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python process_bam.py <bam_folders> <predictions_files> <output_dir>")
        sys.exit(1)

    bam_folders = sys.argv[1]
    predictions_files = sys.argv[2]
    output_dir = sys.argv[3]

    process_bam_files(bam_folders, predictions_files, output_dir)

