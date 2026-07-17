import os
import subprocess
import json
import pandas as pd
from multiprocessing import Pool

# --- CONFIGURATION ---
aln_dir = "/home/insan/CompGenomics/dnds/new_coding_alignment_26_Jan"
tree_file = "/home/insan/CompGenomics/dnds/species_tree.txt"
output_dir = "/home/insan/CompGenomics/dnds/hyphy_meme_results"
final_csv = "HyPhy_MEME_final.csv"

# MEME is slow, so we run fewer OGs at once but give them more cores
threads = 12
cpus_per_run = 10

os.makedirs(output_dir, exist_ok=True)


def run_meme_task(fasta_file):
    og_id = fasta_file.replace(".fasta", "")
    fasta_path = os.path.join(aln_dir, fasta_file)
    json_out = os.path.join(output_dir, f"{og_id}.json")

    # CLI command for MEME
    cmd = [
        "hyphy", "meme",
        "--alignment", fasta_path,
        "--tree", tree_file,
        "--output", json_out
    ]

    env = os.environ.copy()
    env["HYPHY_NUM_CPUS"] = str(cpus_per_run)

    try:
        subprocess.run(cmd, env=env, capture_output=True, check=True)

        with open(json_out, "r") as f:
            data = json.load(f)

        # MEME results are stored per site/codon.
        # We count how many sites have p-value <= 0.05.
        sig_sites = []

        # 'MLE' contains the results for each site
        content = data.get("MLE", {}).get("content", {}).get("0", [])
        headers = data.get("MLE", {}).get("headers", [])

        # Find the index for p-value in the JSON headers
        p_idx = -1
        for i, h in enumerate(headers):
            if h[0] == "p-value":
                p_idx = i

        for site_idx, site_data in enumerate(content):
            p_val = site_data[p_idx]

            if p_val <= 0.05:
                # +1 for 1-based numbering
                sig_sites.append(str(site_idx + 1))

        return {
            "OG_ID": og_id,
            "Num_Sig_Sites": len(sig_sites),
            "Significant_Sites": ", ".join(sig_sites) if sig_sites else "None",
            "Status": "Success"
        }

    except Exception as e:
        return {
            "OG_ID": og_id,
            "Status": f"Error: {str(e)}"
        }


if __name__ == "__main__":
    files = [f for f in os.listdir(aln_dir) if f.endswith(".fasta")]

    print(f"Starting HyPhy MEME Batch Run: {len(files)} files.")

    with Pool(processes=threads) as pool:
        final_results = pool.map(run_meme_task, files)

    df = pd.DataFrame(final_results)
    df.to_csv(final_csv, index=False)

    print(f"MEME Analysis complete! Results saved to {final_csv}")
