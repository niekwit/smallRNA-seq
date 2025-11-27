import os
from snakemake.shell import shell

# Get current directory to return later
current_dir = os.getcwd()

count_file = os.path.join(current_dir, snakemake.output["txt"])
report = os.path.join(current_dir, snakemake.output["report"])
log_file = os.path.join(current_dir, snakemake.log[0])
fastq = os.path.basename(snakemake.input["fastq"])

outdir = os.path.join(current_dir, snakemake.params["outdir"])
os.makedirs(outdir, exist_ok=True)

fastq_dir = os.path.dirname(snakemake.input["fastq"])

os.chdir(fastq_dir)

print(f"Running TEsmall for sample:", snakemake.wildcards["sample"])

shell(
    f"TEsmall -a {snakemake.params.adapter} "
    f"--minlen {snakemake.params.minlen} "
    f"--maxlen {snakemake.params.maxlen} "
    f"--genome {snakemake.params.genome} "
    f"--maxaln {snakemake.params.maxaln} "
    f"--mismatch {snakemake.params.mismatch} "
    f"--parallel {snakemake.threads} "
    f"--fastq {fastq} "
    f"--label {snakemake.wildcards.sample} "
    f"--dbfolder {snakemake.params.dbfolder} "
)

print(f"Moving output files to final destination {outdir}")

# Move all contents of current dir to outdir, except the original fastq file (use python)
for file in os.listdir("."):
    if file != fastq:
        shell(f"mv {file} {outdir}/{file}")

print(f"{snakemake.wildcards['sample']} done!")