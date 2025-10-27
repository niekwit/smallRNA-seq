# import basic packages
import glob
import pandas as pd
#from snakemake.utils import validate


# validate sample sheet and config file
#validate(samples, schema="../../config/schemas/samples.schema.yaml")
#validate(config, schema="../../config/schemas/config.schema.yaml")

def genome_url(genome):
    return f"https://www.dropbox.com/scl/fo/jdpgn6fl8ngd3th3zebap/AEEuUWqO7IggpwHwACz3uEk/TEsmall/{genome}.tar.gz?rlkey=41oz6ppggy82uha5i3yo1rnlx&dl=0"


def samples():

    fastq = glob.glob("reads/*.fastq.gz")
    return [f.split("/")[-1].replace(".fastq.gz", "") for f in fastq]
