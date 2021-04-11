import os
import sys
import glob
import shutil

from snakemake import snakemake

# arguments order
# 

def get_snakefile(root_dir = ".", file = "Snakefile"):
    sf = os.path.join(root_dir, file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file;  tried %s" %sf)
    return sf

def run_snakemake(Snakefile, configfiles, cores):
    status = snakemake(snakefile= Snakefile, 
                      configfiles = configfiles, 
                      cores = cores,
                      force_incomplete = True)
    return status

# 
def main(root_dir, args):
    if len(args) < 3:
        print("python run.py unfinished config.yaml cores")
        sys.exit(1)

    unfinished_dir = args[1]
    config_file    = [args[2]]
    cores = 79 if len(args) <= 3 else int(args[3])

    if not os.path.exists("metadata"):
        os.makedirs("metadata")

    sf = get_snakefile(root_dir, "Snakefile")
    sample_files = glob.glob( os.path.join(unfinished_dir,  "*.txt") )
    while len(sample_files) > 0 :
        for i in range(min(10, len(sample_files))):
            file = sample_files.pop()
            shutil.move(file, "metadata")
            status = run_snakemake(sf, config_file, cores )
            if not status:
                sys.exit(1)

if __name__ == '__main__':
    root_dir = os.path.dirname(os.path.abspath(__file__))
    args = sys.argv
    main(root_dir, args) 


