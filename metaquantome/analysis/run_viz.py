import os
import subprocess

from metaquantome.util.utils import BASE_DIR


def run_viz(plottype, img, infile, mode,
            meancol, nterms='5', width='5',
            height='5',
            target_rank=None):
    r_script_path = os.path.join(BASE_DIR, 'analysis', 'viz.R')
    if plottype == "bar":
        subprocess.run(['Rscript',
                        r_script_path,
                        plottype,
                        img,
                        infile,
                        mode,
                        meancol,
                        nterms,
                        width,
                        height,
                        target_rank])
    return 0
