import os
import subprocess

from metaquantome.util.utils import BASE_DIR


def run_viz(plottype, img, infile,
            mode=None, meancol=None, nterms='5', target_rank=None,
            textannot=None, fc_name=None, gosplit=False,
            width='5', height='5'):
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
    if plottype == "volcano":
        subprocess.run(['Rscript',
                        r_script_path,
                        plottype,
                        img,
                        infile,
                        str(textannot),
                        fc_name,
                        str(gosplit),
                        width,
                        height])
    return 0
