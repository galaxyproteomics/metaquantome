import os
import subprocess
import json

from metaquantome.util.utils import BASE_DIR
from metaquantome.SampleGroups import SampleGroups


def run_viz(plottype, img, infile,
            mode=None, meancol=None, nterms='5', target_rank=None,
            textannot=None, fc_name=None, gosplit=False,
            sinfo=None, filter_to_sig=False, alpha='0.05',
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
    if plottype == "heatmap":
        samp_grps = SampleGroups(sinfo)
        all_intcols_str = ','.join(samp_grps.all_intcols)
        json_dump = json.dumps(samp_grps.sample_names)
        subprocess.run(['Rscript',
                        r_script_path,
                        plottype,
                        img,
                        infile,
                        all_intcols_str,
                        json_dump,
                        str(filter_to_sig),
                        alpha,
                        width,
                        height])
    return 0
