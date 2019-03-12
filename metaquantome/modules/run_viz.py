import os
import subprocess
import json

from metaquantome.util.utils import BASE_DIR
from metaquantome.classes.SampleGroups import SampleGroups


def run_viz(plottype, img, infile, strip=None,
            mode=None, meancol=None, nterms='5', target_rank=None, barcol=6,  # barplot
            textannot=None, fc_name=None, flip_fc=False, gosplit=False,  # volcano
            sinfo=None, filter_to_sig=False, alpha='0.05',  # heatmap
            calculate_sep=False,  # pca
            whichway=None, name=None, id=None, target_onto=None, # ft_dist
            width='5', height='5', tabfile=None):
    """
    Wrapper script for the command-line R visualizations
    The documentation for each of the arguments is in cli.py

    :return: None
    """
    r_script_path = os.path.join(BASE_DIR, 'modules', 'viz.R')
    cmd = ['Rscript', '--vanilla', r_script_path, plottype, img, infile]
    if plottype == "bar":
        cmd += [mode, meancol, nterms, width, height, target_rank, barcol, tabfile]
    elif plottype == "volcano":
        cmd += [str(textannot), fc_name, flip_fc, gosplit, width, height, tabfile]
    elif plottype == "heatmap":
        samp_grps = SampleGroups(sinfo)
        all_intcols_str = ','.join(samp_grps.all_intcols)
        json_dump = json.dumps(samp_grps.sample_names)
        cmd += [all_intcols_str, json_dump, filter_to_sig, alpha, width, height, strip]
    elif plottype == "pca":
        samp_grps = SampleGroups(sinfo)
        all_intcols_str = ','.join(samp_grps.all_intcols)
        json_dump = json.dumps(samp_grps.sample_names)
        cmd += [all_intcols_str, json_dump, calculate_sep, width, height, strip]
    if plottype == "ft_dist":
        cmd += [whichway, name, id, meancol, nterms, width, height,
                target_rank, target_onto, barcol, tabfile]
    else:
        ValueError("Wrong plot type. Must be bar, volcano, heatmap, ft_dist, or pca.")
    # ensure that all elements are strings (even booleans, etc)
    cmd_string = [str(elem) for elem in cmd]

    # run the visualizations, suppressing any output to stdout
    with open(os.devnull, 'w') as fnull:
        subprocess.run(cmd_string, stdout=fnull, check=True)

