import ftplib
import Bio.ExPASy.Enzyme as Enz
import os
import json
import re
import pandas as pd
import logging
import sys


# TODO:
# move logging to runner
# unify with GO term download
# add handling if description is not found in database
# more tests
# move rel abundance to stats.py, generify


LEVEL_NAMES = ['ec0', 'ec1', 'ec2', 'ec3']
ALL_UNKNOWN = ['-']*4


def enzyme_database_handler(overwrite, data_dir):
    logging.basicConfig(level=logging.DEBUG, format='%(message)s', stream = sys.stderr)

    dat_path = os.path.join(data_dir, 'enzyme.dat')
    class_path = os.path.join(data_dir, 'enzclass.txt')

    if (os.path.exists(dat_path) or os.path.exists(class_path)) and not overwrite:
        logging.info('Using ENZYME files in ' + data_dir)
    else:
        logging.info('Downloading enzyme files from ftp.expasy.org to ' + data_dir)
        with ftplib.FTP('ftp.expasy.org') as ftp:
            ftp.login()
            ftp.cwd('/databases/enzyme')

            dat_file = open(dat_path, 'wb')
            class_file = open(class_path, 'wb')
            ftp.retrbinary('RETR enzyme.dat', dat_file.write)
            ftp.retrbinary('RETR enzclass.txt', class_file.write)

            dat_file.close()
            class_file.close()

            ec_path = os.path.join(data_dir, 'ec_id.json')

            # create a simpler file from enzyme.dat
            create_ec_num_enzyme_name_association_file(dat_path, ec_path)

            # get clean enzclass.txt and dump to json
            read_enzyme_class_to_json(data_dir)

    # load both and return
    return load_combined_enzyme_class_ec_id(data_dir)


def create_ec_num_enzyme_name_association_file(enzdat_file, ec_id_file):
    with open(enzdat_file) as file:
        enzyme_dat = Enz.parse(file)
        enz_id = dict()
        for record in enzyme_dat:
            enz_id[record['ID']] = record['DE']

    with open(ec_id_file, 'w') as ec:
        json.dump(enz_id, ec)


def read_enzyme_class_to_json(enzyme_dir):
    class_path = os.path.join(enzyme_dir, 'enzclass.txt')
    enzyme_class = dict()
    with open(class_path) as f:
        for line in f:
            class_match = re.compile(r"^\d\.")
            id_pattern = re.compile(r"(^\d\.[ \a\d-]*\.[ \a\d-]*\.[ \a\d-]*)")
            descript_pattern = re.compile(r"([A-Z]{1}.*)$")
            if re.match(pattern=class_match, string=line[0:2]):
                id=re.search(id_pattern, line)
                clean_id = re.sub(pattern=r' ', repl='', string=id.group(1))
                descript=re.search(descript_pattern, line).group(1)
                enzyme_class[clean_id] = descript

    enz_json_path = os.path.join(enzyme_dir, 'enzclass.json')
    with open(enz_json_path, 'w') as enz_json:
        json.dump(enzyme_class, enz_json)


def load_combined_enzyme_class_ec_id(enzyme_dir):
    enz_class = load_enzyme_class(enzyme_dir)
    ec_id = load_ec_id(enzyme_dir)

    both = {**enz_class, **ec_id}

    return both


def load_enzyme_class(enzyme_dir):
    enz_json_path = os.path.join(enzyme_dir, 'enzclass.json')
    with open(enz_json_path, 'r') as enz_json:
        enz_class = json.load(enz_json)
    return enz_class


def load_ec_id(enzyme_dir):
    ec_path = os.path.join(enzyme_dir, 'ec_id.json')
    with open(ec_path) as ec:
        ec_id = json.load(ec)
    return ec_id


def expand_ec(ecid):
    split_ec = ecid.split('.')
    deepest_known = 0
    for i,j in zip(split_ec, ALL_UNKNOWN):
        if i != j:
            deepest_known += 1
    ec_dict = {LEVEL_NAMES[i]: '.'.join(split_ec[0:(i+1)] + ALL_UNKNOWN[(i+1):]) for i in range(deepest_known)}
    return pd.Series(ec_dict)


def abundance_level(df, level, all_intcols, norm_to_rank=False):
    # sum intensities in each rank
    summed_abund = df.groupby(by=level)[all_intcols].sum(axis=0)

    # normalize to each sample - not currently done
    if norm_to_rank:
        return_df = summed_abund / summed_abund.sum(axis=0)
    else:
        return_df = summed_abund

    return_df['level'] = level
    return_df['id'] = return_df.index
    return return_df