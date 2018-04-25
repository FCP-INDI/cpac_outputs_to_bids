#!/usr/bin/env python
import argparse
import os
import yaml
import subprocess
import json
import sys

import datetime
import time
import urllib.request

import cpac_output_to_bids as cpb


def main(input_arguments=None):
    """
    Provides the command line interface for the cpac_output_to_bids package.

    :param input_arguments: command line arguments, detailed below
    :return: 1 on failure, 0 otherwise
    """
    """
    functionality for command line interface, 

    :param input_arguments: command line argument
    :return: 0 on success, 1 on error
    """

    parser = argparse.ArgumentParser(description='CPAC output to bids {0} Command Line Tool'.format(cpb.__version__))

    parser.add_argument('cpac_data_config',
                        help='CPAC data configuration file used to generate output, helps map derivatives to the '
                             'original data')

    parser.add_argument('cpac_output_path',
                        help='This is the directory of CPAC outputs that will be organized.')

    parser.add_argument('output_dir',
                        help='The directory where the output files should be stored. Output files will be arranged '
                             'according to the BIDS derivative standard in a /derivative subdirectory.')

    parser.add_argument('command',
                        help='Type of renaming to be performed, many of the text files such as nuisance regressors.1d '
                             'must be rewritten and will be unaffected by this command. Imaging files and other files ' 
                             'that do not need to be rewritten can be renamed using either "copy", which will create a '
                             'new copy of the files or "link", which will create a hard link for imaging files. '
                             'Creating hard links requires that the cpac_dir and output_dir are on the same block '
                             'device. "dry_run" will go through and calculate all of the conversions, and write out '
                             'results along with pertinent config files, but will not actually generate or rename any '
                             'files.',
                        choices=['copy', 'link', 'dry_run'])

    parser.add_argument('--derivative_conversion_config',
                        help='Path to a JSON file that contains a list of derivative identifiers and whether they '
                             'should be included in the conversion or not. The "dry_run" command will write out a '
                             'default file that can be modified.')

    parser.add_argument('--debug',
                        help='Print out additional information for debugging.',
                        action='store_true')

    # get the command line arguments
    args = parser.parse_args(input_arguments)

    # create a timestamp for writing config files
    timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d%H%M%S')

    # read in the data configuration file
    with open(args.cpac_data_config, 'r') as yaml_input_stream:
        cpac_data_config_list = yaml.load(yaml_input_stream)

    cpac_data_config_dict = {}
    for cpac_data_config in cpac_data_config_list:
        cpac_data_config_dict[
            "_".join(
                [cpac_data_config['subject_id'], cpac_data_config['unique_id']]).lower()] = cpac_data_config

    # get a list of paths from the directory
    cpac_output_file_paths = []
    for path_root, directories, file_names in os.walk(args.cpac_output_path):
        for file_name in file_names:
            if file_name[0] != '.':
                cpac_output_file_paths.append(os.path.join(path_root, file_name))

    # get the list of derivatives to write:
    if not args.derivative_conversion_config:
        args.derivative_conversion_config = os.path.join(os.path.dirname(__file__),
                                                         'config/cpac_derivatives_dictionary.json')

    with open(args.derivative_conversion_config, 'r') as json_input_stream:
        cpac_derivatives_dictionary = json.load(json_input_stream)

    cpac_write_derivative_list = [key for (key, value) in cpac_derivatives_dictionary.items() if
                                  value['write_derivative'] is True]

    bids_dictionary = cpb.extract_bids_derivative_info_from_cpac_outputs(cpac_output_file_paths,
                                                                         cpac_data_config_dict, {},
                                                                         cpac_write_derivative_list)

    # write out all of the conversion information to json files
    debug_data_config_filename = '/tmp/data_config_{0}.json'.format(timestamp)
    with open(debug_data_config_filename, 'w') as json_out:
        json.dump(cpac_data_config_dict, json_out, indent=4)

    debug_cpac_output_file_paths_filename = '/tmp/cpac_outputs_{0}.json'.format(timestamp)
    with open(debug_cpac_output_file_paths_filename, 'w') as json_out:
        json.dump(cpac_output_file_paths, json_out, indent=4)

    debug_derivative_conversion_config_filename = '/tmp/derivatives_config_{0}.json'.format(timestamp)
    with open(debug_derivative_conversion_config_filename, 'w') as json_out:
        json.dump(cpac_derivatives_dictionary, json_out, indent=4)

    print("cpac_output_to_bids configuration files written to: {0}, {1}".format(debug_data_config_filename,
                                                                                debug_cpac_output_file_paths_filename))
    return 0