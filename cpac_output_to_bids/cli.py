#!/usr/bin/env python
import argparse
import os
import yaml
import json
import sys

import datetime
import time
# import urllib.request

import cpac_output_to_bids as cpb
from . import __version__ as VERSION

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

    parser.add_argument('--derivative_conversion_config',
                        help='Path to a JSON file that contains a list of derivative identifiers and whether they '
                             'should be included in the conversion or not. The "dry_run" command will write out a '
                             'default file that can be modified.')

    parser.add_argument('--debug',
                        help='Print out additional information for debugging.',
                        action='store_true')

    parser.add_argument('-c', '--cpac_output_path',
                        help='This is the directory of CPAC outputs that will be converted to bids. Required.',
                        required=True)

    parser.add_argument('-d', '--cpac_data_config',
                        help='CPAC data configuration file used to generate output, helps map derivatives to the '
                             'original data. Required.',
                        required=True)

    parser.add_argument('-o', '--output_dir',
                        help='The directory where the output files should be stored. Output files will be arranged '
                             'according to the BIDS derivative standard in a /derivative subdirectory. If this is '
                             'omitted the data will be organized into the bids directory containing the original '
                             'raw data as determined from the cpac data config')

    parser.add_argument('command',
                        help='Type of renaming to be performed, many of the text files such as nuisance regressors.1d '
                             'must be rewritten and will be unaffected by this command. Imaging files and other files ' 
                             'that do not need to be rewritten can be renamed by: "copy": create a new copy of the '
                             'files. "hard_link": the new file will be a hard link to the original file. This method '
                             'avoids duplication on your hard disk and the file will remain accessible if the original '
                             'CPAC output is deleted. But, both the hard link and cpac output have to be on the same '
                             'physical disk. "sym_link": the new file is a symbolic link to the original file. Like '
                             'hard links, this method avoids duplicating data on the hard drive. But symbolic links '
                             'can span storage devices. The disadvantage is that symbolic links will no longer work if '
                             'the original cpac output is deleted.\n"dry_run" will go through and calculate all of the '
                             'conversions, and write out results along with pertinent config files, but will not '
                             'actually generate or rename any files.',
                        choices=['copy', 'hard_link', 'sym_link', 'dry_run'])

    # get the command line arguments
    args = parser.parse_args(input_arguments)

    debug_flag = False
    if args.debug:
        debug_flag = True
    if 'dry_run' in args.command:
        debug_flag = True

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
                                                                         cpac_write_derivative_list,
                                                                         args.output_dir)

    cpac_bids_path_mapping = cpb.create_bids_outputs_dictionary(bids_dictionary)

    if debug_flag is True:
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

        debug_cpac_bids_mapping_filename = '/tmp/cpac_bids_mapping_{0}.json'.format(timestamp)
        with open(debug_cpac_bids_mapping_filename, 'w') as json_out:
            json.dump(cpac_bids_path_mapping, json_out, indent=4)

        print("\ncpac_output_to_bids configuration files written to:\n\t{0}\n\t{1}\n\t{2}\n\t{3}\n".format(
            debug_data_config_filename, debug_cpac_output_file_paths_filename,
            debug_derivative_conversion_config_filename,
            debug_cpac_bids_mapping_filename))

    if 'dry_run' in args.command:
        print("This is a dry run, none of the cpac outputs were converted to bids format.")
    else:
        number_of_files_converted = 0
        if 'copy' in args.command:
            print("Copying bidsified files into {}.".format(args.output_dir))
            number_of_files_converted = cpb.copy_cpac_outputs_into_bids(cpac_bids_path_mapping,
                                                                        hard_link_instead_of_copy=False,
                                                                        soft_link_instead_of_copy=False)
        elif 'sym_link' in args.command:
            print("Creating symbolic links for bidsified files in {}.".format(args.output_dir))
            number_of_files_converted = cpb.copy_cpac_outputs_into_bids(cpac_bids_path_mapping,
                                                                        hard_link_instead_of_copy=False,
                                                                        soft_link_instead_of_copy=True)
        elif 'hard_link' in args.command:
            print("Creating hard links for bidsified files in {}.".format(args.output_dir))
            number_of_files_converted = cpb.copy_cpac_outputs_into_bids(cpac_bids_path_mapping,
                                                                        hard_link_instead_of_copy=True,
                                                                        soft_link_instead_of_copy=False)

        print("Converted {0} cpac outputs into bids format.".format(number_of_files_converted))

    return 0
