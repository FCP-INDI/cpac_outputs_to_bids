import os
import re

from scipy.io import loadmat
import numpy as np
import pandas as pd
import json
import shutil

bids_name_converter = {'compcor': 'aCompCorErr',
                       'wm': 'WhiteMatterMean',
                       'csf': 'VentriclesMean',
                       'global': 'GlobalSignalMean',
                       'pc1': 'GlobalSignalPC',
                       'motion6': ['RotZ', 'RotX', 'RotY', 'Z', 'X', 'Y'],
                       'motion24': ['RotZ', 'RotX', 'RotY', 'Z', 'X', 'Y',
                                    'RotZSq', 'RotXSq', 'RotYSq', 'ZSq', 'XSq', 'YSq',
                                    'RotZShiftBack', 'RotXShiftBack', 'RotYShiftBack',
                                    'ZShiftBack', 'XShiftBack', 'YShiftBack',
                                    'RotZShiftBackSq', 'RotXShiftBackSq', 'RotYShiftBackSq',
                                    'ZShiftBackSq', 'XShiftBackSq', 'YShiftBackSq'],
                       'linear': 'Linear',
                       'quadratic': 'Quadratic',
                       'despike': 'Outlier',
                       'constant': 'Constant',
                       'frame_wise_displacement': 'FD',
                       'max_displacement': 'MaxDisplacement'
                       }


def extract_bids_variants_from_cpac_outputs(cpac_output_bids_dictionary):
    """

    Using a list of cpac outputs determine "variants" to use in the bids conversion. CPAC pipelines have forks and
    iterables that create different versions of the outputs by varying the parameters or methods used in their
    generation. The main set of parameters, that are common across all outputs, are properties of the pipeline. The
    parameters that vary between version of the outputs will be encoded as variants. This system allows us to reduce the
    list of parameters in the variant field of the bids file path. This function identifies variants by comparing all of
    the file paths and finding ways in which they differ.

    :param cpac_output_bids_dictionary: directory to add variant information to
    :return: bids dictionary containing variant information

    """

    if not cpac_output_bids_dictionary:
        raise ValueError("Input parameter cpac_output_bids_dictionary is empty!")

    if not isinstance(cpac_output_bids_dictionary, dict):
        raise ValueError("Input parameter cpac_output_bids_dictionary should be a dictionary, not a {0}!".format(
            cpac_output_bids_dictionary.__class__.__name__))

    # format for selector
    cpac_nvr_selector_template = re.compile(
        "_compcor_ncomponents_([0-9]*)_selector_pc1([01]).linear([01]).wm([01]).global([01]).motion([01])."
        "quadratic([01]).gm([01]).compcor([01]).csf([01])")
    cpac_nvr_selector_tags = ['NCompCor',
                              'FirstPC',
                              'LinearTrend',
                              'WhiteMatter',
                              'Global',
                              'QuadraticTrend',
                              'GreyMatter',
                              'CompCor',
                              'CSF']

    # first we need to extract the components of the path that contain parameter information, this entails removing
    # extraneous information such as base path and participant identifiers
    derivative_variant_map = {}
    unique_parameter_strings = {}
    for cpac_output in cpac_output_bids_dictionary:

        cpac_output_parts = cpac_output.split('/')

        # throw away all parts of the path before the pipeline name
        while cpac_output_parts and 'pipeline' not in cpac_output_parts[0]:
            cpac_output_parts.pop(0)

        if not cpac_output_parts or 'pipeline' not in cpac_output_parts[0]:
            raise ValueError("Could not find pipeline in output {0}".format(cpac_output))

        cpac_output_pipeline = cpac_output_parts.pop(0).replace("pipeline_", "")

        # the name of a fork may be appended to the pipeline name with a "__", take the last occurrence
        if "__" in cpac_output_pipeline:
            cpac_output_pipeline = "__".join(cpac_output_pipeline.split("__")[0:-1])

        if cpac_output_pipeline not in unique_parameter_strings:
            unique_parameter_strings[cpac_output_pipeline] = set()

        # discard the next part of the path which should include the participant id and unique id
        cpac_output_parts.pop(0)

        # if the next segment includes scan information ("_scan") discard it
        if "_scan" in cpac_output_parts[0]:
            cpac_output_parts.pop(0)

        # the next should be the derivative name, lets hold on to it for now
        cpac_output_derivative = cpac_output_parts.pop(0)
        if 'path_files' in cpac_output_derivative:
            continue

        if "zstat" in cpac_output_derivative:
            if "vmhc" in cpac_output_derivative:
                cpac_output_parts.insert(0, "VarNorm")
            elif "dr" in cpac_output_derivative or "sca_tempreg" in cpac_output_derivative:
                cpac_output_parts.insert(0, "ZStat")
            else:
                print("Warning! Do not know how to interpret \"zstat\" for derivative {0}.".format(
                    cpac_output_derivative))

        if "fisher_zstd" in cpac_output_derivative:
            cpac_output_parts.insert(0, "Fisher")
        elif "zstd" in cpac_output_derivative:
            cpac_output_parts.insert(0, "ZScore")

        individual_parameter_strings = []
        for cpac_output_part in cpac_output_parts:
            for parameter_key in ['selector', 'hp', 'lp', 'bandpass_freqs', 'fwhm', 'Fisher', 'VarNorm', 'ZScore',
                                  'ZStat', '_threshold_']:
                if parameter_key in cpac_output_part:
                    individual_parameter_strings.append(cpac_output_part)

        if individual_parameter_strings:
            individual_parameter_string = '/'.join(individual_parameter_strings)

            # the way in which filtering is encoded is inconsistent, sometimes it is _bandpass_freqs and
            # other times it is an _hp/_lp pair, convert them all to the _bandpass_freqs
            hp_lp_pair = re.search("_hp_([-+]?\d*\.\d+|\d+)/_lp_([-+]?\d*\.\d+|\d+)", individual_parameter_string)
            if hp_lp_pair:
                individual_parameter_string = individual_parameter_string.replace(hp_lp_pair.group(0),
                                                                                  '_bandpass_freqs_{0}.{1}'.format(
                                                                                      hp_lp_pair.groups()[0],
                                                                                      hp_lp_pair.groups()[1]))
            else:
                lp_hp_pair = re.search("_lp_([-+]?\d*\.\d+|\d+)/_hp_([-+]?\d*\.\d+|\d+)", individual_parameter_string)
                if lp_hp_pair:
                    individual_parameter_string = individual_parameter_string.replace(lp_hp_pair.group(0),
                                                                                      '_bandpass_freqs_{0}.{1}'.format(
                                                                                          lp_hp_pair.groups()[1],
                                                                                          lp_hp_pair.groups()[0]))

            derivative_variant_map[cpac_output] = (cpac_output_pipeline, individual_parameter_string)
            unique_parameter_strings[cpac_output_pipeline].add(individual_parameter_string)

    parameter_variant_mapping = {}
    for (pipeline_name, pipeline_parameter_strings) in unique_parameter_strings.items():

        number_fwhm_dict = {}
        number_bandpass_dict = {}
        number_thresh_dict = {}
        number_selector_dict = {}
        number_standardization_dict = {}

        number_pipeline_parameter_strings = len(pipeline_parameter_strings)

        if pipeline_name not in parameter_variant_mapping:
            parameter_variant_mapping[pipeline_name] = {}

        for parameter_string in pipeline_parameter_strings:

            if parameter_string not in parameter_variant_mapping[pipeline_name]:
                parameter_variant_mapping[pipeline_name][parameter_string] = ''

            for parameter in parameter_string.split('/'):

                if 'threshold' in parameter:
                    thresh_match = re.search('_threshold_([-+]?\d*\.\d+|\d+)', parameter)
                    if thresh_match:
                        thresh = thresh_match.groups()[0]
                        if thresh not in number_thresh_dict:
                            number_thresh_dict[thresh] = 0
                        number_thresh_dict[thresh] += 1

                if 'fwhm' in parameter:
                    fwhm_match = re.search('_fwhm_([-+]?\d*\.\d+|\d+)', parameter)
                    if fwhm_match:
                        fwhm = fwhm_match.groups()[0]
                        if fwhm not in number_fwhm_dict:
                            number_fwhm_dict[fwhm] = 0
                        number_fwhm_dict[fwhm] += 1
                    else:
                        raise ValueError(
                            "FWHM found in parameter, but could not extract width,"
                            " is it formatted correctly? {0}".format(
                                parameter))

                if 'bandpass' in parameter:
                    bandpass_match = re.search('_bandpass_freqs_([-+]?\d*\.\d+|\d+)\.([-+]?\d*\.\d+|\d+)', parameter)
                    if bandpass_match:
                        bandpass_freqs = (bandpass_match.groups()[0], bandpass_match.groups()[1])
                        if bandpass_freqs not in number_bandpass_dict:
                            number_bandpass_dict[bandpass_freqs] = 0
                        number_bandpass_dict[bandpass_freqs] += 1
                    else:
                        raise ValueError(
                            "Bandpass found in parameter, but could not extract frequencies, is it formatted"
                            " correctly? {0}".format(
                                parameter))

                if 'selector' in parameter:
                    if parameter not in number_selector_dict:
                        number_selector_dict[parameter] = 0
                    number_selector_dict[parameter] += 1

                for standardization in ['Fisher', 'VarNorm', 'ZScore', 'ZStat']:
                    if standardization in parameter:
                        if standardization not in number_standardization_dict:
                            number_standardization_dict[standardization] = 0
                        number_standardization_dict[standardization] += 1

        # now lets go through and consolidate the NVR strings, a bit of a pain in the arse
        selector_items_dict = {}
        for (selector, number_selector) in number_selector_dict.items():
            if number_selector < number_pipeline_parameter_strings:
                selector_groups = cpac_nvr_selector_template.match(selector)
                if selector_groups:
                    for (tag_index, tag) in zip(range(0, len(cpac_nvr_selector_tags)), cpac_nvr_selector_tags):
                        if tag not in selector_items_dict:
                            selector_items_dict[tag] = {}
                        if selector_groups.groups()[tag_index] not in selector_items_dict[tag]:
                            selector_items_dict[tag][selector_groups.groups()[tag_index]] = []
                        selector_items_dict[tag][selector_groups.groups()[tag_index]].append(selector)
                else:
                    raise ValueError(
                        "Error! selector {0} does not match template string, is it malformed?".format(selector))

        # go through and remove items with only one value, those are invariant across derivatives.
        selector_string_items = {}
        for tag in cpac_nvr_selector_tags:
            if tag in selector_items_dict:
                selector_item = selector_items_dict[tag]
                if len(selector_item) > 1:
                    for tag_value in selector_item:
                        for selector in selector_items_dict[tag][tag_value]:
                            if selector not in selector_string_items:
                                selector_string_items[selector] = []
                            if tag_value == '1':
                                selector_string_items[selector].append(tag)
                            elif tag_value == '0':
                                selector_string_items[selector].append("No"+tag)
                            else:
                                selector_string_items[selector].append(tag+tag_value)

        for selector in selector_string_items:
            for parameter_string in parameter_variant_mapping[pipeline_name]:
                if selector in parameter_string:
                    parameter_variant_mapping[pipeline_name][parameter_string] += \
                        "".join(selector_string_items[selector])

        for bandpass in number_bandpass_dict:
            # if bandpass was not in every string, it is a variant, add it
            if number_bandpass_dict[bandpass] < number_pipeline_parameter_strings:
                for parameter_string in parameter_variant_mapping[pipeline_name]:
                    bandpass_string = "_bandpass_freqs_{0}.{1}".format(bandpass[0], bandpass[1])
                    if bandpass_string in parameter_string:
                        parameter_variant_mapping[pipeline_name][parameter_string] += \
                            "Bandpassed{0}to{1}".format(bandpass[0].split(".")[-1], bandpass[1].split(".")[-1])

        for standardization in number_standardization_dict:
            if number_standardization_dict[standardization] < number_pipeline_parameter_strings:
                for parameter_string in parameter_variant_mapping[pipeline_name]:
                    if standardization in parameter_string:
                        parameter_variant_mapping[pipeline_name][parameter_string] += \
                            standardization

        for fwhm in number_fwhm_dict:
            # if fwhm was not in every string, it is a variant, add it
            if number_fwhm_dict[fwhm] < number_pipeline_parameter_strings:
                for parameter_string in parameter_variant_mapping[pipeline_name]:
                    if "_fwhm_{0}".format(fwhm) in parameter_string:
                        parameter_variant_mapping[pipeline_name][parameter_string] += \
                            "Smoothed{0}".format(fwhm)

    # for pipeline_name in unique_parameter_strings:
    #     print("Found {0} unique pipelines for pipeline {1}:".format(len(unique_parameter_strings[pipeline_name]),
    #                                                                 pipeline_name))
    #     for parameter_string in unique_parameter_strings[pipeline_name]:
    #         print("{0} -> {1}".format(parameter_string,
    #                                   parameter_variant_mapping[cpac_output_pipeline][parameter_string]))

    for cpac_output, (cpac_output_pipeline, individual_parameter_string) in derivative_variant_map.items():
        if cpac_output in cpac_output_bids_dictionary:

            if "variant" not in cpac_output_bids_dictionary[cpac_output]:
                cpac_output_bids_dictionary[cpac_output]["variant"] = ""
            # elif cpac_output_bids_dictionary[cpac_output]["variant"]:
            #     print("Variant {0} already exists for {1}, appending".format(
            #         cpac_output_bids_dictionary[cpac_output]["variant"], cpac_output))

            cpac_output_bids_dictionary[cpac_output]["variant"] += parameter_variant_mapping[cpac_output_pipeline][
                individual_parameter_string]

    return cpac_output_bids_dictionary


def extract_bids_derivative_info_from_cpac_outputs(cpac_output_list, cpac_data_config_dict,
                                                   cpac_output_bids_dictionary=None,
                                                   derivatives_to_write=None):
    """
    Using a list of cpac outputs determine their BIDs derivative types and other elements to use in the bids conversion.

    :param cpac_output_list: a list of paths to CPAC output files
    :param cpac_data_config_dict: CPAC data configuration file used to generate output, helps map derivatives to the
        original data
    :param cpac_output_bids_dictionary: a dictionary that maps cpac outputs to their BIDs derivative elements, empty by
        default
    :param derivatives_to_write: list of cpac derivatives that are desired, others will be ignored
    :return: a dictionary that maps the output name to the BIDs derivative elements
    """

    if not cpac_output_list:
        raise ValueError("Input parameter cpac_output_list is empty!")

    if not cpac_data_config_dict:
        raise ValueError("cpac_data_config_dict should not be empty")

    if not cpac_output_bids_dictionary:
        cpac_output_bids_dictionary = {}

    if not derivatives_to_write:
        derivatives_to_write = []

    output_path_template = "{bids_root}/derivatives/pipeline-{pipeline}/sub-{sub}/ses-{ses}/{source_type}"

    # go through all of the supported derivatives and match them to the bids derivative elements
    for cpac_output in cpac_output_list:

        # reinitialize the elements dictionary
        derivative_elements = {}

        # begin with the pipeline, subject, session, and task information, things that can be easily extracted using
        # regular expressions
        pipeline_match = re.search(r"pipeline_(.*?)/", cpac_output)
        if pipeline_match:
            derivative_elements["pipeline"] = pipeline_match.groups()[0].split("__")[0]
        else:
            # print a Warning that we could not find pipeline
            print("Warning: Could not extract essential parameter [pipeline] from {0}, is it malformed?".format(
                cpac_output))
            continue

        participant_match = re.search(r"sub-(.*?)[_/]", cpac_output)
        if participant_match:
            derivative_elements["sub"] = participant_match.groups()[0]
        else:
            # print a Warning that we could not find pipeline
            print("Warning: Could not extract essential parameter [participant] from {0}, is it malformed?".format(
                cpac_output))
            continue

        session_match = re.search(r"ses-(.*?)/", cpac_output)
        if session_match:
            derivative_elements["ses"] = session_match.groups()[0]

        cpac_output_derivative = cpac_output[session_match.end():].split('/')[0]
        if not cpac_output_derivative:
            print("Warning! Could not extract derivative from {0}, skipping ...".format(cpac_output))
            continue

        if derivatives_to_write and cpac_output_derivative not in derivatives_to_write:
            continue

        task_match = re.search(r"task-(.*?)[_/]", cpac_output)
        if task_match:
            derivative_elements["task"] = task_match.groups()[0]

        run_match = re.search(r"run-(.*?)[_/]", cpac_output)
        if run_match:
            derivative_elements["run"] = run_match.groups()[0]

        # begin with transforms
        if 'xfm' in cpac_output:
            derivative_elements["derivative"] = "affine"
            derivative_elements["source_type"] = "anat"
            derivative_elements["bids_file_extension"] = "mat"

            if 'functional_to_anat_linear_xfm' in cpac_output:
                derivative_elements["source_type"] = "func"
                derivative_elements["target"] = "T1w"
                derivative_elements["variant"] = "flirtBBR"
            else:
                derivative_elements["target"] = "MNI305"
                if 'symmetric' in cpac_output:
                    derivative_elements["target"] = "MNI305Sym"

                derivative_elements["variant"] = "ANTs"
                if 'affine' in cpac_output:
                    derivative_elements["variant"] += 'Affine'
                elif 'initial' in cpac_output:
                    derivative_elements["variant"] += 'Initial'
                elif 'rigid' in cpac_output:
                    derivative_elements["variant"] += 'Rigid'
                elif 'nonlinear' in cpac_output:
                    derivative_elements["derivative"] = 'warp'
                    derivative_elements["bids_file_extension"] = ""

                if 'mni_to_anatomical' in cpac_output:
                    derivative_elements["source"] = derivative_elements["target"]
                    derivative_elements["target"] = ''
                    derivative_elements["invert_source_target"] = True

        elif 'mni_normalized_anatomical' in cpac_output:
            derivative_elements["source_type"] = "anat"
            derivative_elements["variant"] = "ANTs"
            derivative_elements["derivative"] = "brain"

            derivative_elements["space"] = "MNI305"
            if "symmetric_mni_normalized_anatomical" in cpac_output:
                derivative_elements["space"] = "MNI305Sym"

        elif 'anatomical_' in cpac_output:
            derivative_elements["source_type"] = "anat"
            derivative_elements["variant"] = "ANTs"
            derivative_elements["space"] = "orig"

            if 'anatomical_brain' in cpac_output:
                derivative_elements["derivative"] = "brain"
            elif 'anatomical_gm_mask' in cpac_output:
                derivative_elements["derivative"] = "roi"
                derivative_elements["label"] = "GreyMatter"
            elif 'anatomical_wm_mask' in cpac_output:
                derivative_elements["derivative"] = "roi"
                derivative_elements["label"] = "WhiteMatter"
            elif 'anatomical_csf_mask' in cpac_output:
                derivative_elements["derivative"] = "roi"
                derivative_elements["label"] = "CSF"
            else:
                print("Warning! Do not know how to parse {0}, skipping.".format(cpac_output))

        # handle functional derivatives
        elif '_scan_' in cpac_output:

            derivative_elements["source_type"] = "func"

            if 'to_standard' in cpac_output:
                derivative_elements["space"] = "MNI305"
            elif 'in_anat' in cpac_output:
                derivative_elements["space"] = "T1w"
            else:
                derivative_elements["space"] = "orig"

            if '_network_centrality' in cpac_output or 'centrality_outputs' in cpac_output:

                # network centrality is calculated in template space
                derivative_elements["space"] = "MNI305"

                if 'eigenvector_centrality_weighted' in cpac_output:
                    derivative_elements["derivative"] = "ecw"
                elif 'eigenvector_centrality_binarize' in cpac_output:
                    derivative_elements["derivative"] = "ecb"
                elif 'degree_centrality_weighted' in cpac_output:
                    derivative_elements["derivative"] = "dcw"
                elif 'degree_centrality_binarize' in cpac_output:
                    derivative_elements["derivative"] = "dcb"
                elif 'lfcd_weighted' in cpac_output:
                    derivative_elements["derivative"] = "lfcdw"
                elif 'lfcd_binarize' in cpac_output:
                    derivative_elements["derivative"] = "lfcdb"
                else:
                    print("Warning! Do not know how to parse {0}, skipping.".format(cpac_output))
                    continue
            elif 'falff' in cpac_output:
                derivative_elements['derivative'] = 'falff'
            elif 'alff' in cpac_output:
                derivative_elements['derivative'] = 'alff'
            elif 'reho' in cpac_output:
                derivative_elements['derivative'] = 'reho'
            elif 'vmhc' in cpac_output:
                derivative_elements['derivative'] = 'vmhc'
                derivative_elements["space"] = "MNI305Sym"
            elif 'mean_functional' in cpac_output:
                derivative_elements['derivative'] = 'mean'
            elif 'functional_brain_mask' in cpac_output or 'functional_preprocessed_mask' in cpac_output:
                derivative_elements['derivative'] = 'mask'
            elif 'functional_nuisance_regressors' in cpac_output:
                derivative_elements['derivative'] = 'nuisance'
                derivative_elements["bids_file_extension"] = 'tsv'
            elif 'functional_nuisance_residuals' in cpac_output:
                derivative_elements['derivative'] = 'preprocessed'
            elif 'functional_' in cpac_output:
                derivative_elements['derivative'] = 'preprocessed'
            elif 'motion_correct' in cpac_output:
                derivative_elements['derivative'] = 'preprocessed'
                derivative_elements["variant"] = "MinimalMoCo"
            elif 'slice_time_corrected' in cpac_output:
                derivative_elements['derivative'] = 'preprocessed'
                derivative_elements["variant"] = "SliceTimeCo"
            elif 'preprocessed' in cpac_output:
                derivative_elements['derivative'] = 'preprocessed'
                derivative_elements["variant"] = "Minimal"
            elif 'sca_roi' in cpac_output:
                derivative_elements['derivative'] = 'sca'
                sca_search = re.search('_mask_(.*?)/', cpac_output)
                if sca_search:
                    derivative_elements['atlas'] = sca_search.groups()[0].replace('-', '').replace('_', '')
                else:
                    print("Warning!: could extract atlas from {0}, is it correctly formatted?".format(cpac_output))
                sca_search = re.search('sca_ROI_([0-9]+)', cpac_output)
                if sca_search:
                    derivative_elements['roi'] = sca_search.groups()[0]
                else:
                    print("Warning!: could extract roi from {0}, is it correctly formatted?".format(cpac_output))
            elif 'sca_tempreg' in cpac_output:
                derivative_elements['derivative'] = 'sca'
                derivative_elements['variant'] = 'MultiReg'
                sca_search = re.search('_mask_(.*?)/', cpac_output)
                if sca_search:
                    derivative_elements['atlas'] = sca_search.groups()[0].replace('-', '').replace('_', '')
                else:
                    print("Warning!: could extract atlas from {0}, is it correctly formatted?".format(cpac_output))
                sca_search = re.search('_maps_roi_([0-9]+)', cpac_output)
                if sca_search:
                    derivative_elements['roi'] = sca_search.groups()[0]
                else:
                    print("Warning!: could extract roi from {0}, is it correctly formatted?".format(cpac_output))
            elif 'dr_tempreg_maps' in cpac_output:
                derivative_elements['derivative'] = 'drmap'
                dr_search = re.search('_spatial_map_(.*?)/', cpac_output)
                if dr_search:
                    derivative_elements['atlas'] = dr_search.groups()[0].replace('-', '').replace('_', '')
                else:
                    print("Warning!: could extract atlas from {0}, is it correctly formatted?".format(cpac_output))
                dr_search = re.search('temp_reg_map_z*_*([0-9]{4})', cpac_output)
                if dr_search:
                    derivative_elements['roi'] = dr_search.groups()[0]
                else:
                    print("Warning!: could extract roi from {0}, is it correctly formatted?".format(cpac_output))
            elif 'spatial_map_timeseries' in cpac_output:
                derivative_elements['derivative'] = 'roisdata'
                derivative_elements['bids_file_extension'] = 'tsv'
                derivative_elements['variant'] = 'SpatReg'
                derivative_elements['space'] = ""
                dr_search = re.search('_spatial_map_(.*?)/', cpac_output)
                if dr_search:
                    derivative_elements['atlas'] = dr_search.groups()[0].replace('-', '').replace('_', '')
                else:
                    print("Warning!: could extract atlas from {0}, is it correctly formatted?".format(cpac_output))
            elif 'roi_timeseries' in cpac_output:
                derivative_elements['derivative'] = 'roisdata'
                derivative_elements['variant'] = 'Mean'
                derivative_elements['bids_file_extension'] = 'tsv'
                derivative_elements['space'] = ""
                dr_search = re.search('_mask_(.*?)/', cpac_output)
                if dr_search:
                    derivative_elements['atlas'] = dr_search.groups()[0].replace('-', '').replace('_', '')
                else:
                    print("Warning!: could extract atlas from {0}, is it correctly formatted?".format(cpac_output))

            elif 'power_params' in cpac_output or 'motion_parameters' in cpac_output:
                derivative_elements['derivative'] = 'motionstats'
                derivative_elements['bids_file_extension'] = 'tsv'
            elif 'movement_parameters' in cpac_output:
                derivative_elements['derivative'] = 'motionregressors'
                derivative_elements['bids_file_extension'] = 'tsv'
            elif 'frame_wise_displacement' in cpac_output:
                derivative_elements['derivative'] = 'motionregressors'
                derivative_elements['bids_file_extension'] = 'tsv'
            elif 'max_displacement' in cpac_output:
                derivative_elements['derivative'] = 'motionregressors'
                derivative_elements['bids_file_extension'] = 'tsv'
            elif 'functional_nuisance_regressors' in cpac_output:
                derivative_elements['derivative'] = 'nuisanceregressors'
                derivative_elements['bids_file_extension'] = 'tsv'
            else:
                print("Warning!: do not know how to parse {0}".format(cpac_output))
                continue
        else:
            continue

        bids_input_path_key = '_'.join(
            ['-'.join([key, derivative_elements[key]]) for key in ['sub', 'ses'] if key in derivative_elements]).lower()
        if bids_input_path_key == '' or bids_input_path_key == ' ':
            print("Warning! Could not construct data config key for path {0} : {1}".format(cpac_output,
                                                                                           derivative_elements))

        if bids_input_path_key in cpac_data_config_dict and cpac_data_config_dict[bids_input_path_key]:
            if 'source_type' in derivative_elements:
                if 'anat' in derivative_elements['source_type']:
                    bids_input_path = cpac_data_config_dict[bids_input_path_key][derivative_elements['source_type']]
                elif 'func' in derivative_elements['source_type']:
                    task_key = '_'.join(
                        ['-'.join([key, derivative_elements[key]]) for key in ['task', 'run'] if
                         key in derivative_elements])
                    if task_key == '' or task_key == ' ':
                        print("Warning! Could not construct data config task key for path {0} : {1}".format(
                            cpac_output,
                            derivative_elements))
                    if task_key in cpac_data_config_dict[bids_input_path_key]['rest']:
                        bids_input_path = cpac_data_config_dict[bids_input_path_key]['rest'][task_key]
                    else:
                        print("Warning! Could not find .{0}. in data config {1}.".format(task_key,
                                                                                         cpac_data_config_dict[
                                                                                             bids_input_path_key][
                                                                                             'rest']))
                        print(cpac_data_config_dict)
                        continue
                else:
                    print('Warning! Do not know how to handle source type {0}.'.format(
                        derivative_elements['source_type']))
                    continue
            else:
                print('Warning! bids dictionary does not contain source_type, it must be malformed.')
                continue
        else:
            print("Warning! Could not find {0} in data config".format(bids_input_path_key))
            continue

        bids_input_path_sub_search = re.search('/sub-.*?/', bids_input_path)
        if bids_input_path_sub_search:
            derivative_elements['bids_root'] = bids_input_path[0:bids_input_path_sub_search.start()]
        else:
            print("Warning! Could not find subject directory in output path {0}. Is it malformed?")
            continue

        derivative_elements["bids_derivative_path"] = output_path_template.format(**derivative_elements)

        if "bids_file_extension" not in derivative_elements or not derivative_elements["bids_file_extension"]:
            derivative_elements["bids_file_extension"] = ".".join(os.path.basename(cpac_output).split(".")[1:])

        bids_source_name = os.path.basename(bids_input_path).split(".")[0]
        if "invert_source_target" in derivative_elements and derivative_elements["invert_source_target"] is True:
            derivative_elements["target"] = bids_source_name
        else:
            derivative_elements["source"] = bids_source_name

        if cpac_output not in cpac_output_bids_dictionary:
            cpac_output_bids_dictionary[cpac_output] = {}
        cpac_output_bids_dictionary[cpac_output].update(derivative_elements)

    cpac_output_bids_dictionary = extract_bids_variants_from_cpac_outputs(cpac_output_bids_dictionary)

    return cpac_output_bids_dictionary


def create_bids_outputs_dictionary(cpac_output_bids_dictionary):
    """
    Creates a dictionary that maps cpac paths to bids paths

    :param cpac_output_bids_dictionary: mapping from original file name to bids dictionary, which will be converted into
        a filename
    :return: path mapping dictionary
    """

    cpac_bids_path_map = {}
    for (cpac_output_path, bids_dict) in cpac_output_bids_dictionary.items():
        if 'derivative' in bids_dict:
            bids_path = bids_dict["bids_derivative_path"] + '/' + '_'.join(
                [bids_dict['source']] +
                ['-'.join([key, bids_dict[key]]) for key in
                 ['target', 'space', 'res', 'variant', 'atlas', 'roi', 'label'] if
                 key in bids_dict and bids_dict[key]] +
                [bids_dict['derivative']]) + "." + bids_dict['bids_file_extension']

            cpac_bids_path_map[cpac_output_path] = bids_path

        else:
            print('{0} is missing {1}'.format(bids_dict, 'derivative'))

    return cpac_bids_path_map


def convert_nuisance_regressors_to_bids(cpac_regressor_filepath, bids_regressor_filepath):
    """

    :param cpac_regressor_filepath: path to cpac_regressor file to be converted
    :param bids_regressor_filepath: path to bids regressor file to be written
    :return: bids_regressor_filepath
    """

    if not cpac_regressor_filepath:
        raise ValueError("Error! cpac_regressor_filepath should not be empty!")
    if not bids_regressor_filepath:
        raise ValueError("Error! bids_regressor_filepath should not be empty!")
    if not cpac_regressor_filepath.endswith('mat'):
        raise ValueError("Error! Expecting the regressor file to be a .mat file.")

    # print("convert {0} to {1}".format(cpac_regressor_filepath, bids_regressor_filepath))

    cpac_regressor_dict = loadmat(cpac_regressor_filepath)

    regressor_matrix = []
    regressor_names = []
    for (regressor_name, regressor_data) in cpac_regressor_dict.items():
        if "__" not in regressor_name:

            if 'compcor' in regressor_name:
                for index in range(0, regressor_data.shape[1]):
                    regressor_names.append("aCompCor{0}".format(index))
            elif 'despike' in regressor_name:
                for index in range(0, regressor_data.shape[1]):
                    regressor_names.append("Outlier{0}".format(index))
            else:
                if regressor_data.shape[1] > 1:
                    if 'motion' not in regressor_name:
                        raise ValueError("Only compcor and motion regressors should have more than one component.")

                    if regressor_data.shape[1] == 6:
                        regressor_name = 'motion6'
                    elif regressor_data.shape[1] == 24:
                        regressor_name = 'motion24'

                    regressor_names += bids_name_converter[regressor_name]

                    # convert degrees into radians
                    for index in range(0, regressor_data.shape[1]):
                        if index % 6 < 3:
                            regressor_data[:, index] *= 2.0 * np.pi / 360.0
                        if index % 12 in [6, 7, 8]:
                            regressor_data[:, index] *= 2.0 * np.pi / 360.0

                else:
                    regressor_names.append(bids_name_converter[regressor_name])

            if not isinstance(regressor_matrix, np.ndarray):
                regressor_matrix = regressor_data.reshape(regressor_data.shape[0], -1)
            else:
                regressor_matrix = np.hstack((regressor_matrix, regressor_data.reshape(regressor_data.shape[0], -1)))

    regressor_df = pd.DataFrame(data=regressor_matrix, columns=regressor_names)
    regressor_df.to_csv(bids_regressor_filepath, sep='\t')

    return bids_regressor_filepath


def aggregate_movement_parameters(cpac_movement_parameter_file, bids_movement_parameter_file):
    """

    :param cpac_movement_parameter_file: path to movement parameter file created from CPAC
    :param bids_movement_parameter_file: path to the name the BIDS output file should have
    :return: the name of the output file
    """
    if not cpac_movement_parameter_file:
        raise ValueError("Error: cpac_movement_parameter_file should be not be empty.")

    if not bids_movement_parameter_file:
        raise ValueError("Error: bids_movement_parameter_file should be not be empty.")

    bids_column_names = []
    if 'fristons_twenty_four' in cpac_movement_parameter_file:
        bids_column_names = bids_name_converter['motion24']
    elif 'movement_parameters' in cpac_movement_parameter_file:
        bids_column_names = bids_name_converter['motion6']
    elif 'frame_wise_displacement' in cpac_movement_parameter_file:
        bids_column_names = [bids_name_converter['frame_wise_displacement']]
    elif 'max_displacement' in cpac_movement_parameter_file:
        bids_column_names = [bids_name_converter['max_displacement']]

    cpac_movement_parameters = pd.read_csv(cpac_movement_parameter_file, names=bids_column_names,
                                           delim_whitespace=True, comment='#')

    # if a movement file already exists, we will append the new data to it
    if os.path.exists(bids_movement_parameter_file):
        bids_movement_parameters = pd.read_csv(bids_movement_parameter_file, delimiter='\t', comment='#')

        # only include the components that haven't been included before
        columns_to_append = []
        for column_name in cpac_movement_parameters.columns:
            if column_name not in bids_movement_parameters.columns:
                columns_to_append.append(column_name)

        if columns_to_append:
            bids_movement_parameters = pd.concat(
                (bids_movement_parameters, cpac_movement_parameters[columns_to_append]), axis=1)
    else:
        bids_movement_parameters = cpac_movement_parameters

    bids_movement_parameters.to_csv(bids_movement_parameter_file, sep='\t', index=False)

    return bids_movement_parameter_file


def aggregate_movement_statistics(cpac_movement_statistics_file, bids_movement_statistics_file):
    """

    :param cpac_movement_statistics_file: path to cpac_movement_statistics file 
    :param bids_movement_statistics_file: path of bids_movement_statistics_file to be written
    :return: output file path
    """

    if not cpac_movement_statistics_file:
        return ValueError("cpac_movement_statistics_file must be the path to an existing file, not empty.")

    with open(cpac_movement_statistics_file, 'r') as cpac_stats_fd:
        cpac_stats_lines = cpac_stats_fd.readlines()

    if len(cpac_stats_lines) != 2:
        raise ValueError(
            "The number of rows in {0} ({1}) is not 2, is it correctly formatted?".format(
                cpac_movement_statistics_file, len(cpac_stats_lines)))

    header_line_values = cpac_stats_lines[0].rstrip().split(',')
    data_line_values = cpac_stats_lines[1].rstrip().split(',')

    if len(header_line_values) != len(data_line_values):
        raise ValueError(
            "Number of values in header line does not match number of values in data row, "
            "is motion statistics file {0} formatted correctly?".format(
                cpac_movement_statistics_file))

    cpac_movement_statistics = {statistics_key.replace(' ', ''): statistics_value.replace(' ', '') for
                                (statistics_key, statistics_value) in
                                zip(header_line_values, data_line_values)}

    if os.path.isfile(bids_movement_statistics_file):
        with open(bids_movement_statistics_file, 'r') as in_stream:
            bids_movement_statistics = json.load(in_stream)
    else:
        bids_movement_statistics = {}

    # clean up the values to make them bids compliant
    if 'Subject' in cpac_movement_statistics:
        new_values = {stat_values.split('-')[0]: stat_values.split('-')[1] for stat_values in
                      cpac_movement_statistics['Subject'].split('_')}

        if 'sub' in new_values:
            if 'sub' in bids_movement_statistics and bids_movement_statistics['sub'] != new_values['sub']:
                raise ValueError("Bids movement statistics file {0} already exists, but appears to contain "
                                 "entries for different dataset (sub {1} != {2})".format(bids_movement_statistics_file,
                                                                                         bids_movement_statistics[
                                                                                             'sub'],
                                                                                         new_values['sub']))
            else:
                bids_movement_statistics.update(new_values)
        else:
            raise ValueError("Could not find 'sub' bids key in CPAC movement statistics file {0}".format(
                cpac_movement_statistics_file))

        if 'ses' in new_values:
            if 'ses' in bids_movement_statistics and bids_movement_statistics['ses'] != new_values['ses']:
                raise ValueError("Bids movement statistics file {0} already exists, but appears to contain "
                                 "entries for different dataset (ses {1} != {2})".format(bids_movement_statistics_file,
                                                                                         bids_movement_statistics[
                                                                                             'ses'],
                                                                                         new_values['ses']))
            else:
                bids_movement_statistics.update(new_values)
        else:
            raise ValueError("Could not find 'ses' bids key in CPAC movement statistics file {0}".format(
                cpac_movement_statistics_file))

    if 'Scan' in cpac_movement_statistics:

        new_values = {stat_values.split('-')[0]: stat_values.split('-')[1] for stat_values in
                      cpac_movement_statistics['Scan'].split('_')}

        if 'task' in new_values:
            if 'task' in bids_movement_statistics and bids_movement_statistics['task'] != new_values['task']:
                raise ValueError("Bids movement statistics file {0} already exists, but appears to contain "
                                 "entries for different dataset (task {1} != {2})".format(bids_movement_statistics_file,
                                                                                          bids_movement_statistics[
                                                                                              'task'],
                                                                                          new_values['task']))
            else:
                bids_movement_statistics['task'] = new_values['task']
        else:
            raise ValueError("Could not find 'task' bids key in CPAC movement statistics file {0}".format(
                cpac_movement_statistics_file))

        if 'run' in new_values:
            if 'run' in bids_movement_statistics and bids_movement_statistics['run'] != new_values['run']:
                raise ValueError("Bids movement statistics file {0} already exists, but appears to contain "
                                 "entries for different dataset (run {1} != {2})".format(bids_movement_statistics_file,
                                                                                         bids_movement_statistics[
                                                                                             'run'],
                                                                                         new_values['run']))
            else:
                bids_movement_statistics['run'] = new_values['run']

    # now copy over keys, cleaning each one so that it is BIDs compliant
    for statistics_key, statistics_value in cpac_movement_statistics.items():
        if statistics_key in ['Subject', 'Scan']:
            continue
        bids_movement_statistics[statistics_key.replace('-', '').replace('_', '').replace(' ', '')] = float(
            statistics_value)

    with open(bids_movement_statistics_file, 'w') as out_stream:
        json.dump(bids_movement_statistics, out_stream, indent=4)

    return bids_movement_statistics_file


def convert_cpac_roi_tc_to_bids(cpac_roi_file, bids_roi_file):
    """

    :param cpac_roi_file: roi file output by CPAC
    :param bids_roi_file: roi file to write in bids format
    :return: file name of output file
    """

    if not cpac_roi_file:
        raise ValueError("cpac_roi_file is empty, expected path to existing file.")

    if not bids_roi_file:
        raise ValueError("bids_roi_file is empty, expected path to file to be created.")

    timeseries_df = pd.read_csv(cpac_roi_file, delim_whitespace=True, header=None, comment="#")

    if "spatial_map_timeseries" in cpac_roi_file:
        timeseries_df.columns = ["ROI{0}SpatReg".format(roi_num+1) for roi_num in range(0, timeseries_df.shape[1])]

    elif "roi_stats" in cpac_roi_file:

        with open(cpac_roi_file, 'r') as infile_descriptor:
            header_lines = infile_descriptor.readlines()

        header_line = ''
        for line in header_lines:
            if "Mean" in line:
                header_line = line
                break

        column_names = [column_name.replace("Mean_", "ROI") + "Mean" for column_name in
                        header_line.lstrip("#").rstrip().split(None) if "#" not in column_name]

        timeseries_df.columns = column_names

    else:
        raise ValueError("Do not know how to process {0}".format(cpac_roi_file))

    if isinstance(timeseries_df, pd.DataFrame):
        timeseries_df.to_csv(bids_roi_file, sep='\t', index=False)

    return bids_roi_file


def copy_cpac_outputs_into_bids(cpac_bids_path_map, hard_link_instead_of_copy=False, soft_link_instead_of_copy=False):
    """
    Copy (or link) cpac outputs into the bids structure

    :param cpac_bids_path_map: dictionary that maps cpac output paths to bids paths
    :param hard_link_instead_of_copy: save space by creating links to the files instead of copying, hard links are used
            to retain the data if you delete the CPAC output. Hard links will not work when linking across filesystems
            This will not effect files that are produced, such as TSV
            and JSON files, as a result of the bidsification process, they cannot be created through links.
            default = False.
    :param soft_link_instead_of_copy: save space by creating links to the files instead of copying, symbolic (soft)
            links are used, which unlike hard links can span different filesystems. E.g. you could create links from
            your local drive to files that exist on a network attached storage server.

            to retain the data if you delete the CPAC output. This will not effect files that are produced, such as TSV
            and JSON files, as a result of the bidsification process, they cannot be created through links.
            default = False.

    :return: number of files transferred.
    """

    if not cpac_bids_path_map:
        raise ValueError("Expected cpac_bids_path_map to contain a dictionary, got nothing.")

    if not isinstance(cpac_bids_path_map, dict):
        raise ValueError(
            "Expected cpac_bids_path_map to contain a dictionary, got a {0}".format(type(cpac_bids_path_map)))

    number_of_files_converted = 0

    for (cpac_output_path, bids_output_path) in cpac_bids_path_map.items():

        # first determine whether a converter is required, determine this from the bids path instead of the cpac path
        # since it is a bit cleaner, and multiple cpac files map to the same bids file, so the logic is easier
        if 'nuisance.tsv' in bids_output_path:
            convert_nuisance_regressors_to_bids(cpac_output_path, bids_output_path)

        elif 'motionregressors.tsv' in bids_output_path:
            aggregate_movement_parameters(cpac_output_path, bids_output_path)

        elif 'motionstats.json' in bids_output_path:
            aggregate_movement_statistics(cpac_output_path, bids_output_path)

        elif 'roisdata.tsv' in bids_output_path:
            convert_cpac_roi_tc_to_bids(cpac_output_path, bids_output_path)

        elif hard_link_instead_of_copy is True:

            try:
                os.link(cpac_output_path, bids_output_path)
            except Exception as e:
                print("Error: Could not create hard link. Are {0} and {1} on the same device? {2}".format(
                    cpac_output_path, bids_output_path, e.message))
                raise e

        elif soft_link_instead_of_copy is True:
            os.symlink(cpac_output_path, bids_output_path)

        else:
            shutil.copy(cpac_output_path, bids_output_path)

        number_of_files_converted += 1

    return number_of_files_converted
