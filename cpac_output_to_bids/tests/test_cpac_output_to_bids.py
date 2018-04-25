from unittest import TestCase

# These are mostly smoke tests designed to find major flaws in the code, need some more sensitive
# tests that better validate the outputs


class TestCPACOutputConversion(TestCase):

    def test_extract_path_source_extension(self):

        import cpac_output_to_bids as cpb
        import json
        import yaml
        import os

        cpac_paths_file_path = os.path.dirname(__file__) + '/test_files/sub-M10933594_ses-NFB3.json'
        cpac_derivatives_dictionary_path = os.path.dirname(
            __file__) + '/test_files/test_cpac_derivatives_dictionary.json'
        cpac_data_config_file = os.path.dirname(__file__) + '/test_files/pre_skullstrip_data_config.yml'
        output_reference_dictionary = os.path.dirname(__file__) + '/test_files/sub-M10933594_ses-NFB3_cpac_to_bids.json'

        with open(cpac_derivatives_dictionary_path, 'r') as json_input_stream:
            cpac_derivatives_dictionary = json.load(json_input_stream)

        cpac_write_derivative_list = [key for (key, value) in cpac_derivatives_dictionary.items() if
                                      value['write_derivative'] is True]

        with open(cpac_paths_file_path, 'r') as json_input_stream:
            cpac_output_file_paths = json.load(json_input_stream)

        with open(cpac_data_config_file, 'r') as yaml_input_stream:
            cpac_data_config_list = yaml.load(yaml_input_stream)

        cpac_data_config_dict = {}
        for cpac_data_config in cpac_data_config_list:
            cpac_data_config_dict[
                "_".join([cpac_data_config['subject_id'], cpac_data_config['unique_id']]).lower()] = cpac_data_config

        bids_dictionary = cpb.extract_bids_derivative_info_from_cpac_outputs(cpac_output_file_paths,
                                                                             cpac_data_config_dict, {},
                                                                             cpac_write_derivative_list)

        assert isinstance(bids_dictionary, dict)
        assert bids_dictionary != {}

        cpac_bids_path_mapping = cpb.create_bids_outputs_dictionary(bids_dictionary)
        assert isinstance(cpac_bids_path_mapping, dict)

        # with open(os.path.dirname(__file__)+'/test_files/test_paths.json', 'w') as json_output_stream:
        #     json.dump(cpac_bids_path_mapping, json_output_stream, indent=4)

        # compare against hand verified reference mapping
        with open(output_reference_dictionary, 'r') as json_input_stream:
            reference_cpac_bids_path_mapping = json.load(json_input_stream)

        for cpac_path, bids_path in cpac_bids_path_mapping.items():
            assert bids_path == reference_cpac_bids_path_mapping[cpac_path]

        # eventually could add a test for bids compliance

    def test_convert_nuisance_regressors_to_bids(self):
        import os
        import cpac_output_to_bids as cpb

        cpac_regressor_file_path = os.path.dirname(__file__) + \
            "/test_files/pipeline_cpac/sub-M10933594_ses-NFB3/functional_nuisance_regressors/_scan_task-morald/" \
            "_compcor_ncomponents_5_selector_pc10.linear1.wm0.global1.motion1.quadratic1.gm0.compcor1.csf1/" \
            "nuisance_regressors.mat"

        bids_regressor_file_path = "/tmp/bids_nuisance_regressors.tsv"

        outfile = cpb.convert_nuisance_regressors_to_bids(cpac_regressor_file_path, bids_regressor_file_path)

        assert os.path.isfile(outfile)

    def test_aggregate_movement_parameters(self):
        import cpac_output_to_bids as cpb
        import os
        import yaml
        import json

        cpac_output_path = os.path.dirname(__file__) + '/test_files/pipeline_cpac/'
        cpac_data_config_file = os.path.dirname(__file__) + '/test_files/pre_skullstrip_data_config.yml'
        cpac_derivatives_dictionary_path = os.path.dirname(
            __file__) + '/test_files/test_cpac_derivatives_dictionary.json'

        with open(cpac_derivatives_dictionary_path, 'r') as json_input_stream:
            cpac_derivatives_dictionary = json.load(json_input_stream)

        cpac_write_derivative_list = [key for (key, value) in cpac_derivatives_dictionary.items() if
                                      value['write_derivative'] is True]

        with open(cpac_data_config_file, 'r') as yaml_input_stream:
            cpac_data_config_list = yaml.load(yaml_input_stream)

        cpac_data_config_dict = {}
        for cpac_data_config in cpac_data_config_list:
            cpac_data_config_dict[
                "_".join([cpac_data_config['subject_id'], cpac_data_config['unique_id']]).lower()] = cpac_data_config

        movement_parameter_files = []
        for path_root, directories, file_names in os.walk(cpac_output_path):
            for file_name in file_names:
                full_file_path = os.path.join(path_root, file_name)
                if file_name in ["fristons_twenty_four.1D", "FD.1D", "max_displacement.1D"]:
                    movement_parameter_files.append(full_file_path)

        bids_dictionary = cpb.extract_bids_derivative_info_from_cpac_outputs(movement_parameter_files,
                                                                             cpac_data_config_dict, {},
                                                                             cpac_write_derivative_list)

        cpac_bids_path_mapping = cpb.create_bids_outputs_dictionary(bids_dictionary)

        for (cpac_output_path, bids_output_path) in cpac_bids_path_mapping.items():

            if not os.path.isdir(os.path.dirname(bids_output_path)):
                os.makedirs(os.path.dirname(bids_output_path), exist_ok=True)

            outfile = cpb.aggregate_movement_parameters(cpac_output_path, bids_output_path)
            assert os.path.isfile(outfile)

    def test_aggregate_movement_statistics(self):
        import cpac_output_to_bids as cpb
        import os
        import yaml
        import json

        cpac_output_path = os.path.dirname(__file__) + '/test_files/pipeline_cpac/'
        cpac_data_config_file = os.path.dirname(__file__) + '/test_files/pre_skullstrip_data_config.yml'
        cpac_derivatives_dictionary_path = os.path.dirname(
            __file__) + '/test_files/test_cpac_derivatives_dictionary.json'

        with open(cpac_derivatives_dictionary_path, 'r') as json_input_stream:
            cpac_derivatives_dictionary = json.load(json_input_stream)

        cpac_write_derivative_list = [key for (key, value) in cpac_derivatives_dictionary.items() if
                                      value['write_derivative'] is True]

        with open(cpac_data_config_file, 'r') as yaml_input_stream:
            cpac_data_config_list = yaml.load(yaml_input_stream)

        cpac_data_config_dict = {}
        for cpac_data_config in cpac_data_config_list:
            cpac_data_config_dict[
                "_".join(
                    [cpac_data_config['subject_id'], cpac_data_config['unique_id']]).lower()] = cpac_data_config

        movement_parameter_files = []
        for path_root, directories, file_names in os.walk(cpac_output_path):
            for file_name in file_names:
                full_file_path = os.path.join(path_root, file_name)
                if file_name in ['pow_params.txt', 'motion_parameters.txt']:
                    movement_parameter_files.append(full_file_path)

        bids_dictionary = cpb.extract_bids_derivative_info_from_cpac_outputs(movement_parameter_files,
                                                                             cpac_data_config_dict, {},
                                                                             cpac_write_derivative_list)

        cpac_bids_path_mapping = cpb.create_bids_outputs_dictionary(bids_dictionary)

        for (cpac_output_path, bids_output_path) in cpac_bids_path_mapping.items():

            if not os.path.isdir(os.path.dirname(bids_output_path)):
                os.makedirs(os.path.dirname(bids_output_path), exist_ok=True)

            outfile = cpb.aggregate_movement_statistics(cpac_output_path, bids_output_path)
            assert os.path.isfile(outfile)

    def test_convert_cpac_roi_tc_to_bids(self):
        import cpac_output_to_bids as cpb
        import os
        import yaml
        import json

        cpac_output_path = os.path.dirname(__file__) + '/test_files/pipeline_cpac/'
        cpac_data_config_file = os.path.dirname(__file__) + '/test_files/pre_skullstrip_data_config.yml'
        cpac_derivatives_dictionary_path = os.path.dirname(
            __file__) + '/test_files/test_cpac_derivatives_dictionary.json'

        with open(cpac_derivatives_dictionary_path, 'r') as json_input_stream:
            cpac_derivatives_dictionary = json.load(json_input_stream)

        cpac_write_derivative_list = [key for (key, value) in cpac_derivatives_dictionary.items() if
                                      value['write_derivative'] is True]

        with open(cpac_data_config_file, 'r') as yaml_input_stream:
            cpac_data_config_list = yaml.load(yaml_input_stream)

        cpac_data_config_dict = {}
        for cpac_data_config in cpac_data_config_list:
            cpac_data_config_dict[
                "_".join([cpac_data_config['subject_id'], cpac_data_config['unique_id']]).lower()] = cpac_data_config

        roi_timeseries_files = []
        for path_root, directories, file_names in os.walk(cpac_output_path):
            for file_name in file_names:
                full_file_path = os.path.join(path_root, file_name)
                if 'roi_timeseries' in full_file_path:
                    roi_timeseries_files.append(full_file_path)

        bids_dictionary = cpb.extract_bids_derivative_info_from_cpac_outputs(roi_timeseries_files,
                                                                             cpac_data_config_dict, {},
                                                                             cpac_write_derivative_list)

        cpac_bids_path_mapping = cpb.create_bids_outputs_dictionary(bids_dictionary)

        for (cpac_output_path, bids_output_path) in cpac_bids_path_mapping.items():

            if not os.path.isdir(os.path.dirname(bids_output_path)):
                os.makedirs(os.path.dirname(bids_output_path), exist_ok=True)

            outfile = cpb.convert_cpac_roi_tc_to_bids(cpac_output_path, bids_output_path)
            assert os.path.isfile(outfile)

    def test_everything(self):

        import cpac_output_to_bids as cpb
        import json
        import yaml
        import os

        cpac_paths_file_path = os.path.dirname(__file__) + '/test_files/sub-M10933594_ses-NFB3.json'
        cpac_derivatives_dictionary_path = os.path.dirname(
            __file__) + '/test_files/test_cpac_derivatives_dictionary.json'
        cpac_data_config_file = os.path.dirname(__file__) + '/test_files/pre_skullstrip_data_config.yml'
        output_reference_dictionary = os.path.dirname(__file__) + '/test_files/sub-M10933594_ses-NFB3_cpac_to_bids.json'

        # first thing, create a series of empty cpac output files using "touch", these will mimic nifti files
        # which are either copied or linked to in the process of converting from CPAC to BIDS. Other files that
        # are actually converted in the process can already be found in the test_files directory

        with open(cpac_derivatives_dictionary_path, 'r') as json_input_stream:
            cpac_derivatives_dictionary = json.load(json_input_stream)

        cpac_write_derivative_list = [key for (key, value) in cpac_derivatives_dictionary.items() if
                                      value['write_derivative'] is True]

        with open(cpac_paths_file_path, 'r') as json_input_stream:
            cpac_output_file_paths = json.load(json_input_stream)

        with open(cpac_data_config_file, 'r') as yaml_input_stream:
            cpac_data_config_list = yaml.load(yaml_input_stream)

        cpac_data_config_dict = {}
        for cpac_data_config in cpac_data_config_list:
            cpac_data_config_dict[
                "_".join([cpac_data_config['subject_id'], cpac_data_config['unique_id']]).lower()] = cpac_data_config

        bids_dictionary = cpb.extract_bids_derivative_info_from_cpac_outputs(cpac_output_file_paths,
                                                                             cpac_data_config_dict, {},
                                                                             cpac_write_derivative_list)

        assert isinstance(bids_dictionary, dict)
        assert bids_dictionary != {}

        cpac_bids_path_mapping = cpb.create_bids_outputs_dictionary(bids_dictionary)
        assert isinstance(cpac_bids_path_mapping, dict)

        # with open(os.path.dirname(__file__)+'/test_files/test_paths.json', 'w') as json_output_stream:
        #     json.dump(cpac_bids_path_mapping, json_output_stream, indent=4)

        # compare against hand verified reference mapping
        with open(output_reference_dictionary, 'r') as json_input_stream:
            reference_cpac_bids_path_mapping = json.load(json_input_stream)

        for cpac_path, bids_path in cpac_bids_path_mapping.items():
            assert bids_path == reference_cpac_bids_path_mapping[cpac_path]

        # eventually could add a test for bids compliance

    def test_everything(self):

        import cpac_output_to_bids as cpb
        import os

        cpac_data_config_file = os.path.dirname(__file__) + '/test_files/pre_skullstrip_data_config.yml'
        cpac_output_directory = os.path.dirname(__file__) + '/test_files/pipeline_cpac'

        cpb.main([cpac_data_config_file, cpac_output_directory, 'bids_dir', 'dry_run'])
