import pkg_resources
import os
from shutil import copyfile


def get_path_of_data_file(data_file):
    file_path = pkg_resources.resource_filename("popsynth", "data/%s" % data_file)

    return file_path


def get_path_of_data_dir():
    file_path = pkg_resources.resource_filename("popsynth", "data")

    return file_path


def copy_package_data(data_file):

    data_file_path = get_path_of_data_file(data_file)
    copyfile(data_file_path, "./%s" % data_file)

def copy_template():

    copy_package_data('template_config.yaml')
    
    
