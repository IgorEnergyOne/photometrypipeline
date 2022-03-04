#!/usr/bin/env python
# yaml ConfigParser
import yaml

def read_yml_config(path: str) -> dict:
    """
    reads YAML config into a dictionary
    """
    with open(r'{}'.format(path)) as configfile:
        config = yaml.load(configfile, Loader=yaml.FullLoader)

        for key, value in config.items():
            print (key + " : " + str(value))
    return config

def parse_config(config_dict: dict) -> list:
    """
    """
    config = read_yml_config(path = 'pp_config.yml')
    # pp_run parameters
    prefix = config['pp_run'].get('prefix')
    target = config['pp_run'].get('target')
    filter = config['pp_run'].get('filter')
    fixed_aprad = config['pp_run'].get('fixed_aprad')
    source_tolerance = config['pp_run'].get('source_tolerance')
    solar = config['pp_run'].get('solar')
    rerun_regist = config['pp_run'].get('rerun_registration')
    asteroids = config['pp_run'].get('asteroids')
    reject = config['pp_run'].get('reject')
    keep_wcs = config['pp_run'].get('keep_wcs')


if __name__ == '__main__':
    config = read_yml_config(path = 'pp_config.yml')
    # pp_run parameters
    prefix = config['pp_run'].get('prefix')
    target = config['pp_run'].get('target')
    filter = config['pp_run'].get('filter')
    fixed_aprad = config['pp_run'].get('fixed_aprad')
    source_tolerance = config['pp_run'].get('source_tolerance')
    solar = config['pp_run'].get('solar')
    rerun_regist = config['pp_run'].get('rerun_registration')
    asteroids = config['pp_run'].get('asteroids')
    reject = config['pp_run'].get('reject')
    keep_wcs = config['pp_run'].get('keep_wcs')
    print(prefix)
