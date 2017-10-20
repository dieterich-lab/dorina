#!/usr/bin/env python
# -*- coding: utf-8
import logging
from os import path

import sys
from six.moves import configparser

from dorina.utils import expand_local_path

__all__ = ['config']
_basedir = path.dirname(path.abspath(__file__))
_default_name = 'default.cfg'


log = logging.getLogger(__name__)
logging.captureWarnings(True)
logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s ')


def load_configuration():
    configuration = configparser.ConfigParser()
    with open(path.join(_basedir, _default_name)) as f:
        try:
            configuration.read_file(f)
        except AttributeError:  # py27 support
            configuration.readfp(f)
    return configuration


def validate_configuration(configuration):
    default_path = configuration['DEFAULT'].get('data_path')
    if '~' in default_path:
        configuration['DEFAULT']['data_path'] = expand_local_path(default_path)

    return configuration


config = load_configuration()
config = validate_configuration(config)
#
#
# def load_config(namespace):
#     """Load config, but don't overwrite existing settings"""
#
#     if 'configfile' in namespace:
#         default_file = namespace.configfile
#     else:
#         default_file = path.join(_basedir, _default_name)
#
#     config = configparser.ConfigParser()
#     with open(default_file, 'r') as fp:
#         try:
#             config.read_file(fp)
#         except AttributeError:  # py27 support
#             config.readfp(fp)
#
#     for s in config.sections():
#         if not s in namespace:
#             namespace.__dict__[s] = Namespace()
#         for key, value in config.items(s):
#             key = key.replace('-', '_')
#             if not key in namespace.__dict__[s]:
#                 namespace.__dict__[s].__dict__[key] = value
#
#     # settings from the [DEFAULT] section need to be handled extra
#     for key, value in config.items('DEFAULT'):
#         key = key.replace('-', '_')
#         if not key in namespace:
#             namespace.__dict__[key] = value
#
#
# def set_config(namespace):
#     """Set global configuration"""
#     global _config
#     _config = namespace
#
#
# def get_config():
#     """Get global configuration"""
#     global _config
#     return _config
