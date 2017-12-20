#!/usr/bin/env python
# -*- coding: utf-8
import logging
from copy import copy
from os import path

import sys
from six.moves import configparser

from dorina.utils import validate_data_path

__all__ = ['config']
_basedir = path.dirname(path.abspath(__file__))
_default_name = 'default.cfg'


log = logging.getLogger(__name__)
logging.captureWarnings(True)
logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s ')


def load_configuration(config_path=path.join(_basedir, _default_name)):
    configuration = configparser.ConfigParser()
    with open(config_path) as f:
        try:
            configuration.read_file(f)
        except AttributeError:  # py27 support
            configuration.readfp(f)
    return configuration


def validate_configuration(configuration):
    default_path = configuration.get('DEFAULT', 'data_path')
    if '~' in default_path:
        configuration.set(
            'DEFAULT', 'data_path', value=validate_data_path(default_path))

    return configuration


config = load_configuration()
config = validate_configuration(config)
