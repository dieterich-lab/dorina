#!/usr/bin/env python
# -*- coding: utf-8

from os import path
from argparse import Namespace
try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser

_config = None
_basedir = path.dirname(path.abspath(__file__))
_default_name = 'default.cfg'


def load_config(namespace):
    """Load config, but don't overwrite existing settings"""

    if 'configfile' in namespace:
        default_file = namespace.configfile
    else:
        default_file = path.join(_basedir, _default_name)

    config = ConfigParser.ConfigParser()
    with open(default_file, 'r') as fp:
        try:
            config.read_file(fp)
        except AttributeError:  # py27 support
            config.readfp(fp)

    for s in config.sections():
        if not s in namespace:
            namespace.__dict__[s] = Namespace()
        for key, value in config.items(s):
            key = key.replace('-', '_')
            if not key in namespace.__dict__[s]:
                namespace.__dict__[s].__dict__[key] = value

    # settings from the [DEFAULT] section need to be handled extra
    for key, value in config.items('DEFAULT'):
        key = key.replace('-', '_')
        if not key in namespace:
            namespace.__dict__[key] = value


def set_config(namespace):
    """Set global configuration"""
    global _config
    _config = namespace


def get_config():
    """Get global configuration"""
    global _config
    return _config
