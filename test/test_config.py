#!/usr/bin/env python
# -*- coding: utf-8

from __future__ import unicode_literals
import unittest
try:
    import mock
except ImportError:
    import unittest.mock as mock
from os import path
from dorina.config import load_configuration, validate_configuration


class TestConfig(unittest.TestCase):
    def setUp(self):
        basedir = path.dirname(path.abspath(__file__))
        name = 'test.cfg'
        self.load_configuration = load_configuration
        self.validate_configuration = validate_configuration
        self.config_path = path.join(basedir, name)

    def tearDown(self):
        self.load_configuration = None
        self.validate_configuration = None
        self.config_path = None

    def test_load_config(self):
        """Test load_configuration()"""
        config = self.load_configuration(self.config_path)
        self.assertEqual(config.get('DEFAULT', 'data_path'), '/test/data')

    def test_set_config(self):
        """Test set()"""
        config = self.load_configuration(self.config_path)
        config.set('DEFAULT', 'test', value='testing!')
        self.assertEqual(config.get('DEFAULT', 'test'), 'testing!')

    def test_validate_config(self):
        """Test config.validate_configuration"""
        config = self.load_configuration(self.config_path)
        val_config = self.validate_configuration(config)
        self.assertEqual(config.get('DEFAULT', 'data_path'),
                         val_config.get('DEFAULT', 'data_path'))

        config.set('DEFAULT', 'data_path', value='~')
        val_config = self.validate_configuration(config)
        self.assertNotEqual(val_config.get('DEFAULT', 'data_path'),
                            '~')