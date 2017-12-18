#!/usr/bin/env python
# -*- coding: utf-8

from __future__ import unicode_literals
import unittest
from os import path
from argparse import Namespace
from dorina import config


class TestConfig(unittest.TestCase):
    def setUp(self):
        self.original_config = config
        self.original_basedir = config._basedir
        self.original_default_name = config._default_name

    def tearDown(self):
        config._config = self.original_config
        config._basedir = self.original_basedir
        config._default_name = self.original_default_name

    def test_load_config(self):
        """Test config.load_config()"""
        config._basedir = path.dirname(__file__)
        c = config._default_name = 'test.cfg'
        c = Namespace(testing=True)

    def test_set_config(self):
        """Test config.set_config()"""
        c = Namespace(testing=True)

    def test_get_config(self):
        """Test config.get_config()"""
        config._config = Namespace(testing=True)
