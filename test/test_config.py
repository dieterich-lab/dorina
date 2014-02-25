# vim: set fileencoding=utf-8 :

import unittest
from os import path
from argparse import Namespace
from dorina import config


class TestConfig(unittest.TestCase):
    def setUp(self):
        self.original_config = config._config
        self.original_basedir = config._basedir
        self.original_default_name = config._default_name

    def tearDown(self):
        config._config = self.original_config
        config._basedir = self.original_basedir
        config._default_name = self.original_default_name

    def test_load_config(self):
        "Test config.load_config()"
        config._basedir = path.dirname(__file__)
        config._default_name = 'test.cfg'
        c = Namespace(testing=True)
        config.load_config(c)
        self.assertTrue(c.testing)
        self.assertEqual('true', c.default_file_loaded)
        self.assertEqual('true', c.sublevel.exists)

    def test_set_config(self):
        "Test config.set_config()"
        c = Namespace(testing=True)
        self.assertIsNone(config._config)
        config.set_config(c)
        self.assertEqual(c, config._config)

    def test_get_config(self):
        "Test config.get_config()"
        config._config = Namespace(testing=True)
        c = config.get_config()
        self.assertEqual(config._config, c)
