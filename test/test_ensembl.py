#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 07:41 16/10/2017 2017 

"""
from __future__ import unicode_literals

import unittest

try:
    from unittest import mock
except ImportError:
    import mock

from nose.tools import assert_true, raises
from requests.exceptions import HTTPError

from dorina.ensembl import EnsemblRest


class TestEnsemblRest(unittest.TestCase):
    def setUp(self):
        self.ensembl_rest = EnsemblRest()
        self.ensembl_rest_bad = EnsemblRest(base_url='http://xxx.org')

    def TearDown(self):
        self.ensembl_rest = None
        self.ensembl_rest_bad = None

    @mock.patch('dorina.ensembl.urlretrieve')
    def test_ensembl_ftp(self, mock_urlretrieve):
        urlretrieve = mock_urlretrieve()
        urlretrieve.return_value = []

    @mock.patch('dorina.ensembl.requests.get')
    def test_ensembl_rest_ping(self, mock_get):
        """Test the communication with the Rest server"""
        mock_get.return_value = mock.Mock(ok=True)

        assert_true(self.ensembl_rest.ping())

    @mock.patch('dorina.ensembl.requests.get')
    @raises(HTTPError)
    def test_ensembl_bad_ping(self, mock_get):
        """Test the communication bad server"""
        mock_get.return_value = mock.Mock(ok=False)
        mock_get.side_effect = HTTPError(mock.Mock(status=404), 'not found')
        self.ensembl_rest_bad.ping()

    @mock.patch('dorina.ensembl.requests.get')
    def test_get_info_assembly_w_defaults(self, mock_get):
        normal_response = {'is_chromosome': 1,
                           'length': 156040895,
                           'is_circular': 0,
                           'assembly_exception_type': "REF",
                           'assembly_name': "GRCh38",
                           'coordinate_system': "chromosome"}
        mock_get.return_value = mock.Mock(ok=True)
        mock_get.return_value.json.return_value = normal_response

        assembly = self.ensembl_rest.get_info_assembly('human')
        assert_true(assembly, 'GRCh38')

    def test_get_genetree_members(self):
        pass


class TestEnsemblFTP(unittest.TestCase):
    pass


if __name__ == '__main__':
    import nose

    nose.runmodule()
