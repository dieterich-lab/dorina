#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 09:32 10/10/2017 2017

"""
import logging
from ftplib import FTP
import os
from os import path

import requests

try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve  # py27

from dorina.utils import check_file_extension, uncompress
from dorina.config import config

log = logging.getLogger('dorina.config')


# def unzip_and_write_file(data, data_path, mode='wb'):
#     if data_path.endswith('.gz'):
#         with open(data_path.replace('.gz', ''), mode) as f:
#             f.write(fileobj=data)
#     else:
#         raise NotImplementedError('{} extension is not supported.'.format(
#             data_path.rsplit('.', 1)[-1]
#         ))


class EnsemblFTP(object):
    """
    Base template for retrieving file and information from Ensembl FTP server.
    """

    def __init__(self, version=config['DEFAULT'].get('version'),
                 organism=config['DEFAULT'].get('organism'), *args, **kwargs):
        self.base_url = 'ftp.ensemblorg.ebi.ac.uk'
        self.url = []
        self.local_data = config['DEFAULT'].get('data_path')
        self.available = {}
        self.version = version
        self.organism = organism
        self.assembly = None
        self.file_name = None

    def check_available(self, url=None):
        """
        List the files and directories available at the FTP server.

        :param str url: location to be listed
        """
        if not url.startswith('/'):
            url = '/' + url

        with FTP(self.base_url) as ftp:
            ftp.login()
            # TODO throw error in URL unavailable
            available_dir = ftp.nlst(url)

        self.available[url] = available_dir

    def check_local(self, file_name):

        return path.isfile(path.join(self.local_data, self.assembly, file_name))

    # def ftp_fetch(self, remote_file):
    #     remote_path, file_name = remote_file.rsplit('/', 1)
    #
    #     with FTP(self.base_url) as ftp:
    #         ftp.login()
    #         ftp.cwd(remote_path)
    #
    #         def callback(data):
    #             return unzip_and_write_file(data, self.local_data + file_name)
    #
    #         ftp.retrbinary("RETR {}".format(file_name), callback)

    def retrieve_file_from_url(self, url):
        """
        Retrieve a file from FTP server.
        .. note:: urlretrieve has a better API for error handling than ftp.retrbinary

        :param str url: file url
        """

        remote_path, file_name = url.rsplit('/', 1)
        # if self.check_local(file_name):
        #     pass  # TODO
        if self.assembly is None:
            self.assembly = self.get_assembly()

        local_path = path.join(self.local_data, self.assembly, file_name)
        try:
            os.makedirs(path.join(self.local_data + self.assembly))
        except OSError as e:
            pass  # dir exists
        urlretrieve('ftp://' + self.base_url + url, local_path)
        uncompress(local_path)
        # TODO RETURN PATH

    def fetch_all(self):
        for file in self.url:
            self.retrieve_file_from_url(file)

    def check_extension(self, extension):
        """
        Check for error in the built urls
        :param str extension: file extension
        :raise ValueError: if the file extension is not in url
        """

        self.url = [url for url in self.url
                    if check_file_extension(url, extension)]

    def get_assembly(self, version=None, organism=None):
        """
        This is a hack to get the Assembly for a Ensembl release and organism.
        For currently release, just use the REST .
        :param version:
        :param organism:
        :return:
        """
        if version is None:
            version = self.version
        if organism is None:
            organism = self.organism
        url = u'/pub/{}/gff3/{}/'.format(version, organism)
        self.check_available(url)
        self.url.append(self.available[url][-2])
        return self.url[0].split('.')[1]

    def retrieve_from_tsv_by_index(self, extension='tsv.gz', file_index=5):
        """
        Retrieves files from the TSV directory from Ensembl TSV.
        -  CHECKSUM
        -  Compara ncrna homologies
        -  Compara protein homologies
        -  Ena mappings
        -  Entrez mappings
        -  Karyotype (similar to UCSC genome files)
        -  refseq mappings
        -  Uniprot Uniprot
        -  README Ena
        -  README Entrez
        -  README refseq
        -  README Uniprot

        Default is the karyotype file with chromosome's lengths
        """
        url = u'/pub/{}/tsv/{}/'.format(self.version, self.organism)
        self.check_available(url)
        self.url.append(self.available[url][file_index])
        self.check_extension(extension)
        # self.fetch_all()


class EnsemblGFF(EnsemblFTP):
    """
    Retrieves GFF3 files from Ensembl FTP server.
    For the gff3 directory structure is:
    - CHECKSUMS
    - abinitio file
    - chromosome file
    - a file for each chromosome
    - a complete file with every chromosome
    - README

    Default is the GFF3 file with all chromosomes
    """

    def __init__(self,
                 version=config['DEFAULT'].get('version'),
                 organism=config['DEFAULT'].get('organism'),
                 file_index=-2, *args, **kwargs):
        EnsemblFTP.__init__(self, *args, **kwargs)
        extension = 'gff3'

        url = u'/pub/{}/{}/{}/'.format(version, extension, organism)
        self.check_available(url)
        self.url.append(self.available[url][file_index])
        self.assembly = self.get_assembly(version, organism)
        self.check_extension(extension)
        self.fetch_all()


class EnsemblVariation(EnsemblFTP):
    """
    Retrieves Regulation VCS files from Ensembl Variation FTP server.
    The regulation directory structure is:
    -  CHECKSUM
    -  Peaks
    -  QualityCheck
    -  RegulatoryFeatureActivity
    -  Regulatory features
    -  Motif features

    .. note:
    Only available for human and mice.
    .. todo: better document each data type.
    """

    def __init__(self, version=config['DEFAULT'].get('version'),
                 organism=config['DEFAULT'].get('organism'), extension='vcf.gz',
                 *args, **kwargs):
        EnsemblFTP.__init__(self, *args, **kwargs)

        self.url.append(u'/pub/{}/variation/vcf/{}/{}.{}'.format(
            version, organism, organism.capitalize(), extension))
        # TODO others type of variation files
        self.check_extension(extension)
        self.fetch_all()


class EnsemblRegulation(EnsemblFTP):
    """
    Retrieves Regulation bed files from Ensembl FTP server.
    The regulation directory structure is:
    -  CHECKSUM
    -  Peaks
    -  QualityCheck
    -  RegulatoryFeatureActivity
    -  Regulatory features
    -  Motif features

    .. note:
    Only available for human and mice.
    .. TODO: Document each data type.
    """

    def __init__(self, version=config['DEFAULT'].get('version'),
                 organism=config['DEFAULT'].get('organism'),
                 tissue=config['DEFAULT'].get('tissue'), extension='bed', experiment='*',
                 *args, **kwargs):
        EnsemblFTP.__init__(self, *args, **kwargs)

        url = u'/pub/{}/regulation/{}/Peaks/{}/{}'.format(
            version,
            organism,
            tissue,
            experiment)
        # TODO should we support other type of experiments ?

        if '*' in url:
            url = url.replace('*', '')
        self.check_available(url)
        for _url in self.available[url]:
            self.check_available(_url)
        del self.available[url]
        [self.url.extend(x) for x in self.available.values()]
        self.check_extension(extension)
        self.fetch_all()


class EnsemblRest(object):
    """
    Base class for retrieving data from Ensembl Rest API.

    Limited support to previous Ensembl release using:
    >>> rest = EnsemblRest(base_url='grch37.rest.ensembl.org')
    >>> rest.ping()
    """

    def __init__(self, *args, **kwargs):
        self.base_url = kwargs.get('base_url', 'http://rest.ensembl.org/')
        self.headers = kwargs.get('headers', {'content-type': 'application/json'})
        self.reqs_per_sec = 15
        self.session = requests.Session()

    def ping(self):
        ext = "info/ping?"

        return self.get(ext)

    def get_info_assembly(self, organism):
        ext = 'info/assembly/{}/X?'.format(organism)
        return self.get(ext)

    def get_genetree_members_by_id(self, ensembl_id, params=None):
        """
        Retrieves genetree for a given Ensembl stable id.

        :param str ensembl_id: Ensembl stable id
        :param dict params:

        :return dict: genetree result

        ..note::
        This endpoint supports isoform identification for protein ids, such as
        ENSP00000288602.6, but not for transcripts, as ENST00000288602.10

        """
        # params = {'aligned': 0, 'callback': '-', 'cigar_line': 0, 'clusterset_id': '-',
        #          'db_type': '-', 'nh_format': 'simple', 'object_type': [], 'prune_species': '-',
        #          'prune_taxon': '-', 'sequence': 'protein', 'species': '-'}
        if '.' in ensembl_id:
            ensembl_id = ensembl_id.split('.')[0]
            log.info(
                'Isoform identifier detected, falling back to stable ID {}'.format(ensembl_id))
        ext = 'genetree/member/id/{}?'.format(ensembl_id)
        return self.get(ext, params=params)

    def get_genomic_alignment(self):
        pass  # todo

    def get_phenotype_by_region(self, organism, chromosome, start, end):
        """
        ..Note Maximum region length is 5 Mb.
        :param str organism: Organism scientific name, such as homo_sapiens
        :param int chromosome: Chromosome number, eg: 9
        :param int start: First base pair of the query in genomic region (0 based index)
        :param int end: Last base of the query in the genomic region (inclusive)
        :return dict:
        """
        ext = "/{}/{}:{}-{}?".format(organism, chromosome, start, end)
        return self.get(ext)

    def get(self, ext, json=True, headers=None, params=None):
        """

        :param ext:
        :param json:
        :param headers:
        :param params:
        :return:
        """
        # TODO limit request number and use requests.Session for retry
        if headers is None:
            headers = self.headers
        r = requests.get(self.base_url + ext, headers=headers, params=params)
        if not r.ok:
            r.raise_for_status()

        if json:
            return r.json()
        else:
            return r.content


class ExpressionAtlas(object):
    """
    Retrieve differential gene expression experiment from the EBI Gene atlas
    """
    pass
