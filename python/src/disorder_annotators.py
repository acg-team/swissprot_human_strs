#!/usr/bin/env python
"""
Classes designed to obtain MobiDB-Lite disorder annotations for all proteins in a specified fasta file.
Annotations are either calculated locally (LocalDisorderAnnotator), or retrieved from MobiDB
through the API (APIDisorderAnnotator). The former is written specifically to function in the MobiDB-Lite Docker created
from the image at https://github.com/BioComputingUP/MobiDB-lite_docker.

Author: Max Verbiest
Contact: max.verbiest@zhaw.ch
"""

import subprocess
import os
import json
import urllib.request
import urllib.error

try:
    from Bio import SeqIO
except ModuleNotFoundError:
    print("WARNING: Could not load Bio.SeqIO -> can only perform local disorder annotations")

__all__ = [
    "DisorderAnnotator",
    "APIDisorderAnnotator",
    "LocalDisorderAnnotator",
]


class DisorderAnnotator(object):
    """Abstract class to outline functionality.
    Subclasses implement methods for annotating disorder in protein sequences, extracted from
    a specified fasta file.
    """

    def __init__(self, fasta, out_file):
        """
        Parameters
        fasta (str):        path to fasta file to extract protein sequences from
        out_file (str):     path where output file will be created (directory structure must exist already)
        """

        self.fasta = fasta
        self.out_file = self.check_output_dir(out_file)

    def check_output_dir(self, out_file):
        """Check whether the directory where the output file will be generated exists. If not -> FileNotFoundError"""

        out_dir = os.path.join(*out_file.split("/")[:-1])
        if not os.path.isdir(out_dir):
            raise FileNotFoundError("Directory for output file does not exist")
        return out_file

    def get_disorder_annotations(self):
        """Implemented in child classes"""

        pass


class APIDisorderAnnotator(DisorderAnnotator):
    """
    Class to query MobiDB for all proteins in a specified fasta file. Protein identifiers must be in standard
    UniProt/ SwissProt format. Protein IDs are queried to MobiDB and predicted disordered regions (from MobiDB-Lite)
    are extracted and written to an output file.

    For information on MobiDB API, see:
            https://mobidb.bio.unipd.it/help/apidoc
        and
            https://mobidb.bio.unipd.it/about/mobidb
    """

    def __init__(self, fasta, out_file, curated=False):
        """
        Parameters
        curated (bool): (Not implemented!) Should curated information from
                        MobiDB be incorporated in output? (default:False)
        """

        super().__init__(fasta, out_file)
        self.fasta = SeqIO.parse(fasta, "fasta")
        self.base_url = "https://mobidb.bio.unipd.it/api/download?acc={}{}&format=json"
        self.filter = "&projection=prediction-disorder-mobidb_lite"
        if curated:
            # self.filter += ",curated-disorder-priority"
            raise NotImplementedError("Handling of curated MobiDB annotations is not yet supported")

    def query_mobidb(self, prot_id):
        """Construct url for current protein id and filter, query MobiDB

        Parameters
        prot_id (str):  UniProt identifier of protein of interest

        Returns
        disorder_dict (dict):
                        Dictionary containing all information available for protein entry in MobiDB after filter.
                        If only part of this information is of interest, it can be extracted from the dictionary later
        """

        prot_url = self.base_url.format(prot_id, self.filter)
        with urllib.request.urlopen(prot_url) as response:
            disorder = response.read().decode("utf-8")
            disorder_dict = json.loads(disorder)
        return disorder_dict

    def get_disorder_annotations(self):
        """ Retrieve MobiDB-Lite annotations through MobiDB web API.
        Implements method DisorderAnnotator.get_disorder_annotations(). Output is generated to mimic MobiDB-Lite output
        when run with the '-f interpro' and '-sf' options, and is written to an output file.
        """

        with open(self.out_file, "w") as o:
            for record in self.fasta:
                prot_id = record.id.split("|")[1]
                print("Started work on sequence '{}'...".format(prot_id))
                try:
                    disorder_dict = self.query_mobidb(prot_id)
                except urllib.error.URLError:
                    print("{} was not found in MobiDB".format(prot_id))
                    continue

                try:
                    disorder_regions = [(i[0], i[1]) for i in
                                        disorder_dict["prediction-disorder-mobidb_lite"]["regions"]]
                except KeyError:
                    print("No disorder predicted".format(prot_id))
                    continue
                print("Predicted disordered regions: {}".format(disorder_regions))
                for coords in disorder_regions:
                    o.write("{}\t{}\t{}\n".format(record.id, coords[0], coords[1]))


class LocalDisorderAnnotator(DisorderAnnotator):
    """
    Class to run MobiDB-lite disorder annotations locally. Will only work using the docker image provided
    on https://github.com/BioComputingUP/MobiDB-lite_docker
    """

    def __init__(self, fasta, out_file, out_format="interpro", silence_sf=True):
        """
        Parameters
        out_format (str):   Output format, define output detail (MobiDB-lite gitbub for details). Accepted values
                            are: "interpro (default)", "fasta", "caid", "mobidb4"
        silence_sf (bool):  If True (default), sequence features (e.g. polar, polyampholyte) will be silenced in output
        """

        super().__init__(fasta, out_file)
        # path to mobidb-lite in the docker
        self.out_format = self.check_out_format(out_format)
        self.silence_sf = silence_sf
        self.mobidb_path = "/usr/src/mobidb/mobidb_lite.py"

    def check_out_format(self, format):
        """Check whether the specified output format is supported by MobiDB-lite"""

        if not format in {"interpro", "fasta", "caid", "mobidb4"}:
            raise ValueError("Please pick one of the following output formats: interpro, fasta, caid, mobidb4")
        return format

    def run_mobidb(self, threads=1, env=os.environ.copy()):
        """Run MobiDB-Lite script, method taken from:
                            https://github.com/BioComputingUP/MobiDB-lite_docker/blob/master/test.ipynb

        Parameters:
        threads (int):     Number of threads involved in disordered regions computation range between 1(default) and 7
        env (dict):        Environmental variables (e.g. python path) which must be set in order
                           to correctly run MobiDB Lite script
        """

        # Call subprocess
        return subprocess.run(
            check=True,  # Check command execution
            encoding='utf-8',  # Set stdout/stderr encoding
            env=env,  # Set environmental variables
            stdout=subprocess.PIPE,  # Capture stdout
            stderr=subprocess.PIPE,  # Capture stderr
            # Command line arguments
            args=[
                # Run script with Python 3
                'python3', self.mobidb_path,
                # Set output file format
                '-f', '{:s}'.format(self.out_format),
                # Set output file
                '-o', '{:s}'.format(self.out_file),
                # Set number of threads, if any
                *(['-t', '{:d}'.format(threads)] if threads else []),
                # Silence sequence features, if desired
                *(['{:s}'.format("-sf")] if self.silence_sf else []),
                # Set input file path
                '{:s}'.format(self.fasta)
            ]
        )

    def get_disorder_annotations(self, threads=1):
        """Run MobiDB-Lite locally and report output/ errors
        Implements method DisorderAnnotator.get_disorder_annotations() Annotations are written to output file.

        Parameters
        threads (int):   How many threads MobiDB-Lite should use: range between 1(default) and 7
        """

        try:
            # Run MobiDB Lite
            ran = self.run_mobidb(threads)
            print('MobiDB Lite succesfully completed')
            # If annotations were not written to file: retrieve stdout and stderr
            if ran.stdout:
                print('stdout:')
                print(ran.stdout)
            if ran.stderr:
                print('stderr:')
                print(ran.stderr)
            print()
        except subprocess.CalledProcessError as err:
            # Show error
            print('MobiDB Lite exited with code', err.returncode)
            print('command:')
            print(err.cmd)
            print('stdout:')
            print(err.stdout)
            print('stderr:')
            print(err.stderr)
            print()
