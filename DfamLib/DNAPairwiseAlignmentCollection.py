#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    Usage:

    A collection of DNAPairwiseAlignment objects and related functions

SEE ALSO:
          Dfam: http://www.dfam.org

AUTHOR(S):
    Sean Rice <srice@systemsbiology.org>

LICENSE:
    This code may be used in accordance with the Creative Commons
    Zero ("CC0") public domain dedication:
    https://creativecommons.org/publicdomain/zero/1.0/

DISCLAIMER:
  This software is provided ``AS IS'' and any express or implied
  warranties, including, but not limited to, the implied warranties of
  merchantability and fitness for a particular purpose, are disclaimed.
  In no event shall the authors or the Dfam consortium members be
  liable for any direct, indirect, incidental, special, exemplary, or
  consequential damages (including, but not limited to, procurement of
  substitute goods or services; loss of use, data, or profits; or
  business interruption) however caused and on any theory of liability,
  whether in contract, strict liability, or tort (including negligence
  or otherwise) arising in any way out of the use of this software, even
  if advised of the possibility of such damage.

  NOTE: This is an early draft - more functions to come

"""

#
# Module imports
import sys
import os
import math
import string
import subprocess
import datetime
import json
import re
import itertools
import tempfile

# Dfam Library
sys.path.append(os.path.dirname(__file__))
from DNAPairwiseAlignment import DNAPairwiseAlignment


class DNAPairwiseAlignmentCollection:
    """ """

    def __init__(self, *args, **kwargs):
        allowed_keys = set(["debug"])
        self.__dict__.update((k, None) for k in allowed_keys)
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)
        self._results = []

    def __iter__(self):
        """ """
        return iter(self._results)

    def append(self, result):
        if not isinstance(result, DNAPairwiseAlignment):
            raise TypeError("New result is not a DNAPairwiseAlignment object")
        self._results.append(result)

    def pop(self):
        popped = self._results.pop()
        return popped

    def extend(self, results):
        for result in results:
            if not isinstance(result, DNAPairwiseAlignment):
                raise TypeError("New result is not a DNAPairwiseAlignment object")
        self._results.extend(results)

    def __len__(self):
        return len(self._results)

    def len(self):
        return len(self)

    def sort(self, key=None, reverse=True):
        """
        In Place Sort
            e.g resColl.sort(key=lambda r: r.bit_score, reverse=True)
        """
        if key is None:
            self._results = sorted(
                self._results, key=lambda r: r.e_value, reverse=reverse
            )
        else:
            self._results = sorted(self._results, key=key, reverse=reverse)

    def get(self, index):
        return self._results[index]

    def __getitem__(self, index):
        return self._results[index]

    def getCrossmatchOutput(
        self,
        show_alignment=False,
        out_file_format=False,
        reference="model",
        outputDir=None,
    ):
        """
        This function simply creates a string of all the collection's crossmatch
        output data. It accepts the same args as the individual alignment crossmatch
        output function, with the addition of an optional outputDir, which would
        save the output in a file in the specified directory
        """
        outputStr = ""
        for alignment in self._results:
            outputStr += alignment.crossmatch_encode(
                show_alignment=show_alignment,
                out_file_format=out_file_format,
                reference=reference,
            )
        if outputDir != None:
            if os.path.isdir(outputDir) == False:
                raise ValueError(
                    "Error: output directory was specified but could not be found: "
                    + outputDir
                )
            newFileName = "DNAPAC" + str(datetime.datetime.now()) + ".out"
            newFile = os.path.join(outputDir, newFileName)
            newOutputFile = open(newFile, "x")
            newOutputFile.write(outputStr)
            newOutputFile.close()
        return outputStr

    def __repr__(self):
        """
        __repr__() - Generic representation of an object in JSON format
        """
        return json.dumps(self.__dict__, indent=4)
