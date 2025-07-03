#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    Usage: from Sequence import Sequence

    A generic sequence object

SEE ALSO: Dfam: http://www.dfam.org

AUTHOR(S):
    Robert Hubley <rhubley@systemsbiology.org>

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

"""

#
# Module imports
#
import sys
import os
import math
import string
import json
import re
import itertools
import tempfile
import numpy as np


def _xstr(s):
    """
    Cast variables to strings with Null handling

    This is a helper function for casting variables as
    strings. Null values are converted to the string '(Null)'.
    """
    return "(Null)" if s is None else str(s)


class SequenceEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return str(obj.view("c"))
        return json.JSONEncoder.default(self, obj)


class Sequence:
    """
    A generic sequence class

    A generic object to hold an immutable biological sequence
    ( DNA, Amino Acid, or RNA ) and optional metadata.

    Example:
            Sequence('ACCTTAAAGGGC')
        or
            Sequence('ACCTTAAAGGGC', id='seq1')
        or
            Sequence('ACCTTAAAGGGC', id='seq1', start_pos=1, end=12)

        or, if applicable:

            Sequence('ACCTTAAAGGGC', id='seq1', start_pos=1, end=12, orient='+')

    Parameters
    ----------

        sequence  :
              id  :
           start  :
             end  :
          orient  :

    """

    comp_map = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "Y": "R",
        "R": "Y",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
        "X": "X",
        "-": "-",
        ".": ".",
        " ": " ",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
        "y": "r",
        "r": "y",
        "s": "s",
        "w": "w",
        "k": "m",
        "m": "k",
        "b": "v",
        "d": "h",
        "h": "d",
        "v": "b",
        "n": "n",
        "x": "x",
    }

    _sequence = None
    _left_flanking_sequence = None
    _right_flanking_sequence = None
    id = None
    start_pos = None
    end_pos = None
    orient = None

    def __init__(self, sequence, *args, **kwargs):
        allowed_keys = set([])
        self.__dict__.update((k, None) for k in allowed_keys)
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)

        self._sequence = self._to_byte_str(sequence)

        if "left_flanking_sequence" in kwargs:
            self._left_flanking_sequence = self._to_byte_str(
                kwargs.get("left_flanking_sequence")
            )

        if "right_flanking_sequence" in kwargs:
            self._right_flanking_sequence = self._to_byte_str(
                kwargs.get("right_flanking_sequence")
            )

        if "id" in kwargs:
            self.id = kwargs.get("id", None)

        if "start_pos" in kwargs and "end_pos" in kwargs:
            self.start_pos = int(kwargs.get("start_pos", None))
            self.end_pos = int(kwargs.get("end_pos", None))
        elif "start_pos" in kwargs or "end_pos" in kwargs:
            raise Exception(
                "Optional Sequence range attributes 'start_pos' "
                "and 'end_pos' must both be specified.  Only one found."
            )
        if "orient" in kwargs:
            self.orient = kwargs.get("orient", None)

    def _to_byte_str(self, seq):
        if isinstance(seq, np.ndarray):
            if seq.dtype == "|S1":
                seq = seq.view(np.uint8)
                if seq.shape == ():
                    seq = np.array([seq], dtype=np.uint8)
            # elif seq.dtype == '|U1' or seq.dtype == '<U1':
            #    np.array([ch.astype(int) for ch in seq])
            elif seq.dtype != np.uint8:
                raise TypeError(
                    "Cannot automatically convert numpy.ndarray of dtype: " + seq.dtype
                )
            if not seq.flags["C_CONTIGUOUS"]:
                seq = np.ascontiguousarray(seq)
        elif isinstance(seq, str):
            seq = seq.encode("ascii")
            seq = np.frombuffer(seq, dtype=np.uint8)
        else:
            raise TypeError(
                "Required numpy.ndarray of dtype " "np.uint8, or python 'str' data."
            )
        return seq

    @property
    def sequence(self):
        return self._sequence.tostring().decode()

    @property
    def sequence_as_byte_array(self):
        return self._sequence.view("|S1")

    @property
    def sequence_as_unicode_array(self):
        return self._sequence.view("|S1").astype("U")

    def complement(self, reverse=True):
        lookup = np.zeros(256, dtype=np.uint8)
        for key, value in self.comp_map.items():
            lookup[ord(key)] = ord(value)
        c_seq = lookup[self._sequence]
        if reverse:
            c_seq = c_seq[::-1]
        self._sequence = c_seq
        if self.orient == "+":
            self.orient = "-"
        elif self.orient == "-":
            self.orient = "+"

    def __str__(self):
        """
        Generate a readable summary for the Sequence

        Args:
            None

        Returns:
          A string or "" if the object does not contain
          data.
        """
        m_str = "Sequence:\n"
        if self.id is not None:
            m_str = m_str + "  id       = '" + self.id + "'\n"
        if self.start_pos is not None or self.end_pos is not None:
            m_str = (
                m_str
                + "  coords   = "
                + _xstr(self.start_pos)
                + "-"
                + _xstr(self.end_pos)
                + "\n"
            )
        if self.orient is not None:
            m_str = m_str + "  orient   = " + self.orient + "\n"
        if self._sequence.size:
            m_str = m_str + "  length   = " + _xstr(self._sequence.size) + "\n"
            m_str = (
                m_str
                + "  alphabet = "
                + ",".join([chr(item) for item in set(self._sequence)])
                + "\n"
            )

            # Determine character counts -- TODO: this could be case
            # sensitive. E.g it might be useful to know the number of 'c's vs
            # 'C's in a sequence.
            ucase_seq = [ord(chr(item).upper()) for item in self._sequence]
            freqs = np.bincount(ucase_seq)
            (indices,) = np.nonzero(freqs)
            chars = indices.astype(np.uint8).tostring().decode("ascii")
            obs_counts = freqs[indices]
            char_cnts = dict(zip(chars, obs_counts.tolist()))
            m_str = m_str + "  counts   = " + str(char_cnts) + "\n"

            # Display the sequence or sequence ends
            if self._sequence.size > 20:
                m_str = (
                    m_str
                    + "  sequence = '"
                    + "".join([chr(item) for item in self._sequence[1:10]])
                    + "'..'"
                    + "".join([chr(item) for item in self._sequence[-10:]])
                    + "'"
                )
            else:
                m_str = (
                    m_str
                    + "  sequence = '"
                    + "".join([chr(item) for item in self._sequence])
                    + "'"
                )
        return m_str

    def __repr__(self):
        """
        __repr__() - Generic representation of an object in JSON format
        """
        return json.dumps(self.__dict__, indent=4, cls=SequenceEncoder)
