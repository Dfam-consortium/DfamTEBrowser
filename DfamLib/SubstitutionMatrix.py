#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    Usage: from SubstitutionMatrix import SubstitutionMatrix

    A generic substitution matrix object

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


class SubstitutionMatrix:
    """
    A generic substitition matrix class

    A minimal set of functions and structures for holding a
    generic ( nucleotide, amino-acid, RNA ) scoring matrices.

    Example:
            SubstitutionMatrix( matrix = [ [  5, -1, -2, -3 ],
                                           [ -2,  3, -1, -1 ],
                                           [ -1, -1,  3, -2 ],
                                           [ -1, -2, -1,  4 ] ],
                                alphabet = [ 'A', 'C', 'G', 'T' ] )

       or
            SubstitutionMatrix( file = "25p41g.matrix" )

    Additional options (optional):

       Background Frequencies
           background_frequences : Assumed background composition for
                                   the scoring matrix.

            e.g. background_frequencies = {'A': 0.232, 'C': 0.213, ... }

       Gap Penalties
           This object stores all penalties for scoring gap events.  The
           convention here is that all values are considered a penalty
           and therefore they are stored as negative number regardless of
           how they are passed to the constructor.

           gap_penalty            : A single penalty to apply equally to each gap
                                    position.
        or
           gap_open_penalty       : Affine gap penalty. The penalty to apply to
                                    the first position of a contiguous gap.

           gap_extension_penalty  : Affine gap penalty. The penalty to apply to
                                    the 2nd and remaining positions of a gap.
        or
           ins_open_penalty       : The penalty to apply to the first position
                                    of a contiguous insertion.
           ins_extension_penalty  : The penalty to apply to the 2nd and remaining
                                    positions of an insertion.
           del_open_penalty       : The penalty to applly to the first position
                                    of a contiguous deletion.
           del_extension_penalty  : The penalty to apply to the 2nd and remaining
                                    positions of a deletion.

    TODO: Add in the ability to store matrix KA parameters calculated by ALP
          Should the alphabet be stored as unicode or byte strings?

    """

    matrix = []
    alphabet_r = []
    alphabet_h = {}
    background_freqs = None
    ins_gap_init = None
    ins_gap_extn = None
    del_gap_init = None
    del_gap_extn = None

    # Class static attributes
    freq_data_RE = re.compile(
        r"^#?\s*FREQS\s+(\S+)\s+([\d\.]+)\s+(\S+)\s+([\d\.]+)\s+(\S+)\s+([\d\.]+)\s+(\S+)\s+([\d\.]+)\s*$"
    )
    alphabet_data_RE = re.compile(r"^\s*([A-Za-z\s]+)$")
    matrix_data_RE = re.compile(r"^\s*[A-Za-z]?\s*([\-\d\s]+)$")

    def __init__(self, *args, **kwargs):
        allowed_keys = set([])
        self.__dict__.update((k, None) for k in allowed_keys)
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)

        # Instance attributes
        self.matrix = []
        self.alphabet_r = []
        self.alphabet_h = {}
        self.background_freqs = {}
        self.ins_gap_init = None
        self.ins_gap_extn = None
        self.del_gap_init = None
        self.del_gap_extn = None
        self._cached_lambda = None

        ## Required construction
        # Construct from matrix and alphabet arrays
        if "matrix" in kwargs and "alphabet" in kwargs:
            self.matrix = kwargs.get("matrix", None)
            # TODO: Weak validation - only the first row column count
            #       is validated against the number of rows.  Is there
            #       a quick way to assess if this is a N x N matrix?
            if len(self.matrix) != len(self.matrix[0]):
                raise Exception("ERROR: Not a N x N matrix!")
            self.alphabet_r = kwargs.get("alphabet", None)
            idx = 0
            for ch in self.alphabet_r:
                self.alphabet_h[ch] = idx
                idx += 1
            if len(self.alphabet_h) != len(self.matrix):
                raise Exception("ERROR: Alphabet size and matrix size mismatch!")
        # Construct from numpy/pandas data
        elif "matrix" in kwargs:
            print("Matrix")
        # Construct from file
        elif "file" in kwargs:
            self.readMatrixFromFile(kwargs.get("file", None))
        else:
            raise Exception("ERROR: 'matrix' or 'alphabet' required for construction")

        ## Handle optional data
        if "background_frequencies" in kwargs:
            if isinstance(kwargs.get("background_frequencies"), dict):
                self.background_freqs = kwargs.get("background_frequencies")
            else:
                print(
                    "ERROR: expected dictionary for background_frequencies, got "
                    + type(kwargs.get("background_frequencies"))
                )

        #
        if "gap_penalty" in kwargs:
            penalty = -abs(int(kwargs.get("gap_penalty")))
            self.ins_gap_init = penalty
            self.ins_gap_extn = penalty
            self.del_gap_init = penalty
            self.del_gap_extn = penalty
        elif "gap_open_penalty" in kwargs and "gap_extension_penalty" in kwargs:
            self.ins_gap_init = -abs(int(kwargs.get("gap_open_penalty")))
            self.ins_gap_extn = -abs(int(kwargs.get("gap_extension_penalty")))
            self.del_gap_init = -abs(int(kwargs.get("gap_open_penalty")))
            self.del_gap_extn = -abs(int(kwargs.get("gap_extension_penalty")))
        elif (
            "ins_open_penalty" in kwargs
            and "ins_extension_penalty" in kwargs
            and "del_open_penalty" in kwargs
            and "del_extension_penalty" in kwargs
        ):
            self.ins_gap_init = -abs(int(kwargs.get("ins_open_penalty")))
            self.ins_gap_extn = -abs(int(kwargs.get("ins_extension_penalty")))
            self.del_gap_init = -abs(int(kwargs.get("del_open_penalty")))
            self.del_gap_extn = -abs(int(kwargs.get("del_extension_penalty")))
        elif (
            "ins_open_penalty" in kwargs
            or "ins_extension_penalty" in kwargs
            or "del_open_penalty" in kwargs
            or "del_extension_penalty" in kwargs
            or "gap_open_penalty" in kwargs
            or "gap_extension_penalty" in kwargs
        ):
            print("ERROR: strange combination of gap penalty attributes!")

    def getGapPenalty(self):
        """
        Get the basic per-position gap penalty

        Return the average per-position gap penalty or 'None'.
            None

        Returns:
          The average of the stored gap parameters:
            average( ins_open_penalty, ins_extension_penalty,
                     del_open_penalty, del_extension_penalty )
          or None if gap penalties are not defined.
        """
        if (
            self.ins_gap_init != None
            and self.ins_gap_extn != None
            and self.del_gap_extn != None
            and self.del_gap_extn != None
        ):
            return int(
                (
                    self.ins_gap_init
                    + self.ins_gap_extn
                    + self.del_gap_init
                    + self.del_gap_extn
                )
                / 4
            )
        else:
            return None

    def getAffineGapPenalty(self):
        """
        Get the affine gap penalties

        Return the average (insertion/deletion) affine gap penalties or 'None'.

        Args:
            None

        Returns:
          The average open penalty and average extension penalty of the
          stored gap parameters:
            ( average( ins_open_penalty + del_open_penalty ),
              average( ins_extension_penalty, del_extension_penalty ) )
          or None if gap penalties are not defined.
        """
        if (
            self.ins_gap_init != None
            and self.ins_gap_extn != None
            and self.del_gap_extn != None
            and self.del_gap_extn != None
        ):
            return (
                int(((self.ins_gap_init + self.del_gap_init) / 2)),
                int(((self.ins_gap_extn + self.del_gap_extn) / 2)),
            )
        else:
            return None

    def getAffineInsDelPenalty(self):
        """
        Get the affine gap penalties

        Return the insertion and deletion affine gap penalties or 'None'.

        Args:
            None

        Returns:
          The insertion open/extension and deletion open/extension
          penalties or None if gap penalties are not defiled.
          The value returned is a list in the following order:
            ( ins_open_penalty, ins_extension_penalty,
              del_open_penalty, del_exentension_penalty )
        """
        if (
            self.ins_gap_init != None
            and self.ins_gap_extn != None
            and self.del_gap_extn != None
            and self.del_gap_extn != None
        ):
            return (
                self.ins_gap_init,
                self.ins_gap_extn,
                self.del_gap_init,
                self.del_gap_extn,
            )
        else:
            return None

    def getMatrixValue(self, query_char, subject_char):
        """
        Get the matrix value for a pair of characters

        Given the query character (typically the ancestral
        state) and the subject character (typically the
        derived or current state) for an aligned pair,
        return the matrix score.

        Args:
            query_char   :
            subject_char :

        Returns:
        """
        if query_char == None or subject_char == None:
            raise ("ERROR: Missing query or subject characters!")

        if not query_char in self.alphabet_h:
            raise (
                "ERROR: query_char '"
                + query_char
                + "' is not a valid symbol in this matrix!"
            )

        if not subject_char in self.alphabet_h:
            raise (
                "ERROR: subject_char '"
                + subject_char
                + "' is not a valid symbol in this matrix!"
            )
        return self.matrix[self.alphabet_h[query_char]][self.alphabet_h[subject_char]]

    def toString(self):
        """
        Generate a formatted matrix string from object

        This function converts the Matrix object into
        common output formats for use with external tools.

        Args:
            None

        Returns:
          A string or "" if the object does not contain
          data.
        """
        m_str = ""
        if self.background_freqs is not None:
            m_str = m_str + "FREQS "
            for ch in self.background_freqs:
                m_str = (
                    m_str
                    + "{:>2}".format(ch)
                    + " {:0.3f}".format(self.background_freqs[ch])
                )
            m_str = m_str + "\n"
        for ch in self.alphabet_r:
            m_str = m_str + "{:>5}".format(ch)
        m_str = m_str + "\n"
        for row in self.matrix:
            for col in row:
                # TODO: pre-calculate the max/min of the matrix
                #       and use that to set the format width
                m_str = m_str + "{:5.0f}".format(col)
            m_str = m_str + "\n"

        return m_str

    def saveMatrixToFile(self, matrixFilename):
        """
        Generate a standard matrix file from the Object

        Args:
            matrixFilename   :  The name of the file to
                                create.
        """
        with open(matrixFilename, "w") as matrix_file:
            matrix_file.write(self.toString())

    def readMatrixFromFile(self, matrixFilename):
        """
        Read matrix data from a file

        Matrix format conventions:

        Ancestral State (Rows) -> Current State (Columns)

        E.g given the matrix

               A   C   G   T
           A   8 -18 -10 -21
           C -17  12 -16  -7
           G  -7 -16  12 -17
           T -21 -10 -18   8

        the score for an assumed ancestral "A" mutating to a "G"
        would be -10.

        Crossmatch, RMBlast, WUBlast, ABBlast etc. all store
        matrices in this assumed order. DNA matrices may be formated
        minimally as:

               A   C   G   T
               8 -18 -10 -21
             -17  12 -16  -7
              -7 -16  12 -17
             -21 -10 -18   8

        Optionally, the format may include the row's designator although
        it is treated as a convenience label and is ignored.  The object
        requires at that row order matches the column order.

        In addition a single line with background base frequencies may
        be provided as in

            FREQS A 0.325 C 0.175 G 0.175 T 0.325
        or
            # FREQS A 0.285 C 0.215 G 0.215 T 0.285

        Finally, any number of "#" prefixed comment lines may be
        present.  These are ingored unless they contain a valid
        background frequency line.

        For example the following is a typical matrix used
        with RMBlast:

            # FREQS A 0.285 C 0.215 G 0.215 T 0.285
                A   R   G   C   Y   T   K   M   S   W   N   X
            A   9   3  -7 -18 -19 -21 -14  -4 -12  -5  -1 -30
            R   0   3   1 -18 -18 -19  -8  -9  -8 -10  -1 -30
            G -10  11  11 -18 -18 -18  -3 -14  -3 -14  -1 -30
            C -18 -18 -18  11  11 -10 -14  -3  -3 -14  -1 -30
            Y -19 -18 -18   2   0   0  -9  -8  -8 -10  -1 -30
            T -21 -19 -18  -7   3   9  -4 -14 -12  -5  -1 -30
            K -16  -9  -3 -12  -8  -4  -3 -14  -8 -10  -1 -30
            M  -4  -8 -12  -3  -9 -16 -14  -3  -8 -10  -1 -30
            S -14  -9  -3  -3  -9 -14  -9  -9  -3 -14  -1 -30
            W  -5  -9 -12 -12  -9  -5  -9  -9 -12  -5  -1 -30
            N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
            X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
        """
        # Clear current values
        self.matrix = []
        self.alphabet_r = []
        self.alphabet_h = {}
        self.background_freqs = None
        self.ins_gap_init = None
        self.ins_gap_extn = None
        self.del_gap_init = None
        self.del_gap_extn = None
        with open(matrixFilename, "r") as matrix_file:
            for line in matrix_file:
                freqs_match = SubstitutionMatrix.freq_data_RE.match(line)
                if freqs_match:
                    self.background_freqs = {}
                    # TODO: this is hardcoded for 4 freqs
                    for g_idx in range(1, 8, 2):
                        self.background_freqs[freqs_match.group(g_idx)] = float(
                            freqs_match.group(g_idx + 1)
                        )
                    continue

                # Ignore remaining comments
                if line.startswith("#"):
                    continue

                # Match alphabet line
                alpha_chars = self.alphabet_data_RE.match(line)
                if alpha_chars:
                    self.alphabet_r = alpha_chars.group(1).split()
                    idx = 0
                    for ch in self.alphabet_r:
                        self.alphabet_h[ch] = idx
                        idx += 1
                    continue

                mat_data = self.matrix_data_RE.match(line)
                if mat_data:
                    if self.alphabet_r is None or len(self.alphabet_r) == 0:
                        raise Exception("Missing matrix column labels!")
                    row_data = mat_data.group(1).split()
                    row_data_ints = []
                    for score_str in row_data:
                        row_data_ints.append(int(score_str))
                    self.matrix.append(row_data_ints)

    def getLambda(self):
        """
        Calculates the value of lambda for this substitution matrix, used
        for calculating complexity-adjusted scores.

        The calculation is based on Matrix.pm from RepeatMasker, which
        was itself based on Phil Green's swat/cross_match programs.

        Hazard: getLambda() is relatively expensive to calculate, so it is cached.
        If values on the matrix are set manually after calculation,
        the old lambda value will be returned.
        """

        if self._cached_lambda is not None:
            return self._cached_lambda

        # "lambda" is a keyword in python!
        lmbda_upper = 0
        lmbda_lower = 0
        lmbda = 0.5

        freqs = self.background_freqs
        alph_h = self.alphabet_h
        matrix = self.matrix

        # S = The sum, over all pairs of bases A and B in the frequencies list, of:
        #         freq(A) * freq(B) * e^(lambda * score(A, B))
        def calc_S(try_lmbda):
            nonlocal freqs, alph_h, matrix

            S = 0
            check = 0
            for (base_a, freq_a) in freqs.items():
                idx_a = alph_h[base_a]
                for (base_b, freq_b) in freqs.items():
                    idx_b = alph_h[base_b]
                    S += freq_a * freq_b * math.exp(try_lmbda * matrix[idx_a][idx_b])
                    check += freq_a * freq_b

            if check < 0.999 or check > 1.001:
                # This sum-of-products of frequencies should equal 1, since
                # the sum of frequencies should have equaled 1
                raise Exception(
                    "failed to calculate matrix lambda: sanity check failed"
                    + " ({} should have been close to 1)".format(check)
                )

            return S

        # Goal: Find a value for lambda such that S = 1.

        # First, we will find a pair of numbers (lmbda_lower, lmbda_upper) such
        # that the desired lambda value most be somewhere between these two
        # values.

        # Then, we use binary search to find the closest possible value of
        # lambda within a tolerance threshold.

        S = 0
        while S < 1:
            S = calc_S(lmbda)
            if S < 1:
                lmbda_lower = lmbda
                lmbda = lmbda * 2
        lmbda_upper = lmbda

        while lmbda_upper - lmbda_lower > 0.00001:
            lmbda = (lmbda_lower + lmbda_upper) / 2

            S = calc_S(lmbda)
            if S >= 1:
                lmbda_upper = lmbda
            else:
                lmbda_lower = lmbda

        self._cached_lambda = lmbda
        return lmbda

    def __repr__(self):
        """
        __repr__() - Generic representation of an object in JSON format
        """
        return json.dumps(self.__dict__, indent=4)
