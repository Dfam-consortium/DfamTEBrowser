#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    Usage: from MultAlign import MultAlign

    A multiple sequence alignment object based on
    MultAln.pm and built up from MultipleSeqAlign.py

SEE ALSO: Dfam-js Project - MultAlign
             ( https://github.com/Dfam-consortium/Dfam-js )
          Dfam: http://www.dfam.org

AUTHOR(S):
    Robert Hubley <rhubley@systemsbiology.org>
    Sean Rice <srice@systemsbiology.org> (It was a pleasure!)
    Jeb Rosen <jrosen@systemsbiology.org>
    Original perl implementation and concepts by Arnie Kas and Arian Smit

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
import warnings
import tempfile
import numpy as np
from Sequence import Sequence
from DNAPairwiseAlignmentCollection import DNAPairwiseAlignmentCollection
from SubstitutionMatrix import SubstitutionMatrix
import datetime
import ConsensusCaller


def _xstr(s):
    """
    Cast variables to strings with Null handling

    This is a helper function for casting variables as
    strings. Null values are converted to the string '(Null)'.
    """
    return "(Null)" if s is None else str(s)


class StockholmParseError(Exception):
    """Exception class for errors in parsing stockholm files."""

    def __init__(self, message):
        self.message = message


class MultAlign:
    """
    A generic Multiple Sequence Alignment (MSA) class

    An optimized class for holding the minimal set of information
    needed to store, interrogate, and call the consensus on a
    MSA.  The intention here is to use this as a base class for
    the DfamSeedAlignment class **eventually**.  For now it's
    a utility class.

    Example:
            MultAlign( sequences = [ '  AACT-TTGACC---CC',
                                            'GAAAGT-TTCTCCAGTCC',
                                            'G-AACTATTC-CCA--CC',
                                            'GAAACT-TTG-CC     '],
                              reference =   'xxxxxx.xxx.xx...xx' )

    NEW:
            MultAlign( sequences = [ '  AACT-TTGACC---CC',
                                            'GAAAGT-TTCTCCAGTCC',
                                            'G-AACTATTC-CCA--CC',
                                            'GAAACT-TTG-CC     '])

        or


            MultAlign( sequences = [
                 ['  AACT-TTGACC---CC', 'seq1'],
                 ['GAAAGT-TTCTCCAGTCC', 'seq2'],
                 ['G-AACTATTC-CCA--CC', 'seq3'],
                 ['GAAACT-TTG-CC     ', 'seq4']
                                          ],
                              reference =   'xxxxxx.xxx.xx...xx' )

        or

            MultAlign( sequences = [
                 ['  AACT-TTGACC---CC', 'seq1', 10, 21],
                 ['GAAAGT-TTCTCCAGTCC', 'seq2', 15, 31],
                 ['G-AACTATTC-CCA--CC', 'seq3', 1, 14],
                 ['GAAACT-TTG-CC     ', 'seq4', 71, 81]
                                          ],
                              reference =   'xxxxxx.xxx.xx...xx' )

         or

            MultAlign( sequences = [
                 ['  AACT-TTGACC---CC', 'seq1', 10, 21, '+'],
                 ['GAAAGT-TTCTCCAGTCC', 'seq2', 15, 31, '-'],
                 ['G-AACTATTC-CCA--CC', 'seq3', 1, 14, '-'],
                 ['GAAACT-TTG-CC     ', 'seq4', 71, 81, '-']
                                          ] )

      Where sequences must all be the same length, contain upper
      case DNA IUB codes, '-' gap alignment symbols within sequences
      and use ' ' (space) padding at start/end unaligned regions.

      Sequnce metadata is currently stored in a parallel array of dicts, the keys
      to which can be found below. Some data (such as id) will be required,
      or result in an error if omitted, and other data may be omitted and None
       will be used as a placeholder in the dict

    Note: src_start/end = alignedSeqStart (from MultAln):
            where in the original sequence the sequence is located: values range
            from 0 to the length of the original sequence
          aligned_start/end = aligned_start/end (from MultAln)
            the relative position (to the reference) possible values range from
            0 to length of the reference

    TODO: Extend to support a generic Sequence Object based on np arrays
          Add optional - MSA ID
          Add Move A2M/A3M routines here from DfamSeedAlignment
          Create generic file import/export for alignment formats
          Arian's LINUP format

    """

    alphabet_r = ["A", "R", "G", "C", "Y", "T", "K", "M", "S", "W", "N", "X", "Z", "-"]
    alphabet_h = {
        "A": 0,
        "R": 1,
        "G": 2,
        "C": 3,
        "Y": 4,
        "T": 5,
        "K": 6,
        "M": 7,
        "S": 8,
        "W": 9,
        "N": 10,
        "X": 11,
        "Z": 12,
        "-": 13,
    }
    # Static Regular Expressions
    alignDataRE = re.compile("^(\S+)\s+([\.\-ACGTRYKMSWBDHVNXacgtrykmswbdhvnx]+)\s*$")

    # store/update consensus
    reference = None

    def __init__(self, *args, **kwargs):
        """
        Optional Args:

            DNAPAC: accepts a DNAPairwiseAlignmentCollection object to be parsed
                    into multAlign object. DNAPAC MUST be comprised of DNA Pairwise
                    Alignments where many sequences are run agains 1 reference,
                    see _alignFromDNAPAC() for more info.*

            MSF: Accepts a string path to an MSF file to be parsed into the MultAlign
                    oject. See _alignFromMSF() for more info.*

            Stockholm: Accepts a string path to a stockholm file to be parsed into
                        the MultAlign object. See _alignFromMSF() for more info.*

            FASTA: Accepts a string path to a FASTA file to be parsed into the
                    MultAlign object. See _alignFromFASTA() for more info.*

            sequences: accepts a list of sequences, which must be pre-aligned by
                        the caller, and which must all be of type str or Sequence.
                        These sequences must be of the same length as the reference.
                        This allows the user to directly provide sequences to the
                        object rather than parsing, etc.*

        ->  !NOTE!: only 1 of the above (starred) arguments can be used upon construction
                    You cannot pass a FASTA and a stockholm, the object will
                    only be constructed from 1 file. If multiple are passed, the
                    object will be constructed from whichever arg appears first
                    in this list: DNAPAC, MSF, Stockholm, FASTA, sequences

            seqMetaData: accepts a list of python dicts which must be the same length
                            as the list of sequences provided, and the metadata
                            dicts must be indexed parallel in their lists to the
                            sequences they correspond to. Any metadata can be stored,
                            (though it may not necessarily be used by the object.)
                            The metadata that is (currently) stored and used by the
                            object includes:
                                left_flanking: str
                                right_flanking: str
                                name: str
                                src_start: int
                                src_end: int
                                src_orientation: str (MUST be either "+" or "-")
                                src_transitions: float
                                src_transversions: float
                                src_gc_background: float TODO: verify
                                aligned_start: int
                                aligned_end: int

                            Providing a list of metadata dicts that use these name
                            but do not have the types listed above will not end well.
                            NOTE: this argument will only be accepted if the
                            'sequences' argument was also provided.

            reference: accepts a single sequence of type str or Sequence. Must be
                        same length as and pre-aligned to the sequences provided.
                        NOTE: this argument will only be accepted if 'sequences'
                        argument was ALSO provided.

            refMetaData: accepts a single dict that will be used to hold the reference
                        meta data. Same rules apply to this dict as to those in
                        seqMetaData argument - recall it's just 1 dict this time
                        though.
                        NOTE: this argument will only be accepted if the 'reference'
                        arg was ALSO provided

            checkForIllegalChars: function accepts a boolean (default true) that
                                    will determine if reference and aligned sequences
                                    will be checked for illegal chars which include
                                    chars that are non IUB codes or unrecognized
                                    gap chars. (recognized gaps are ".", " ", "-")

            internalGapChar: accepts a single character that will be used to represent
                            gaps that are internal to the sequence (aka that occur
                            after the first and before the last non-gap base pair).
                            Recognized gaps include space, dash, and period
                            (" ", "-", "."). Any operations that involve sequence
                            or reference creation/storage/manipulation will use this
                            character to represent internal gaps.
                            Default is dash ("-")
                            Having this set once for the entire object prevents
                            performing operations that would cause different seqs
                            to have different gap characters

            externalGapChar:similar to internalGapChar -
                            accepts a single character that will be used to represent
                            gaps that are external to the sequence (aka that occur
                            before the first and after the last non-gap base pair).
                            These gaps are solely used to align the sequences and
                            reference within the underlying 2d array and are not
                            part of the original parsed sequence.
                            Recognized gaps are the same as for internalGapChar.
                            Default is space (" ")
                            ________________________________________________
                            |    ACGTCGTA--TCTCTC-CTCTC--CTCTCT----C       |
                            |^external gaps             ^internal gaps     |
                            |______________________________________________|

        TODO:
        --  currently sequence metadata (like flank seqs) are NOT stored as numpy
            arrays of chars while it's tbd if they should be or not :
            if they are changed to chars, serializeInput and output will need to
            be altered to reflect that

        --  Add additional explanatory comments - many could probably just be copied
            and pasted from MultAln.pm

        --  Currently seq and ref meta data args are only interpreted in case where
            user supplies the ref and seqs directly - keep?
        """
        # Here we feed the allowed optional args into the object...
        allowed_keys = set(
            [
                "checkForIllegalChars",
                "internalGapChar",
                "externalGapChar",
                "verifyMetaData",
            ]
        )
        self.__dict__.update((k, None) for k in allowed_keys)
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)

        # Robert has requested that minimum data/metadata requirements be removed,
        # so we will set these here and if they are never set down the line by
        # constructor / user, at least they exist.
        self.refernce = None
        self.aligned_sequences = None
        self.refMetaData = None
        self.seqMetaData = None

        # check for internal / external gap arg validity - or set to default
        if self.internalGapChar == None:
            self.internalGapChar = "-"
        elif self.internalGapChar not in ("-", ".", " "):
            raise ValueError(
                "MultAlign Error: currently only period, space, and dash"
                + " ('.', ' ', '-') are acccepted as intenalGapChar"
            )
        if self.externalGapChar == None:
            self.externalGapChar = " "
        elif self.externalGapChar not in ("-", ".", " "):
            raise ValueError(
                "MultAlign Error: currently only period, space, and dash"
                + " ('.', ' ', '-') are acccepted as externalGapChar"
            )

        if "DNAPAC" in kwargs:
            self._alignFromDNAPAC(
                kwargs.get("DNAPAC"),
                alignment_reference=kwargs.get("alignment_reference", None),
                debug=kwargs.get("debug", False),
                reference=kwargs.get("reference"),
            )
        elif "MSF" in kwargs:
            self._alignFromMSF(
                kwargs.get("MSF"),
                debug=kwargs.get("debug", False),
                firstAsRef=kwargs.get("firstAsRef", False),
            )

        elif "Stockholm" in kwargs:
            self._alignFromStockholm(
                kwargs.get("Stockholm"), debug=kwargs.get("debug", False)
            )

        elif "FASTA" in kwargs:
            self._alignFromFASTA(kwargs.get("FASTA"), debug=kwargs.get("debug", False))

        # if user gives us nothing to parse, manual entry is also allowed...
        elif "sequences" in kwargs:
            # Sequences should already be padded to be the same length
            seqs = kwargs.get("sequences", None)

            # It's faster to build out a python 2d array and then convert it to a
            # ndarray than to call ndarray.append()
            aligned_data = []
            if "seqMetaData" not in kwargs:
                seqMetaData = []
            else:
                seqMetaData = kwargs.get("seqMetaData")
            for seq in seqs:
                if isinstance(seq, Sequence):
                    lFlank = None
                    rFlank = None
                    aligned_data.append(
                        [ch.upper() for ch in seq.sequence_as_unicode_array]
                    )

                    if "seqMetaData" not in kwargs:
                        if isinstance(seq._left_flanking_sequence, np.ndarray):
                            lFlank = "".join(
                                seq._left_flanking_sequence.view("|S1")
                                .astype("U")
                                .tolist()
                            )
                        if isinstance(seq._right_flanking_sequence, np.ndarray):
                            rFlank = "".join(
                                seq._right_flanking_sequence.view("|S1")
                                .astype("U")
                                .tolist()
                            )
                        metaDict = {
                            "left_flanking": lFlank,
                            "right_flanking": rFlank,
                            "name": seq.id,
                            "src_start": int(seq.start_pos),
                            "src_end": int(seq.end_pos),
                            "src_orientation": seq.orient,
                            "src_transitions": None,
                            "src_transversions": None,
                            "src_gc_background": None,
                            "aligned_start": None,
                            "aligned_end": None,
                        }

                elif isinstance(seq, str):
                    aligned_data.append(list(seq.upper()))

                    if "seqMetaData" not in kwargs:
                        metaDict = {
                            "left_flanking": None,
                            "right_flanking": None,
                            "name": None,
                            "src_start": None,
                            "src_end": None,
                            "src_orientation": None,
                            "src_transitions": None,
                            "src_transversions": None,
                            "src_gc_background": None,
                            "aligned_start": None,
                            "aligned_end": None,
                        }

                else:
                    raise ValueError(
                        "MultAlign Error: Currently only string and"
                        + " sequence types are supported for the 'sequences'"
                        + " argument."
                    )

                if "seqMetaData" not in kwargs:
                    seqMetaData.append(metaDict)

            self.aligned_sequences = np.asarray(aligned_data)
            self.seqMetaData = seqMetaData
            # Toss out our garbage
            aligned_data = None
            seqs = None

        if (
            "reference" in kwargs
            and "DNAPAC" not in kwargs
            and "MSF" not in kwargs
            and "FASTA" not in kwargs
            and "Stockholm" not in kwargs
        ):
            referenceHolder = kwargs.get("reference", None)
            if isinstance(referenceHolder, str):
                self.reference = np.array(list(kwargs.get("reference", None).upper()))
                if "refMetaData" not in kwargs:
                    self.refMetaData = {
                        "left_flanking": None,
                        "right_flanking": None,
                        "name": None,
                        "src_start": None,
                        "src_end": None,
                        "src_orientation": None,
                        "src_transitions": None,
                        "src_transversions": None,
                        "src_gc_background": None,
                    }
            elif isinstance(referenceHolder, Sequence):
                lFlank = None
                rFlank = None
                self.reference = np.array(
                    [ch.upper() for ch in referenceHolder.sequence_as_unicode_array]
                )
                if "refMetaData" not in kwargs:
                    if isinstance(referenceHolder._left_flanking_sequence, np.ndarray):
                        lFlank = "".join(
                            seq._left_flanking_sequence.view("|S1").astype("U").tolist()
                        )
                    if isinstance(referenceHolder._right_flanking_sequence, np.ndarray):
                        rFlank = "".join(
                            referenceHolder._right_flanking_sequence.view("|S1")
                            .astype("U")
                            .tolist()
                        )
                    self.refMetaData = {
                        "left_flanking": lFlank,
                        "right_flanking": rFlank,
                        "name": referenceHolder.id,
                        "src_start": int(referenceHolder.start_pos),
                        "src_end": int(referenceHolder.end_pos),
                        "src_orientation": referenceHolder.orient,
                        "src_transitions": None,
                        "src_transversions": None,
                        "src_gc_background": None,
                    }

            if "refMetaData" in kwargs:
                self.refMetaData = kwargs.get("refMetaData")

        # Here we will validate that all sequences and the reference have the same
        # length (if we actually got them). They may or may not provide the reference
        # and aligned sequences which will inform our comparison
        if self.aligned_sequences is not None:
            if self.reference is not None:
                comparator = len(self.reference)
                compName = "reference"
            else:
                comparator = len(self.aligned_sequences[0])
                compName = (
                    "0th aligned Sequence"  # we may not have the name so this will do
                )

            for i, seq in enumerate(self.aligned_sequences):
                if len(seq) != comparator:
                    raise ValueError(
                        "Error: MultAlign Construction Error: \n"
                        + "After construction, the Length of aligned sequence "
                        + str(i)
                        + ", length "
                        + str(len(seq))
                        + "did not "
                        + "match the length of the "
                        + compName
                        + ", of "
                        + "length "
                        + str(comparator)
                        + ".\n"
                    )

        # We do not allow any non-IUB characters (unless user explicitly specifies to)
        if self.checkForIllegalChars != False:
            self.illegalCharacterChecker()

        # The reason for different truth equality checks here is that we always
        # check for illegal characters unless explicitly told not to, and never
        # check the metadata unless explicitly told to. These values will be
        # None if not provided by the user
        if self.verifyMetaData == True:
            self.verifyMetaData()

    # ===========================================================================

    def reverseComplement(self, seq):
        """
        Simple (helper) function to reverse complement a single sequence
        """
        complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
        bases = list(seq)
        bases = reversed([complement.get(base, base) for base in bases])
        bases = "".join(bases)
        return bases

    # ===========================================================================

    def metaDataVerifier(self, warn=True):
        """
        TODO: formally test...

        A helper function that verifies that all metadata is of the correct type.
        Function will NOT be called by default, but users may use argument
        (verifyMeta = True) to call the function on the object after its construction
        is completed, or at any point during MultAlign use / manipulation.
        The function will warn the user if no metadata is present regardless of
        argument passed to warn.


        Args:
            warn :  accepts a boolean that determines if mistyped metadata throws
                    warnings or errors. Warnings will be used by default.
        """

        typeDict = {
            "left_flanking": str,
            "right_flanking": str,
            "name": str,
            "src_start": int,
            "src_end": int,
            "src_orientation": str,
            "src_transitions": int,
            "src_transversions": int,
            "src_gc_background": int,
            "aligned_start": int,
            "aligned_end": int,
        }

        if self.seqMetaData != None:
            for i, metaDict in enumerate(self.seqMetaData):
                if metaDict["src_orientation"] not in ("+", "-", None):
                    if warn:
                        warnings.warn(
                            "\nMultAlign Warning: src orientation at seq "
                            + "index "
                            + str(i)
                            + " was not valid orientation. "
                            + "only '+', '-', and None are valid orientations!\n"
                        )
                    else:
                        raise ValueError(
                            "\nMultAlign Error: src orientation at seq "
                            + "index "
                            + str(i)
                            + " was not valid orientation. "
                            + "only '+', '-', and None are valid orientations!\n"
                        )

                for key in typeDict:
                    if type(metaDict[key]) != typeDict[key] and metaDict[key] != None:
                        if warn:
                            warnings.warn(
                                "\nMultAlign Warning: seq meta data for seq at "
                                + "index "
                                + str(i)
                                + " has type mismatch: "
                                + "meta data for '"
                                + key
                                + "' should be type "
                                + str(typeDict[key])
                                + ", but got type "
                                + str(type(metaDict[key]))
                                + " instead!\n"
                                + "This may cause some functions to fail "
                                + "or give unexpected results!\n"
                            )
                        else:
                            raise ValueError(
                                "\nMultAlign Error: seq meta data for seq at "
                                + "index "
                                + str(i)
                                + " has type mismatch: "
                                + "meta data for '"
                                + key
                                + "' should be type "
                                + str(typeDict[key])
                                + ", but got type "
                                + str(type(metaDict[key]))
                                + " instead!\n"
                                + "This may cause some functions to fail "
                                + "or give unexpected results!\n"
                            )
        else:
            warnings.warn(
                "\nMultAlign Warning: absolutely no seq meta data currently"
                + " stored in MultAlign object - this may cause some functions"
                + " to fail or give unexpected results!\n"
            )

        if self.refMetaData != None:
            if self.refMetaData["src_orientation"] not in ("+", "-", None):
                if warn:
                    warnings.warn(
                        "\nMultAlign Warning: ref orientation "
                        + "was not valid orientation. "
                        + "only '+', '-', and None are valid orientations!\n"
                    )
                else:
                    raise ValueError(
                        "\nMultAlign Warning: ref orientation "
                        + "was not valid orientation. "
                        + "only '+', '-', and None are valid orientations!\n"
                    )
            for key in typeDict:
                if (
                    type(self.refMetaData[key]) != typeDict[key]
                    and self.refMetaData[key] != None
                ):
                    if warn:
                        warnings.warn(
                            "\nMultAlign Warning: Ref Meta Data"
                            + " has type mismatch: "
                            + "meta data for '"
                            + key
                            + "' should be type "
                            + str(typeDict[key])
                            + ", but got type "
                            + str(type(self.refMetaData[key]))
                            + " instead!\n"
                            + "This may cause some functions to fail "
                            + "or give unexpected results!\n"
                        )
                    else:
                        raise ValueError(
                            "\nMultAlign Error: Ref Meta Data"
                            + " has type mismatch: "
                            + "meta data for '"
                            + key
                            + "' should be type "
                            + str(typeDict[key])
                            + ", but got type "
                            + str(type(self.refMetaData[key]))
                            + " instead!\n"
                            + "This may cause some functions to fail "
                            + "or give unexpected results!\n"
                        )
        else:
            warnings.warn(
                "\nMultAlign Warning: absolutely no ref meta data currently"
                + " stored in MultAlign object - this may cause some functions"
                + " to fail or give unexpected results!\n"
            )

    # ===========================================================================

    def _alignFromMSF(self, MSF, debug=False, firstAsRef=True):

        """
        Helper function to create a MultAlign object from an msf file.
        Arguments:
            MSF: string path to file

            firstAsRef: boolean arg (default True) to determine if the sequence
            appearing 1st in the file should be stored as the reference. If false,
            the 1st sequence scanned will become the 0th aligned sequence, and there
            will be no reference

        TODO:
        --  As with all of these _alignFromX(), we could implement additional
            metadata capturing, or we could just continue ditching it all....
        --  Are all guarenteed to be forward (+) or...?
        --  Current regex for capturing sequence lines does NOT expect white space
            to be used as a space character - keep?
        --  I assumed that the gaps would be periods or dashes, other gaps will
            not be properly replaced with the user's desired gap chars - do
            we need any other gap chars covered?
        --  Use 'identifier' instead of 'name' below in the code and warnings.
        """

        if not os.path.isfile(MSF):
            raise ValueError(
                "Error: MSF argument to _alignFromMSF is not a valid"
                + " File path! Path given: \n"
                + str(MSF)
                + "\n"
            )

        # matches !!AA_MULTIPLE_ALIGNMENT 1.0
        fileTypeLine = r"!!AA_MULTIPLE_ALIGNMENT (?P<vNum>\d.\d)"

        # matches: Name: S11448 Len: 743 Check: 3635 Weight: 1.00
        # currently only captures name, as that's all we keep from said line...
        # nameLine = r"^\s*Name\s*:\s*(?P<name>\S+)."
        nameLine = "\s*Name\s*:\s*(?P<name>[\S ]+)Len:"

        # matches: MIR#SINE/MIR  A.BCD..E.FGH...Ix..J.KxxL......... etc.
        # Note: not currently designed to expect spaces as gap chars
        seqLine = r"^\s*(?P<name>\S+)\s+(?P<seq>[a-zA-Z. -]+)\s+"

        correctFileType = False
        firstTimeAround = True
        dupesFound = False
        foundNames = False
        seqMetaData = []
        nameTracker = {}
        alignedSeqs = []
        index = 0

        with open(MSF, "r") as MSF:
            for strLine in MSF:
                if re.match(fileTypeLine, strLine):
                    if re.match(fileTypeLine, strLine).groupdict()["vNum"] != "1.0":
                        warnings.warn(
                            "MultAlign alignFromMSF warning: version detected"
                            + "other than 1.0 - parser was designed for 1.0 and "
                            + " problems may arise parsing other versions\n"
                        )
                        correctFileType = True

                elif re.match(nameLine, strLine):
                    if firstTimeAround and not correctFileType:
                        warnings.warn(
                            "MultAlign alignFromMSF warning: No fileType"
                            + " line detected ('!!!!AA_MULTIPLE_ALI...') - "
                            + "file may be wrong type...\n"
                        )
                    firstTimeAround = False
                    foundNames = True

                    parsedName = (
                        re.match(nameLine, strLine).groupdict()["name"].rstrip(" ")
                    )
                    if parsedName not in nameTracker:
                        nameTracker[parsedName] = 1
                        seqName = parsedName
                    else:
                        nameTracker[parsedName] += 1
                        dupesFound = True
                        seqName = parsedName + str(nameTracker[parsedName])

                    metaDict = {
                        "left_flanking": None,
                        "right_flanking": None,
                        "name": seqName,
                        "src_start": None,
                        "src_end": None,
                        "src_orientation": "+",
                        "src_transitions": None,
                        "src_transversions": None,
                        "src_gc_background": None,
                        "aligned_start": None,
                        "aligned_end": None,
                    }
                    seqMetaData.append(metaDict)
                    alignedSeqs.append("")

                elif re.match(seqLine, strLine) and foundNames:
                    alignedSeqs[index % len(alignedSeqs)] += re.match(
                        seqLine, strLine
                    ).groupdict()["seq"]
                    index += 1

        if dupesFound:
            warnings.warn(
                "Multalign _alignFromMSF() Warning: duplicates names "
                + "detected - all duplicate names will have a counter appended "
                + "starting at 2, i.e. 'Rob, Jeb, Rob2, Arian, Rob3...'"
                + " within the MA metadata\n"
            )

        # FIXME: places that have similar gap char change code may not function if
        # they were written as for seq in seqs instead of using indices
        for i in range(len(alignedSeqs)):
            length = len(alignedSeqs[i])
            if length != len(alignedSeqs[0]):
                raise ValueError(
                    "MultAlign _alignFromMSF() Error: sequence #"
                    + i
                    + ", named "
                    + seqMetaData[i]["name"]
                    + ", was not"
                    + " the same length as seq 0"
                )

            alignedSeqs[i] = (
                alignedSeqs[i]
                .replace("-", self.internalGapChar)
                .replace(".", self.internalGapChar)
            )
            alignedSeqs[i] = (
                alignedSeqs[i]
                .rstrip(self.internalGapChar)
                .ljust(length, self.externalGapChar)
            )
            alignedSeqs[i] = (
                alignedSeqs[i]
                .lstrip(self.internalGapChar)
                .rjust(length, self.externalGapChar)
            )
            alignedSeqs[i] = list(alignedSeqs[i])

        if firstAsRef:
            reference = alignedSeqs.pop(0)
            # print("\n=========================")
            # print(reference)
            refMetaData = seqMetaData.pop(0)
            reference = np.asarray(list(reference))
            # print(reference)
            self.reference = reference
            self.refMetaData = refMetaData

        self.aligned_sequences = np.asarray(alignedSeqs)
        self.seqMetaData = seqMetaData

    # ===========================================================================

    def toMSF(
        self,
        fileName=None,
        debug=False,
        replaceGapsWith=".",
        includeReference=False,
        descriptionLine=None,
    ):

        """
        Function to output contents of MultAlign object as file in msf format.
        Args:
            fileName: accepts desired name of file (will be named 'noNameGiven'
                        if no arg given)

            replaceGapsWith: character to replace all sequence gap characters with -
                            generally '.' by convention but others are allowable

            includeReference: if the reference is a sequence (or having the consensus
                                is desirable) user can opt to include the reference
                                in the MSF output as 1st sequence

            descriptionLine: accepts a string that will be inserted into the 'description'
                                section of the msf file

        """
        g = replaceGapsWith

        msfOutput = "!!AA_MULTIPLE_ALIGNMENT 1.0\n"
        if descriptionLine != None:
            msfOutput += descriptionLine.rstrip("\n")  # just to prevent doubles
        msfOutput += "\n\n"

        length = len(self.reference)
        now = datetime.datetime.now()
        now = now.strftime("%m/%d/%Y_%H:%M:%S")
        maxNameLen = 0

        if fileName == None:
            fileName = "NoNameGiven"

        if type(fileName) == str:
            msfOutput += (
                fileName
                + "  MultAlign.msf MSF: "
                + str(length)
                + "   Type: N   "
                + now
                + " Check: 0 ..\n\n"
            )
        else:
            raise ValueError("MultAlign Error: fileName arg to toMSF() must be str!\n")

        if includeReference:
            refName = self.refMetaData["name"]
            if type(refName) != str:
                warnings.warn(
                    "MultAlign toMSF() warning: includeReference arg set"
                    + " to true, but ref name not given or of incorrect type -"
                    + " 'Reference' will be substituted for name"
                )
                refName = "Reference"
            msfOutput += (
                " Name: "
                + refName
                + "    Len: "
                + str(length)
                + "    Check: 0    Weight: 1.00\n"
            )
            maxNameLen = len(refName)

        egc = self.externalGapChar
        for i, seq in enumerate(self.aligned_sequences):
            msfOutput += (
                " Name: "
                + self.seqMetaData[i]["name"]
                + "    Len: "
                + str(len("".join(seq).lstrip(egc).rstrip(egc)))
                + "    Check: 0    Weight: 1.00\n"
            )
            if len(self.seqMetaData[i]["name"]) > maxNameLen:
                maxNameLen = len(self.seqMetaData[i]["name"])

        msfOutput += "\n//\n\n"

        chunkSize = 50
        lineStart = 0
        lineEnd = chunkSize

        while lineEnd < length:
            if includeReference:
                seq = (
                    "".join(self.reference[lineStart:lineEnd])
                    .replace(" ", g)
                    .replace("-", g)
                )
                name = ("ref:" + refName)[0:maxNameLen]
                msfOutput += name + "  " + seq + "\n"

            for i, seq in enumerate(self.aligned_sequences):
                msfOutput += self.seqMetaData[i]["name"].ljust(maxNameLen) + "  "
                msfOutput += (
                    "".join(seq)[lineStart:lineEnd].replace(" ", g).replace("-", g)
                    + "\n"
                )
            lineStart += chunkSize
            lineEnd += chunkSize
            msfOutput += "\n"

        if includeReference:
            seq = "".join(self.reference[lineStart:]).replace(" ", g).replace("-", g)
            name = ("ref:" + refName)[0:maxNameLen]
            msfOutput += name + "  " + seq + "\n"

        for i, seq in enumerate(self.aligned_sequences):
            msfOutput += self.seqMetaData[i]["name"].ljust(maxNameLen) + "  "
            msfOutput += "".join(seq)[lineStart:].replace(" ", g).replace("-", g) + "\n"

        # TODO: after you test, add the file options
        if fileName != None:
            outputFile = open(fileName, "w")
            outputFile.write(msfOutput)
            outputFile.close()
        return msfOutput

    # ===========================================================================

    def _alignFromFASTA(self, FASTA, debug=False):

        """
        Populates the multAlign from a FASTA file. Currently the reference is simply
        not set.
        Note: we will assign a src_start and src_end that simply goes from bp 1 -
        <length of sequences> when parsing from fasta files...

        Args:
            FASTA: string path to the desired file
                    reference

            internal/externalGapChar: see init description.

        TODO:
            -Is there any additional metadata we could grab?
            -Is the orientation guarenteed +? (fix if not!!!!)


        """

        nameCapture = r"^>(?P<name>\S+)"  # captures '>[string here]' at line start ONLY
        seqCapture = r"[\w \-]+"  # captures sequence data [variou]
        leadingWhite = r"^\s+\w+"
        trailingWhite = r"w+\s+$"

        seqMetaData = []
        seqNameCount = {}
        seqList = []
        currentIndex = -1  # this way we just add 1 every name capture and start at 0
        firstSeqLine = True
        firstLineEver = True

        with open(FASTA, "r") as FA:
            for strLine in FA:
                if re.match(nameCapture, strLine):
                    seqName = re.match(nameCapture, strLine).groupdict()["name"]

                    if seqName not in seqNameCount:
                        seqNameCount["seqName"] = 1
                    else:
                        warnings.warn(
                            "\nMultALign Warning: In fasta file import, "
                            + "the seq id "
                            + seqName
                            + "is duplicated- ints (starting "
                            + "at 1) will be appended to name of subsequent appearances\n"
                        )
                        seqName = seqName + str(seqNameCount)
                        seqNameCount["seqName"] = seqNameCount["seqName"] + 1

                    metaDict = {
                        "left_flanking": None,
                        "right_flanking": None,
                        "name": seqName,
                        "src_start": 1,
                        "src_end": None,
                        "src_orientation": "+",
                        "src_transitions": None,
                        "src_transversions": None,
                        "src_gc_background": None,
                        "aligned_start": None,
                        "aligned_end": None,
                    }
                    seqMetaData.append(metaDict)

                    if not firstLineEver:
                        seq = seqList[-1]
                        length = len(seq)
                        seq = seq.replace("-", self.internalGapChar)
                        seq = seq.lstrip(self.internalGapChar).rjust(
                            length, self.externalGapChar
                        )
                        seq = seq.rstrip(self.internalGapChar).ljust(
                            length, self.externalGapChar
                        )
                        seqList[-1] = list(seq)
                        if length != len(seqList[0]):
                            seqNum = len(seqList) - 1
                            raise ValueError(
                                "MultAlign _alignFromFASTA() Error: sequence #"
                                + seqNum
                                + ", named "
                                + seqMetaData[i]["name"]
                                + ", was not"
                                + " the same length as seq 0"
                            )

                    currentIndex += 1
                    seqList.append("")
                    firstSeqLine = True
                    firstLineEver = False

                elif re.match(seqCapture, strLine):
                    seqList[currentIndex] += strLine.rstrip("\n")

        seq = seqList[-1]
        length = len(seq)
        seq = seq.replace("-", self.internalGapChar)
        seq = seq.rstrip(self.internalGapChar).ljust(length, self.externalGapChar)
        seq = seq.lstrip(self.internalGapChar).rjust(length, self.externalGapChar)
        seqList[-1] = list(seq)

        for metaDict in seqMetaData:  # This prevents issues with printing
            metaDict["src_end"] = len(seqList[0])

        self.aligned_sequences = np.asarray(seqList)
        self.seqMetaData = seqMetaData

        self.reference = np.asarray(list(self.consensus()))
        self.refMetaData = {
            "left_flanking": None,
            "right_flanking": None,
            "name": "consensus",
            "src_start": 1,
            "src_end": len(self.reference),
            "src_orientation": "+",
            "src_transitions": None,
            "src_transversions": None,
            "src_gc_background": None,
            "aligned_start": None,
            "aligned_end": None,
        }

    # ===========================================================================
    def toFASTA(
        self,
        fileName=None,
        lineWidth=50,
        addFlankSeqs=False,
        removeExternalGaps=False,
        removeAllGaps=False,
    ):
        """
        Function to create a FASTA file of the MSA data stored in the
        MultAlign object.By default the function will just return a string, but
        options to create a file are included as well.

        Args:
            fileName: if no arg given, function will output a string. Otherwise,
                        accepts a string that will be the name of the new file
                        created

            lineWidth: width of the sequence lines in the output

            addFlankSeqs: accepts boolean value - determines if flanking sequences
                            for aligned sequences should be output as well
                            Default = False
                            NOTE: in order to maintain alignment, additional
                            external gap characters will be appended / prepended
                            if all flanking sequence lengths are not the same.

            removeExternalGaps: accepts a boolean - determines if external gaps
                                should be removed from the output, compromising the
                                aligned nature but reducing the file size
                                Default = False

            removeAllGaps: accepts a boolean value: determines if all gaps, both
                            internal and external, should be omitted from the output
                            file, leaving only base pair data
                            Default = False


        TODO:
            -See from FASTA TODO for not-yet-implemented considerations
            -Can flanking sequences have gaps? If so, must alter the addFlanking
                    section to remove them if user desires...

        """

        # Check to see if user wants output to file - if so, pick the name
        toFile = True
        if fileName == None:
            toFile = False
            if self.refMetaData is not None:
                if self.refMetaData["name"] != None:
                    fileName = self.refMetaData["name"] + "_fromMA.fa"
            elif (
                self.seqMetaData is not None
            ):  # necessary bc ref metadata may not be present for FA, tbd
                fileName = self.seqMetaData[0]["name"] + "_fromMA.fa"
            else:
                fileName = "NoFileNameGivenOrImplied"

        if type(fileName) != str:
            raise TypeError(
                "MultALign Error: 'filename' arg to toFasta" + " must be of type str"
            )

        # a bit of preprocessing for the flanks...
        if addFlankSeqs:
            flanksExist = False
            flanksHaveValue = False
            largestLFlank = 0
            largestRFlank = 0
            for metaSet in self.seqMetaData:
                if "left_flanking" in metaSet:
                    flanksExist = True
                    if metaSet["left_flanking"] != None:
                        flanksHaveValue = True
                        if len(metaSet["left_flanking"]) > largestLFlank:
                            largestLFlank = len(metaSet["left_flanking"])

            if "right_flanking" in metaSet:
                flanksExist = True
                if metaSet["right_flanking"] != None:
                    flanksHaveValue = True
                    if len(metaSet["right_flanking"]) > largestRFlank:
                        largestRFlank = len(metaSet["right_flanking"])

            if flanksExist == False:
                addFlankSeqs = False
                warnings.warn(
                    "MultAlign Warning: toFasta(): 'addFlankSeqs' arg was "
                    + "set to True, but seq meta dicts contain no "
                    + "'left_flanking' or 'right_flanking' keys\n"
                )
            elif flanksHaveValue == False:
                addFlankSeqs = False
                warnings.warn(
                    "MultAlign Warning: toFasta(): 'addFlankSeqs' arg was "
                    + "set to True, but all sequence meta dicts contain "
                    + "None as their flank values"
                )

            elif largestLflank == 0 and largestRFlank == 0:
                addFlankSeqs = False
                warnings.warn(
                    "MultAlign Warning: toFasta(): 'addFlankSeqs' arg was"
                    + " set to True, but all flanking seqs have length 0"
                )

        FAOutput = ""
        length = len(self.aligned_sequences[0])

        for i, seq in enumerate(self.aligned_sequences):
            seqToOutput = "".join(seq)
            startIndex = 0
            endIndex = lineWidth

            # get rid of gaps if they please...
            if removeAllGaps:
                seqToOutput = seqToOutput.replace(internalGapChar, "").replace(
                    externalGapChar, ""
                )

            elif removeExternalGaps:
                seqToOutput = seqToOutput.lstrip(externalGapChar).rstrip(
                    externalGapChar
                )

            # Add flanks if they please...
            if addFlankSeqs:
                seqToOutput = (
                    self.seqMetaData[i]["left_flanking"]
                    + seqToOutput
                    + self.seqMetaData[i]["right_flanking"]
                )

                # In case of variable length flanking seqs, we must add additional
                # gaps to maintain alignment, unless they don't care about that...
                if removeAlignmentGaps == False and removeAllGaps == False:
                    lDif = largestLFlank - len(self.seqMetaData[i]["left_flanking"])
                    rDif = largestRFlank - len(self.seqMetaData[i]["right_flanking"])
                    seqToOutput = ("-" * lDif) + seqToOutput + ("-" * rDif)

            FAOutput += ">"
            FAOutput += self.seqMetaData[i]["name"] + "\n"

            while endIndex < length:
                FAOutput += (
                    seqToOutput[startIndex:endIndex].replace(" ", "-").replace(".", "-")
                    + "\n"
                )
                startIndex += lineWidth
                endIndex += lineWidth

            FAOutput += seqToOutput[startIndex:].replace(" ", "-").replace(".", "-")
            if i != len(self.aligned_sequences):
                FAOutput += "\n"

        # TODO: after you test, add the file options
        if toFile:
            outputFile = open(fileName, "w")
            outputFile.write(FAOutput)
            outputFile.close()
        return FAOutput

    # ===========================================================================

    def _alignFromStockholm(self, stockholm, debug=False, desiredIndex=0):

        """
        Helper function to populate MultAlign object with data from a given stockholm
        file.

        Args:
            stockholm: accepts string path to the desired stockholm file

            internal/externalGapChar: see init

            desiredIndex: accepts an integer, default 0. Stockholm files may
                            contain data for multiple MSAs, but a single multAlign
                            object can only hold 1. In the event that the desired
                            stockholm contains multiple MSAs, this arg allows the
                            user to select which one should be parsed into the
                            MultAlign object.

            For users wanting to import multiple MSAs from a Stockholm file into
            several MultAligns, see function of the same name in MultAlignCollection.py



        TODO:
        NOTE: the current implementation only allows the user to select 1 desired
            index from the file and create a multalign from that multiple alignment:
            to get all alignments, the function must be called from
            MultAlignCollection to prevent circular imports - tbd if we keep that
            If you would like to translate all MSAs in a stockholm file to MultAligns,
            please use the function of similar name in MultAlignCollection.py
        -Implement internal and external gap char args
        -trim fat (especially if we keep the 1 alignment rule)
        -Should we use the stk 'reference' as the reference for a multAlign, or
            another sequence (how do we pick?)?
        -Don't worry about additional metadata - wait for meeting with Jeb and Rob
        """
        # Smitten V1  e.g. seq1:1-100 or seq1:203-20 (reverse strand)
        smittenV1 = re.compile("^(\d+)-(\d+)$")
        # Smitten V2  e.g. seq1:1-100_+ or seq1:1-100_- (reverse strand)
        smittenV2 = re.compile("^(\d+)-(\d+)_([+-])$")

        if not os.path.isfile(stockholm):
            raise ValueError(
                "Error: Stockholm argument to _alignFromStockholm() is not a valid"
                + " File path! Path given: \n"
                + str(stockholm)
                + "\n"
            )

        # Code below is based on implementation in DfamSeedAlignment.py - fromStockholmFile()
        # Every record starts with a "STOCKHOLM" version tag
        stkHeaderRE = re.compile("^#\s+STOCKHOLM\s+[\d\.]+")
        seqPrefixRE = re.compile("^([\-\.]+).*")
        seqSuffixRE = re.compile(".*[^\.\-]([\-\.]+)$")
        idNomenclatureRE = re.compile(
            "^(?:([^:\s]+):)?(?:([^:\s]+):)?(\S+):(\d+)-(\d+)$"
        )
        fieldCodeDataRE = re.compile(r"#=GF\s*(\S*)\s*(.*)")
        alignDataRE = re.compile(
            "^(\S+)\s+([\.\-ACGTRYKMSWBDHVNXacgtrykmswbdhvnx]+)\s*$"
        )
        metadata = {}
        alignments = []
        uniq_ident = {}
        reference = None
        interpret_ids = True
        currentIndex = 0

        def set_field_once(field, value):
            nonlocal metadata
            if metadata[field] is None:
                metadata[field] = value
            else:
                raise StockholmParseError(
                    "Only one value is allowed for the '{}' field".format(field)
                )

        with open(stockholm, "r") as stk:
            for strLine in stk:
                # print("line: " + strLine)
                if strLine.startswith("//"):
                    uniq_ident = {}
                    if alignments:
                        # Adjust model_start and model_end for alignments
                        if reference:
                            for align in alignments:
                                seq = align[3]
                                if len(seq) != len(reference):
                                    raise Exception(
                                        "Reference and sequence lengths differ! "
                                        + metadata["name"]
                                        + " ["
                                        + str(len(seq))
                                        + "] and reference ["
                                        + str(len(reference))
                                        + "]"
                                    )

                                model_start = None
                                model_end = None

                                i = 0
                                for (r, s) in zip(reference, seq):
                                    if r not in ["-", "."]:
                                        i += 1

                                        if s not in ["-", " "]:
                                            # shared position, update model_end
                                            model_end = i
                                            # update model_start but only for the first shared position
                                            if model_start is None:
                                                model_start = i

                                if model_end is None:
                                    # TODO: handle this case (sequence shares no columns with reference)
                                    print(
                                        "!!! Sequence "
                                        + metadata["name"]
                                        + "shares no columns with reference!!"
                                    )
                                    print("ref:" + reference)
                                    print("seq:" + seq)
                                    pass

                                align[5] = model_start
                                align[6] = model_end

                        """
                        print("\n======== I N Q U I S I T I O N =========================\n")
                        print("type of alignments: " + str(type(alignments)))
                        print("len of alignments: " + str(len(alignments)))
                        print("1st alignments: \n" + str(alignments[0]))

                        print("\nThe reference is: " + reference)
                        print("\nThe name is: " + metadata['name'])
                        print("\n======== E N D   I N Q U I S I T I O N =========================\n")
                        """
                        # Recall that in a file of potentially n of these, we only
                        # want this one...
                        if currentIndex == desiredIndex:
                            toBeRef = reference  # Here we replace gaps with those desired by user...
                            length = len(toBeRef)
                            toBeRef = toBeRef.replace(
                                ".", self.internalGapChar
                            ).replace("-", self.internalGapChar)
                            toBeRef = toBeRef.rstrip(self.internalGapChar).ljust(
                                length, self.externalGapChar
                            )
                            toBeRef = toBeRef.lstrip(self.internalGapChar).rjust(
                                length, self.externalGapChar
                            )

                            self.reference = np.asarray(list(reference))
                            self.refMetaData = {
                                "left_flanking": None,
                                "right_flanking": None,
                                "name": metadata["name"],
                                "src_start": None,
                                "src_end": None,
                                "src_orientation": None,
                                "src_transitions": None,
                                "src_transversions": None,
                                "src_gc_background": None,
                            }

                            alignedSequences = []
                            seqMetaDictList = []
                            for sequence in alignments:
                                toBeSeq = sequence[3]
                                toBeSeq = toBeSeq.replace(
                                    " ", self.internalGapChar
                                ).replace("-", self.internalGapChar)
                                toBeSeq = toBeSeq.rstrip(self.internalGapChar).ljust(
                                    length, self.externalGapChar
                                )
                                toBeSeq = toBeSeq.lstrip(self.internalGapChar).rjust(
                                    length, self.externalGapChar
                                )
                                alignedSequences.append(list(toBeSeq))

                                seqMetaDict = {
                                    "left_flanking": None,
                                    "right_flanking": None,
                                    "name": sequence[0],
                                    "src_start": sequence[1],
                                    "src_end": sequence[2],
                                    "src_orientation": sequence[4],
                                    "src_transitions": None,
                                    "src_transversions": None,
                                    "src_gc_background": None,
                                }
                                seqMetaDictList.append(seqMetaDict)

                            self.seqMetaData = seqMetaDictList
                            self.aligned_sequences = np.asarray(alignedSequences)

                        """
                        NOTE: If you are adding additional metadata to the MultAlign
                        in order to capture more of the metadata from stockholm files,
                        check out DfamSeedAlignment.py and look at the construction
                        of the seeds as a shortcut to know what metadata is being
                        captured by the parser and what they are!

                        seed = DfamSeedAlignment( reference = reference,
                                        alignments = alignments,
                                        accession = metadata["accession"],
                                        name = metadata["name"], title = metadata["title"],
                                        description = metadata["description"],
                                        authors = metadata["authors"],
                                        classification = metadata["classification"],
                                        seed_source = metadata["seed_source"],
                                        target_site_cons = metadata["target_site_cons"],
                                        curation_notes = metadata["curation_notes"],
                                        citations = metadata["citations"], clades = metadata["clades"],
                                        aliases = metadata["aliases"])

                        seeds.append(seed)
                        """
                    currentIndex += 1
                    continue

                if stkHeaderRE.match(strLine):
                    # Start of a Stockholm record
                    metadata = {
                        "name": None,
                        "accession": None,
                        "title": None,
                        "authors": None,
                        "classification": None,
                        "seed_source": None,
                        "citations": [],
                        "clades": [],
                        "aliases": [],
                        "description": None,
                        "curation_notes": None,
                        "target_site_cons": None,
                    }
                    alignments = []
                    reference = None
                    continue

                if strLine.startswith("#=GF"):
                    fields = fieldCodeDataRE.match(strLine)
                    if fields:
                        field_key = fields.group(1)
                        field_value = fields.group(2)

                        if field_key == "ID":
                            set_field_once("name", field_value)
                        elif field_key == "AC":
                            set_field_once("accession", field_value)
                        elif field_key == "DE":
                            set_field_once("title", field_value)
                        elif field_key == "AU":
                            set_field_once("authors", field_value)
                        elif field_key == "TP":
                            set_field_once("classification", field_value)
                        elif field_key == "SE":
                            set_field_once("seed_source", field_value)
                        elif field_key == "OC":
                            metadata["clades"] += [field_value]
                        elif field_key == "DR":
                            metadata["aliases"] += [field_value]
                        elif field_key == "RN":
                            metadata["citations"] += [{"title": None, "authors": None}]
                        elif field_key == "RM":
                            metadata["citations"][-1]["pmid"] = field_value
                        # Title allowed to span multiple lines
                        elif field_key == "RT":
                            if metadata["citations"][-1]["title"] is None:
                                metadata["citations"][-1]["title"] = ""
                            else:
                                metadata["citations"][-1]["title"] += "\n"
                            metadata["citations"][-1]["title"] += field_value
                        # Authors allowed to span multiple lines
                        elif field_key == "RA":
                            if metadata["citations"][-1]["authors"] is None:
                                metadata["citations"][-1]["authors"] = ""
                            else:
                                metadata["citations"][-1]["authors"] += "\n"
                            metadata["citations"][-1]["authors"] += field_value
                        elif field_key == "RL":
                            metadata["citations"][-1]["journal"] = field_value
                        elif field_key == "TD":
                            set_field_once("target_site_cons", field_value)
                        elif field_key == "CC":
                            if metadata["description"] is None:
                                metadata["description"] = ""
                            else:
                                metadata["description"] += "\n"

                            metadata["description"] += field_value
                        elif field_key == "**":
                            if metadata["curation_notes"] is None:
                                metadata["curation_notes"] = ""
                            else:
                                metadata["curation_notes"] += "\n"

                            metadata["curation_notes"] += field_value
                    continue

                if strLine.startswith("#=GC"):
                    # Consensus Line
                    tokens = strLine.split()
                    if len(tokens) >= 3:
                        reference = tokens[2]
                    continue

                mats = alignDataRE.match(strLine)
                if mats is not None:
                    full_id = mats.group(1)
                    align_seq = mats.group(2).upper()

                    # Warn about duplicate identifiers
                    if full_id in uniq_ident:
                        # raise StockholmParseError("Non unique identifier: " + full_id)
                        print(
                            "Warning: skipping duplicate sequence " + full_id,
                            file=sys.stderr,
                        )
                        continue
                    uniq_ident[full_id] = 1

                    tmp = align_seq
                    tmp = tmp.replace("-", "")
                    tmp = tmp.replace(".", "")
                    seq_len = len(tmp)

                    align_seq = align_seq.replace(".", "-")

                    pre = seqPrefixRE.match(align_seq)
                    if pre is not None:
                        align_seq = (
                            self.externalGapChar * len(pre.group(1))
                            + align_seq[len(pre.group(1)) :]
                        )
                    suf = seqSuffixRE.match(align_seq)
                    if suf is not None:
                        align_seq = align_seq[
                            : -len(suf.group(1))
                        ] + self.externalGapChar * len(suf.group(1))

                    # Handle Smitten V1/V2 sequence identifiers (without normalization)
                    # otherwise treat the entire string as the sequence id
                    seq_id = full_id
                    seq_start = None
                    seq_end = None
                    strand = None
                    if interpret_ids:
                       idflds = full_id.split(":")
                       #  Must have at least one colon otherwise we cannot interpret it
                       if len(idflds) > 1:
                           match = smittenV2.match(idflds[-1])
                           if match:
                               seq_id = ":".join(idflds[:-1])
                               seq_start = int(match.group(1))
                               seq_end = int(match.group(2))
                               strand = match.group(3)
                           else:
                               match = smittenV1.match(idflds[-1])
                               if match:
                                   seq_id = ":".join(idflds[:-1])
                                   seq_start = int(match.group(1))
                                   seq_end = int(match.group(2))
                                   if seq_start <= seq_end:
                                       strand = "+"
                                   else:
                                       tmp = seq_start
                                       seq_start = seq_end
                                       seq_end = tmp
                                       strand = "-"

                    # Note: these positions are adjusted afterwards once we have a reference
                    model_start = 0
                    model_end = 0

                    #print("seq_id: " + seq_id)
                    alignments.append(
                        [
                            seq_id,
                            seq_start,
                            seq_end,
                            align_seq,
                            strand,
                            model_start,
                            model_end,
                        ]
                    )

        return None

    # ===========================================================================
    def toStockholm(self, fileName=None, refX=True, version="1.0", description=None):
        """
        Function to create a stockholm file of the MSA data stored in the
        MultAlign object.By default the function will just return a string, but
        options to create a file are included as well. By default, the reference
        will be converted to the consensus and all non gap chars changed to x's,
        though this can be prevented by changing refX to False

        Description is optional and will create the GC DE line in the stk
        Including toFile but excluding fileName will auto name file
        '[referenceName]_fromMA.stk'
        If the filename given exists, the data stk will be appended. If not,
        the file will be created


        TODO:
        -Lots of stk metadata is not stored in current MultAlign - fix?
            -see todos below
        -implement additional argument functionality
            -Reerence Selector
            -Optionally store metadata NOT found in stockholm in the
            -Directory
        -find out what metadata is compulsory
        -Should we add ability to append to a pre-exiting file (since they can
         contain more than 1)?
        """
        # Check to see if user wants output to file - if so, pick the name
        # TODO: cleanup after quick fix about toFIle
        toFile = False
        if fileName != None:
            toFile = True

        if toFile:
            if fileName == None:
                fileName = self.refMetaData["name"] + "_fromMA.stk"

            if type(fileName) != str:
                raise TypeError(
                    "MultALign Error: 'filename' arg to toStockholm"
                    + " must be of type str"
                )
        stkOutput = ""
        versionLine = "# STOCKHOLM " + str(version) + "\n"
        name = re.sub(r"_fromSTK", "", self.refMetaData["name"])
        idLine = "#=GF ID    " + name + "\n"

        defLine = ""
        if description != None:
            defLine = "#=GF DE    " + description + "\n"  # TODO: implement metadata
        numLine = "#=GF SQ    " + str(len(self.aligned_sequences)) + "\n"

        # TODO: Because we can't get a raw string from the gap characters, things could
        # break if someone were to use some unsavory gap characters like backslash,
        # but perhaps we will just let people shoot themselves in foot in this way
        # if they want....
        charsToReplace = (
            "[" + self.internalGapChar + self.externalGapChar + "]"
        )  # produces '[xy]', list of chars to substitte below
        ref = re.sub(
            charsToReplace, r".", "".join(self.reference)
        )  # replace gap chars with '.'
        if refX:
            ref = re.sub(r"[a-zA-Z]", r"x", ref)
        refLine = "#=GC RF    " + ref + "\n"

        stkOutput += versionLine + idLine + defLine + numLine + refLine

        for i, seq in enumerate(self.aligned_sequences):
            metaInfo = self.seqMetaData[i]
            seqLine = (
                metaInfo["name"]
                + ":"
                + str(metaInfo["src_start"])
                + "-"
                + str(metaInfo["src_end"])
                + "    "
            )
            seqLine += re.sub(r"[ -]", r".", "".join(seq)) + "\n"
            stkOutput += seqLine

        stkOutput += "//"

        # TODO: after you test, add the file options
        if toFile:
            outputFile = open(fileName, "a")
            outputFile.write(stkOutput)
            outputFile.close()
        return stkOutput

    # ===========================================================================

    def _alignFromDNAPAC(
        self,
        DNAPAC,
        reference=None,
        alignment_reference=None,
        flankingSequenceDatabase=None,
        maxFlankingSequenceLen=50,
        debug=False,
    ):
        """
        DNAPAC = DNAPairwiseAlignmentCollection
        This is a 'private' method for generating a multAlign from
        a search result collection containing a search of one sequence
        against many others. Using this function on MSAs NOT created in this manner
        will most likely result in a crash

        Args:
            reference: accepts string, 'reference' or 'query'. If none given,
                        default is 'query'. Should match 'referece' arg of the
                        DNAPAC in question

            alignmentReference: accepts same args as 'reference' arg, FIXME:
                                what is actually supposed to be happening with this?

            maxFlankingSequenceLen: accepts int that determins the maximum number
                                    of flanking characters to store

            internal/externalGapChar: see init

        TODO:
            -- default max flank of 50 sound good?
            -- figure out what to expect from 'flanking seq databases'
                    -"Don't worry about it now" -Robert

            --FIXME! what is supposed to be done with reference arg?

            --Do we want to make it so that changing the 'internalGapChar' argument
              changes what the funciton looks for as a gap in reconstructing the
              reference, or should it just always be '-' expected as a rule?

        """

        # Firstly, we must ensure we are recieving valid input...
        if not isinstance(DNAPAC, DNAPairwiseAlignmentCollection):
            raise ValueError(
                "Error: DNAPAC must be of type DNAPairwiseAlignmentCollection: "
                + str(type(DNAPAC))
                + " was given instead!"
            )

        if DNAPAC.len() < 1:
            raise ValueError("Error: DNAPAC must have at least 1 alignment!")

        """
        Similar to the original, we will define functions for use here that will
        allow us to change which member data we access based on what the user dictates
        as the reference sequence.

        #TODO: somethings up with orientation
        #Could do 'smart' tight loop at top to determine 'common' one and see
        #which is query / target
        """
        if alignment_reference == None or alignment_reference == "query":

            def getRefStart(align):
                return align.query_start

            def getRefEnd(align):
                return align.query_end

            def getRefLen(align):
                return align.query_len

            def getRefSeq(align):
                return align.aligned_query_seq

            def getRefName(align):
                return align.query_name

            def getAlignedSeqStart(align):
                return align.target_start

            def getAlignedSeqEnd(align):
                return align.target_end

            def getAlignedSeqLen(align):
                return align.target_len

            def getAlignedSeqSeq(align):
                return align.aligned_target_seq

            def getAlignedSeqName(align):
                return align.target_name

            def getAlignedSeqOrientation(align):
                return align.orientation

        elif alignment_reference == "target":

            def getRefStart(align):
                return align.target_start

            def getRefEnd(align):
                return align.target_end

            def getRefLen(align):
                return align.target_len

            def getRefSeq(align):
                if align.orientation != "+":
                    return self.reverseComplement(align.aligned_target_seq)
                return align.aligned_target_seq

            def getRefName(align):
                return align.target_name

            def getRefOrientation(align):
                return align.orientation

            def getAlignedSeqStart(align):
                return align.query_start

            def getAlignedSeqEnd(align):
                return align.query_end

            def getAlignedSeqLen(align):
                return align.query_len

            def getAlignedSeqSeq(align):
                if align.orientation != "+":
                    return self.reverseComplement(align.aligned_query_seq)
                return align.aligned_query_seq

            def getAlignedSeqName(align):
                return align.query_name

            def getAlignedSeqOrientation(align):
                return align.orientation

            # In this case, we will have to do a bit of reverse complementing...

        else:
            raise ValueError(
                "Error: 'alignment_reference' Argument to _alignFromDNAPAC() MUST "
                + "be either 'query' or 'target' - "
                + str(alignment_reference)
                + " was given..."
            )

        # OLDFIXME: here in the original multalign there is an 'if orientation matters'
        # statement that reverse complements some sequences, but the comment states
        # that 'Since we are only concerned about character
        #    matching we do not need to correct the sequence before
        #    proceding....but we do save the orientation in case anyone
        #    else wants to know.'
        # What is the truth?     -Robert said not to worry about it

        """
        Here, we must reconstruct the reference sequence without any gaps, as this
        is representative of the original sequence before alignment, and gaps only
        occur in the pieces of the reference where insertions occurred in the
        sequence aligned to the reference. The positions of gaps in the reconstructed
        reference will be computed below.
        """

        combRefSeq = ""  # This will be our reference
        tRemaining = 0
        tMin = 0  # This will be the lowest coordinate of our reconstructed ref
        tMax = 0  # this will be the highest coordinate of the reconstructed ref
        gapChars = r"[\s\-.]"

        if reference == None or reference == "":
            tMin = getRefStart(DNAPAC.get(0))
            tMax = getRefEnd(DNAPAC.get(0))

            # find our tmin and tmax
            for align in DNAPAC:
                if getRefStart(align) < tMin:
                    tMin = getRefStart(align)
                if getRefEnd(align) > tMax:
                    tMax = getRefEnd(align)
                    tRemaining = getRefLen(align) - getRefEnd(align)

            combRefSeq = list("".ljust(tMax - tMin + 1))
            # print("\nLen of combRefSeq (at top): " + str(len(combRefSeq))) Not this...

            # loop through and reconstruct the reference a chunk at a time based
            # on coordinates, and removing "-" gap characters
            for align in DNAPAC:
                startIndex = getRefStart(align) - tMin
                basesToInsert = getRefSeq(align).replace("-", "")
                endIndex = startIndex + len(basesToInsert)
                combRefSeq[startIndex:endIndex] = list(basesToInsert)

            # We throw an error if we don't get 'full coverage', which is to say we
            # cannot reconstruct an entire contiguous chunk of the reference.
            # A reconstructed ref of "EFGHIJKL" would be allowable as it's 1
            # contiguous chunk, but "ABCDEF     MNOPQRS" would not.
            if re.search(gapChars, "".join(combRefSeq)):
                raise ValueError(
                    "Error: Reference sequence is incomplete, no"
                    + "coverage availalbe for at least one region."
                )

        else:  # If the user gives us the reference, we can skip all that above
            combRefSeq = list(reference)
            tMin = 1
            tMax = len(combRefSeq)

        if debug:
            print(
                "\nAlignmentFromDNAPAC: Reconstructed Ref Seq is: "
                + "".join(combRefSeq)
                + "\n"
            )

        """
        Here we create 'gap patterns' for each alignment. The patterns will be used
        to determine where gaps should be placed in the reconstructed reference,
        as well as in the aligned sequences to maintain their alignment .If two
        sequences are aligned to the same region of the reference, but have
        differing insertions (aka their references have differing gaps,) we will
        use this data to add gaps to the sequences as well, e.g.:

        ref1    ABCDEFGHI-JKL   ref2    HI---J-KLMNOP   ref3    MNOPQRSTUVWXYZ
        seq1    ABCDEFGHIxJKL   seq2    HIxxxJxKLMNOP   seq3    MNOPQRSTUVWXYZ

        which /should/ come together to form a multAlign like this:
        ______________________________________
        Refe    ABCDEFGHI---J-KLMNOPQRSTUVWXYZ      Note: for alphabet examples,
        seq1    ABCDEFGHIx--J-KL                    the index of each 'base' is
        seq2           HIxxxJxKLMNOP                just its postition in the
        seq3                    MNOPQRSTUVWXYZ      alphabet (A = 1, B=2, ...)
        ______________________________________
        Notice how additional gaps must be 'added' to seq1, and the larger
        gaps found in ref2 reign supreme over the gaps at that position in the
        reconstructed reference
        """

        gapPatterns = [None] * (DNAPAC.len())
        index = 0
        notGaps = r"[^-]"
        seqMetaDictList = []
        for align in DNAPAC:

            # Here we deal with ref gaps
            unGaps = re.split(notGaps, getRefSeq(align))
            # this makes a list of gaps in the sequence, e.g.: ['', '', '---', ...]
            startDif = getRefStart(align) - tMin
            prependedZeros = [0] * startDif
            gapPatterns[index] = prependedZeros + list(map(len, unGaps))
            # The above line may be a bit confusing... the prepended zeros assure
            # that all gap patterns are aligned to eachother and to the reconstructed
            # reference. list(map(len, unGaps)) is making a list of the mapping of
            # the length (hence len() function as an argument) of the unGaps list,
            # resulting in something like [0, 0, 3, ....]. It can be seen as a list
            # of the length of gaps that come before and after a given base in the
            # sequence

            # Here we grab sequence metadata:
            # TODO: tranistions and transversions, or no?
            seqMetaDict = {
                "left_flanking": None,
                "right_flanking": None,
                "name": getAlignedSeqName(align),
                "src_start": getAlignedSeqStart(align),
                "src_end": getAlignedSeqEnd(align),
                "src_orientation": align.orientation,
                "src_transitions": None,
                "src_transversions": None,
                "src_gc_background": None,
            }
            seqMetaDictList.append(seqMetaDict)
            index += 1

        # Deal with some metadata...
        # TODO: implement reference transition and transversion collection
        self.seqMetaData = seqMetaDictList
        self.refMetaData = {
            "left_flanking": None,
            "right_flanking": None,
            "name": getRefName(DNAPAC[0]),
            "src_start": tMin,
            "src_end": tMax,
            "src_orientation": "+",
            "src_transitions": None,
            "src_transversions": None,
            "src_gc_background": None,
        }

        """
        Now we must find the maximum of the gap patterns. As seen in the example
        above, if two alignments covering the same portion of the reference have
        gaps of different sizes, the reconstructed reference will have a gap of the
        larger size at that position.

        """
        maxPattern = [0] * (len(combRefSeq) + 1)
        maxPatternOriginNum = [None] * (len(combRefSeq) + 1)

        index = 0
        alignNum = 0
        for pattern in gapPatterns:
            # print(pattern)
            index = 0
            for count in pattern:
                if count > maxPattern[index]:
                    maxPattern[index] = count
                    maxPatternOriginNum[index] = alignNum
                index += 1
            alignNum += 1

        if debug:
            print("_alignFromDNAPAC: max gap pattern: ")
            print(maxPattern)

        """
        Now that we have the maxPattern, we can simply iterate over the pattern
        and reconstructed refernce to create the final version of the reference,
        which will have gaps in the appropriate places
        """
        finalRef = ""
        for i in range(len(maxPattern) - 1):
            finalRef += "".ljust(maxPattern[i], self.internalGapChar)
            finalRef += combRefSeq[i]

        if debug:
            print("_alignFromDNAPAC: reconstructed ref with gaps: ")
            print(finalRef)

        """
        Now (as also seen in the example above) we must add additional gaps to the
        aligned sequences so that each one maintains its alignment to the reference
        now that the reconstructed reference has a myriad of gaps inserted in positions
        where other alignments that cover the same region of the reference may not
        have had such gaps, or such long gaps.
        """
        # first computing 'totalGaps' allows us to compute the total number of
        # spaces before a given index in the reference, and thus compute the
        # number of spaces that need to be prepended to a sequence to align it
        # to the reconstructed ref
        totalGaps = [0]
        for i in range(len(maxPattern) - 1):
            totalGaps.append(totalGaps[i] + maxPattern[i])
        """
        Then for each alignment, we must go through and add the appropriate
        gaps. Notice how gaps will not be added for a sequence if the portion of
        the reference it aligns to has a '-' at that position, as well as how when
        we do insert gaps, we insert the number of gaps called for by the maxPattern
        /MINUS/ the number of gaps that occur in that alignments specific gap pattern
        to prevent redundently adding gaps and resulting in failure to align
        """
        finalOldFashionedSeqs = []
        for l in range(len(DNAPAC)):
            currentAlign = DNAPAC[l]
            start = getRefStart(currentAlign) - tMin
            start += totalGaps[start]
            seq = ""
            k = getRefStart(currentAlign) - tMin
            for i in range(len(getRefSeq(currentAlign))):
                if getRefSeq(currentAlign)[i] != "-":
                    numGaps = maxPattern[k]
                    if k < len(gapPatterns[l]):
                        numGaps -= gapPatterns[l][k]
                    seq += "".ljust(numGaps, self.internalGapChar)
                    k += 1
                seq += getAlignedSeqSeq(currentAlign)[i]
            seq = "".ljust(start, " ") + seq
            seq = seq.ljust(len(finalRef), " ")
            length = len(seq)
            seq = seq.replace(" ", self.internalGapChar).replace(
                "-", self.internalGapChar
            )
            seq = seq.rstrip(self.internalGapChar).ljust(length, self.externalGapChar)
            seq = seq.lstrip(self.internalGapChar).rjust(length, self.externalGapChar)
            finalOldFashionedSeqs.append(list(seq))

            if len(finalOldFashionedSeqs[-1]) != len(finalOldFashionedSeqs[0]):
                raise ValueError(
                    "MultAlign _alignFromDNAPAC() Error: sequence #"
                    + l
                    + ", named "
                    + seqMetaDictList[l]["name"]
                    + ", was not"
                    + " the same length as seq 0"
                )

        self.reference = np.asarray(list(finalRef))
        self.aligned_sequences = np.asarray(finalOldFashionedSeqs)

        if len(self.reference) != len(self.aligned_sequences[0]):
            raise ValueError(
                "MultAlign _alignFromDNAPAC() Error: reference"
                + " was not the same length as seq 0"
            )

        if debug:
            print("\n")
            print("|".join(list(finalRef)))
            print("-".ljust(len(finalRef) * 2, "-"))
            for alignedSeq in finalOldFashionedSeqs:
                print("|".join(alignedSeq))

    # ===========================================================================

    def reverseComplementMultAlign(self):
        """
        Reverse Complements the MultAlign by doing so to each sequence
        (Transcribing base pairs, reversing sequences, swapping start and end)
        #TODO: find out what other metaData values need switched
        #Fix to allow for lacking metadata - have to start remembering none is required!
        """
        for i, seqMeta in enumerate(self.seqMetaData):
            self.seqMetaData[i]["src_start"], self.seqMetaData[i]["src_end"] = (
                seqMeta[i]["src_end"],
                seqMeta[i]["src_start"],
            )  # swap without holder
            self.aligned_sequences[i] = np.asarray(
                list(self.reverseComplement(self.aligned_sequences[i]))
            )
            if self.seqMetaData[i]["src_orientation"] is not None:
                if self.seqMetaData[i]["src_orientation"] == "+":
                    self.seqMetaData[i]["src_orientation"] = "-"
                else:
                    self.seqMetaData[i]["src_orientation"] = "+"

        self.refMetaData[i]["src_start"], self.refMetaData[i]["src_end"] = (
            self.refMetaData[i]["src_end"],
            self.refMetaData[i]["src_start"],
        )  # swap without holder
        self.reference[i] = np.asarray(list(self.reverseComplement(self.reference[i])))

    # ===========================================================================

    def illegalCharacterChecker(self):
        """
        Function to ensure that no 'illegal' characters are included in sequences or
        references - can be called on object construction with optional arg OR any
        time after

        Currently used to force users to only use accepted gaps ('.', ' ', '-') and
        IUB codes upon construction, unless arg checkForIllegalChars is set to False
        by user
        """
        allowedCharacters = "ACGTRYKMSWBDHVNXacgtrykmswbdhvnx-. "
        counter = 0
        if self.aligned_sequences is not None:
            for x in self.aligned_sequences.flat:
                if x not in allowedCharacters:
                    seqNum = str(int(counter / len(self.aligned_sequences[0])))
                    bpLocation = str(counter % len(self.aligned_sequences[0]))
                    raise ValueError(
                        "Error: illegal character in Sequence "
                        + seqNum
                        + " at bp #"
                        + bpLocation
                        + ". Char is: "
                        + str(x)
                    )
                counter += 1

        counter = 0
        if self.reference is not None:
            for x in self.reference:
                if x not in allowedCharacters:
                    raise ValueError(
                        "Error: illegal character in Reference at index " + str(counter)
                    )
                counter += 1

    # ===========================================================================

    def trim_off(self, left=0, right=0):
        """
        Function to trim base pairs off the left or right of the alignment - trimming
        can occur on one or both sides by the desired number of base pairs. This
        function directly alters the data structure, and may even change the number
        of sequences if a sequence is trimmed down to nothing but gaps, as such
        sequences will be removed. Errors will be thrown if 1. all sequences are
        removed for being trimmed to nothing but gaps or 2. the reference is trimmed
        to nothing but gaps.

        Arguments:
            left: number of bp to trim off left side of alignment
            right: number of bp to trim off right side of alignment

        TODO: Should we implement safety measures to ensure that people do not
        trim the thing into non-existence (currently just throw error)?
        TODO: Test seq removal feature
        """
        rightCoord = len(self.reference) - right
        toKeep = list(range(left, rightCoord))
        self.aligned_sequences = self.aligned_sequences[:, toKeep]
        self.reference = self.reference[left:rightCoord]

        # Remove sequences trimmed to just gaps...
        indicesToRemove = []
        length = len(self.aligned_sequences)  # Iterate backwards to prevent
        for i in range(length - 1, -1, -1):  # indexing issues with deletion
            isAllInts = np.all(self.aligned_sequences[i] == self.internalGapChar)
            isAllExts = np.all(self.aligned_sequences[i] == self.externalGapChar)
            if isAllInts or isAllExts:
                indicesToRemove.append(i)
                self.aligned_sequences = np.delete(self.aligned_seuqences, i, 0)
                self.seqMetaData = np.delete(self.seqMetaData, i)
                print(
                    "\nMultAlign Notice: sequence "
                    + str(i)
                    + " was trimmed "
                    + "from the MultAlign for being trimmed down to nothing but "
                    + "gap characters by trim_off()."
                )

        if len(self.aligned_sequences) == 0:
            raise ValueError(
                "MultAlign Error: All sequences were trimmed by trim_off()"
                + " down to nothing but gap characters and were removed"
            )

        isAllInts = np.all(self.reference == self.internalGapChar)
        isAllExts = np.all(self.reference == self.externalGapChar)
        if isAllInts or isAllExts:
            raise ValueError(
                "MultAlign Error: The reference sequence was trimmed "
                + "by trim_off() down to nothing but gap characters."
            )

    # ===========================================================================

    def getUngappedRef(self, asNumpyArray=False):
        """
        A function that returns the ungapped reference sequence as a string by
        defualt, or optionally as a numpy array.

        The array (shouldn't) be a reference to the original, but rather a
        completely new numpy array

        Argument asNumpyArray: Accepts a boolean value that will determine if the
                                reference is returned as a string or as a
                                numpy array. Default is false (returns string)
        """
        ungappedRef = "".join(self.reference)
        ungappedRef = ungappedRef.replace(" ", "").replace("-", "").replace(".", "")

        if asNumpyArray:
            ungappedRef = np.array(list(ungappedRef))

        return ungappedRef

    # ===========================================================================

    def getUngappedSeq(self, index, asNumpyArray=False):
        """
        A function that returns an ungapped sequence as a string by
        defualt, or optionally as a numpy array.

        The array (shouldn't) be a reference to the original, but rather a
        completely new numpy array

        Argument asNumpyArray: Accepts a boolean value that will determine if the
                                reference is returned as a string or as a
                                numpy array. Default is false (returns string)
        """
        ungappedSeq = "".join(self.aligned_sequences[index])
        ungappedSeq = ungappedSeq.replace(" ", "").replace("-", "").replace(".", "")

        if asNumpyArray:
            ungappedSeq = np.array(list(ungappedSeq))

        return ungappedSeq

    # ===========================================================================

    def getCoverage(self):
        """
        A function that returns the coverage depth of the MSA

        While the double for loop may make you cringe at a glance, it is (I believe)
        an improvement upon the previous design by cutting down on overhead and
        excess data storage (that is captured by profile but not used here) as well
        as doing the same amount of work just in a different order (adding to
        depth as we move through the columns, rather than moving through all
        columns, then adding up all depths after the fact)

        """
        coverage = []
        for i in range(len(self.aligned_sequences[0])):
            depth = 0
            uniq_chars, counts = np.unique(
                self.aligned_sequences[:, i], return_counts=True
            )
            for j, char in enumerate(uniq_chars):
                if char not in (".", " ", "-"):
                    depth += counts[j]
            coverage.append(depth)

        return coverage

    # ===========================================================================

    def getAlignedStart(self, seqIndex):
        """
        Returns where (relative to the 0th position in the gapped sequence) the
        sequence actually begins (as opposed to the gap chars prepended to it
        for alignment purposes). Output varies from 0 to length of ungapped seq -1.

        While no longer necessary in many places due to the new underlying data
        structure, we have decided to keep this function (from the original
        MultAln.pm) for the few places it may still be needed and where it is not
        convenient to compute the data using the aligned_sequences themselves

        Argument seqIndex: accepts and int that determines the index of the sequence
                            you would like the aligned start for
        """
        if (seqIndex > self.length - 1) or (seqIndex < 0):
            raise valueError(
                "Error:getAlignedStart: invalid index provided, index must be between"
                + " 0 and "
                + str(self.length - 1)
            )

        startPos = 0
        for base in self.aligned_sequences[seqIndex]:
            if str(base) in ("-", ".", " "):
                startPos += 1
                continue
            break

        return startPos

    # ===========================================================================

    def getAlignedEnd(self, seqIndex):
        """
        Similar to get aligned start, but with the end. See getAlignedStart()
        above.
        """
        if (seqIndex > self.length - 1) or (seqIndex < 0):
            raise valueError(
                "Error: invalid index provided, index must be between"
                + " 0 and "
                + str(self.length - 1)
            )

        endPos = len(self.reference)
        reverseSeq = self.aligned_sequences[seqIndex][::-1]
        for base in reverseSeq:
            if str(base) in (" ", "-", "."):
                endPos -= 1
                continue
            break

        return endPos

    # ===========================================================================

    def spaceTrimmer(self, front=True, back=True):
        """
        A function that when called alters the contents of the underlying data
        structure by evaluating (optionally) the front and back of the nucleotides
        in the multalign (both in sequences and reference) and deletes any columns
        that hold no bases (only spaces ' '). If the reference or a single sequence
        has a base pair in a column, the column will not be deleted.

        Args:
            front: accepts a boolean that determines whether spaces should be
                    trimmed off the front (0 to nth) position of the multAlign

            back: accepts a boolean that determines whether spaces should be
                    trimmed off the back (nth to length - 1) position of the multAlign
        """
        # colEvals is a variable containing a list of bools - if all values in a
        # column are the same, there will be a true for that column's index in the list
        colEvals = np.all(
            self.aligned_sequences == self.aligned_sequences[0, :], axis=0
        )
        counter = 0
        colsFromFront = 0
        colsFromBack = 0

        # For each base in the reference, if that base is a " ", and the corresponding
        # index in the 0th sequence is " ", and every value in that column is the
        # same, we will add 1 to the # of columns to chop off the front
        if front:
            for base in np.nditer(self.reference):
                if (
                    base == " "
                    and self.aligned_sequences[0][counter] == " "
                    and colEvals[counter]
                ):
                    colsFromFront += 1
                    counter += 1
                else:
                    break

        if back:
            counter = -1
            # Worry not, this || linear time op is only creating a view of the
            # original array  V
            reversedRef = self.reference[::-1]
            # FIXME: For some reason using np.nditer here will un-reverse the
            # bases and screw everything up - how do we optimize this loop?
            # if you do do the same to getAlignedEnd()
            for base in reversedRef:
                if (
                    base == " "
                    and self.aligned_sequences[0][counter] == " "
                    and colEvals[counter]
                ):
                    colsFromBack += 1
                    counter -= 1
                else:
                    break

        # See trim_off(), it's the same code but no func call to reduce overhead
        rightCoord = len(self.reference) - colsFromBack
        toKeep = list(range(colsFromFront, rightCoord))
        self.aligned_sequences = self.aligned_sequences[:, toKeep]
        self.reference = self.reference[colsFromFront:rightCoord]

    # ===========================================================================

    def profile(self):
        colFreqData = []
        for i in range(len(self.aligned_sequences[0])):
            uniq_chars, counts = np.unique(
                self.aligned_sequences[:, i], return_counts=True
            )
            freqDict = dict(zip(uniq_chars, counts))
            colFreqData.append([freqDict, i])
        return colFreqData

    # ===========================================================================

    def normalizeSeqRefs(self):
        """
        In certain instances, additional metadata such as sequence start, end,
        orientation, and name may be parsed/stored in the 'name' metadata field.
        This function detects the most common structure of this convention
        (the name of which I am blanking) and removes this metadata from the name
        field and moves it to its proper storing location.

        TODO: there is typo in the original, fix if you like..
        Should we add more specific errors, or just leave it since this
        is not widely used?
        """
        # e.g.: chrUn_KK085329v1_11857_12355_R
        formatToLookFor = r"^(?P<name>\S+)_(?P<start>\d+)_(?P<end>\d+)(?P<orient>_R)?"
        for refMeta in self.seqMetaData:
            if re.match(formatToLookFor, refMeta["name"]):
                curStart = refMeta["src_start"]
                curEnd = refMeta["src_end"]
                if (
                    curStart == None
                    or curEnd == None
                    or refMeta["src_orientation"] == None
                ):
                    raise valueError(
                        "Error: normalizeSeqRefs() requires that all "
                        + "sequence data have start, end and orientation metaData "
                        + "supplied - some of these sequences do not!\n"
                    )
                nameData = re.match(formatToLookFor, refMeta["name"]).groupdict()
                if nameData["orient"] != None:
                    refMeta["src_end"] = int(nameData["end"]) - curStart + 1
                    refMeta["src_start"] = int(nameData["end"]) - curEnd + 1

                    if refMeta["src_orientation"] == "+":
                        refMeta["src_orientation"] = "-"
                    else:
                        refMeta["src_orientation"] = "+"
                else:
                    refMeta["src_start"] = int(nameData["start"]) + curStart + 1
                    refMeta["src_end"] = int(nameData["start"]) + curEnd - 1
                refMeta["name"] = nameData["name"]

    # ===========================================================================

    def printAlignments(
        self,
        blockSize=50,
        includeRef=False,
        showScore=False,
        showCons=False,
        origOrder=False,
    ):

        """
        Function that simply prints the alignments to the screen.

        Args:
            blocksize: accepts int that determines the number of bp to print before
                        printing a newline (default 50)

            includeRef: accepts boolean that determines whether or not the reference
                        should be printed (default false)

            showScore: accepts boolean that determines whether or not the score
                        will be printed (scores currently generated from
                        getLowScoringAlignmentColumns)

            showCons: accepts a boolean that determines whether or not the consensus
                        should be printed (default false)

            origOrder: accepts a boolean that determines whether or not the sequences
                        should be printed in order of 'src_start' (ascending)
        TOD0:
        The line marked 'BEASTLINE' is the love child of these two:
        https://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list
        https://www.geeksforgeeks.org/ways-sort-list-dictionaries-values-python-using-lambda-function
        I am open to other suggestions for it!
        """

        gotCoordinates = True
        gotRef = True
        gotRefCoordinates = True
        consensus = None

        if showCons:
            consensus = self.consensus()

        sortedIndexes = list(range(0, len(self.aligned_sequences)))
        if self.seqMetaData[0]["src_start"] is not None:
            if origOrder == False:
                # BEASTLINE here VVVV
                sortedIndexes = [
                    x
                    for _, x in sorted(
                        zip(self.seqMetaData, sortedIndexes),
                        key=lambda pair: pair[0]["src_start"],
                    )
                ]
        else:
            warnings.warn(
                "MultAlign Warning: printAlignments(): print function "
                + "was called without aligned sequence start and end "
                + "metadata: start and end bp coordinates will be replaced "
                + "with counts...\n"
            )
            gotCoordinates = False

        maxCoordLen = 1
        maxIDLen = 1
        if self.reference is not None:
            maxIDLen = len(self.refMetaData["name"])
            if self.refMetaData["src_end"] is not None:
                maxCoordLen = len(str(self.refMetaData["src_end"]))
            else:
                gotRefCoordinates = False

        else:
            warnings.warn(
                "MultAlign Warning: printAlignments: print function was "
                + "called without reference or ref meta data - ref will be "
                + "replaced with consensus...\n"
            )
            gotRef = False

        for seqData in self.seqMetaData:
            if len(seqData["name"]) > maxIDLen:
                maxIDLen = len(seqData["name"])
            if gotCoordinates:
                if len(str(seqData["src_start"])) > maxCoordLen:
                    maxCoordLen = len(str(seqData["src_start"]))
                if len(str(seqData["src_end"])) > maxCoordLen:
                    maxCoordLen = len(str(seqData["src_end"]))
                    print(len(str(seqData["src_end"])))
            else:
                maxCoordLen = len(self.aligned_sequences[0])

        # we have to get length of maximum score to adjust bp printing
        # (and the scores themselves...)
        scores = None
        maxScoreLen = 0
        if showScore == True:
            scores = self.getLowScoringAlignmentColumns(justScores=True)
            maxScoreLen = 1
            for i, score in enumerate(scores):
                score = str(score)
                scores[i] = score
                if len(score) > maxScoreLen:
                    maxScoreLen = len(str(score))

            for i, score in enumerate(scores):
                scores[i] = score.ljust(maxScoreLen)

        lineStart = 0
        lineEnd = blockSize - 1
        if gotRef and gotRefCoordinates:
            refBaseStartPos = self.refMetaData["src_start"]
        else:
            refBaseStartPos = 1
        consBaseStartPos = 1
        if gotRef:
            reference = self.reference_seq
        else:
            reference = self.consensus()
        counter = 0

        # TODO: this may be subject to optimazation, but it's the best I've got on the fly
        # without having aligned_start and end...
        endings = [None] * len(
            self.aligned_sequences
        )  # we will use this to stop printing sequences that have ended...
        while lineStart < len(reference):
            # Show the score if requested:
            if showScore == True:
                seq = " ".join(scores[lineStart : lineStart + blockSize])
                name = "Scores"[0:maxIDLen]
                end = lineStart + blockSize - 1

                outStr = (
                    name.ljust(maxIDLen)
                    + " ".rjust(maxCoordLen)
                    + "  "
                    + seq
                    + "".ljust(blockSize - len(seq))
                    + "    "
                    + "\n"
                )
                print(outStr)
                scoreBaseStartPos = end + 1

            # print consensus if requested (unless we have no ref, in which case
            # the consensus will be printed below in the ref's place...)
            if showCons == True and gotRef == True:
                seq = "".ljust(maxScoreLen, " ").join(
                    consensus[lineStart : lineStart + blockSize]
                )
                name = "Consensus"[0:maxIDLen]
                numLetters = len(seq) - seq.count(" ")
                end = consBaseStartPos + numLetters - 1

                outStr = (
                    name.ljust(maxIDLen)
                    + " "
                    + str(consBaseStartPos).rjust(maxCoordLen)
                    + " "
                    + seq
                    + "".ljust(blockSize - len(seq))
                    + "    "
                    + str(end)
                    + "\n"
                )
                print(outStr)
                consBaseStartPos = end + 1

            # print out the reference:
            seq = "".ljust(maxScoreLen, " ").join(
                reference[lineStart : lineStart + blockSize]
            )
            if gotRef:
                name = ("ref: " + self.refMetaData["name"])[0:maxIDLen]
            else:
                name = "Cons"
            numLetters = (
                len(seq)
                - seq.count(self.internalGapChar)
                - seq.count(self.externalGapChar)
            )
            end = refBaseStartPos + numLetters - 1

            outStr = (
                name.ljust(maxIDLen)
                + " "
                + str(refBaseStartPos).rjust(maxCoordLen)
                + " "
                + seq
                + "".ljust(blockSize - len(seq))
                + "    "
                + str(end)
                + "\n"
            )
            print(outStr)
            refBaseStartPos = end + 1

            # Let's print the seqs
            for index in sortedIndexes:
                if endings[index] == None:
                    if gotCoordinates:
                        seqBaseStartPos = self.seqMetaData[index]["src_start"]
                    else:
                        seqBaseStartPos = 1
                else:
                    seqBaseStartPos = endings[index]
                seq = "".ljust(maxScoreLen, " ").join(
                    self.aligned_sequences[index][lineStart : lineStart + blockSize]
                )
                name = self.seqMetaData[index]["name"]
                numLetters = (
                    len(seq)
                    - seq.count(self.internalGapChar)
                    - seq.count(self.externalGapChar)
                )
                end = seqBaseStartPos + numLetters - 1

                if end < seqBaseStartPos:
                    continue

                outStr = (
                    name.ljust(maxIDLen)
                    + " "
                    + str(seqBaseStartPos).rjust(maxCoordLen)
                    + " "
                    + seq
                    + "".ljust(blockSize - len(seq))
                    + "    "
                    + str(end)
                    + "\n"
                )
                print(outStr)
                counter += 1
                endings[index] = end + 1

            print("\n")
            lineStart = lineEnd + 1
            lineEnd = lineStart + blockSize - 1

    # ===========================================================================

    def getAlignmentBlock(self, start, end, rawSequences=False):
        """
        Function that returns a 'block' of the alignment in form of a 2d numpy array.

        Args:
            start: accepts int that determines the desired starting index of the
                    block to be returned

            end: accepts int that determines the desired end index of the block
                    to be returned

            rawSequences: TODO: what is this?

        #TODO: The original uses a for loop and grabs any sequence with index >= start
        #and <= the end, grabbing 1 more than a simple slice. Is this grandfathered in,
        #or should we change to just be standard slice like everything else - this
        #would also allow for negative numbers to be used (like -1 for last bp)
        #FIXME: This currently returns a numpy array of lists - because you can't have
        #2d numpy arrays with funky dimensions - only when rawSequences = True, but
        #returns a 2d numpy array otherwise - not good!
        #ALSO the original only checked if rawSequences was defined, I have changed it
        #to a boolean value - keep?
        """
        endSlice = end + 1
        results = self.aligned_sequences[start:endSlice]

        if rawSequences != False:
            seqLength = len(results[0])
            holder = []
            seq = []
            counter = 0
            for base in np.nditer(results):

                if base != " ":
                    seq.append(str(base))
                counter += 1
                if counter % seqLength == 0:
                    holder.append(seq)
                    seq = []
            results = np.array(holder)
        return results

    # ===========================================================================
    def consensus(
        self,
        parameter_set="Linup",
        cpg_adjustment=True,
        with_inserts=False,
        only_ref_sites=False,
    ):
        """
        Call the consensus of this multiple alignment. See ConsensusCaller.py for
        implementation details and the allowed values for paramter_set.

        Parameters:
            parameter_set    : Only "Linup" is supported at this time.

            cpg_adjustment   : Whether or not to apply the CpG adjustment algorithm.

            with_inserts     : When True, returns a consensus sequence padded with
                               gaps ('-') such that it lines up with the rest
                               of the alignment

            only_ref_sites   : Rarely used. When True, only attempt to call sites
                               that were already called as non-gaps in the reference.
                               In other words, do not try to call existing gaps even if
                               there is evidence for calling a base.

        Returns:
            The new consensus sequence, as a string.
        """

        if only_ref_sites:
            if self.reference is not None:
                # A reference sequence, in this context, is either a single sequence
                # from which all sequences have been aligned and padded to correspond
                # column-for-column with the MSA, or a special sequence which denotes
                # the "true" columns vs the insertion columns.  A good example of the
                # latter is the ref line in a HMMER stockholm file which uses "X"/"x" or
                # "-"/"." to characterize the match or insertion columns within the MSA.
                # The match states may also be denoted using the consensus base.
                only_positions = [
                    i for (i, n) in enumerate(self.reference[:]) if n not in "-. "
                ]
            else:
                raise Exception(
                    "Asked to call only_ref_sites, but no reference was provided"
                )
        else:
            only_positions = None

        consensus, cons_positions = ConsensusCaller._call_consensus_ndarray(
            self.aligned_sequences,
            parameter_set,
            cpg_adjustment=cpg_adjustment,
            only_positions=only_positions,
        )

        if with_inserts:
            new_consensus = ["-"] * self.aligned_sequences.shape[1]
            for (c, i) in zip(consensus, cons_positions):
                new_consensus[i] = c
        else:
            new_consensus = consensus

        return "".join(new_consensus)

    # ===========================================================================

    def getLowScoringAlignmentColumns(
        self,
        matrix="default",
        gapInitiationPenalty=-40,
        gapExtensionPenalty=-15,
        threshold=1,
        justScores=False,
    ):
        """
        Generate a list of low blocks within the multiple alignment.

        matrix               : A reference to a SequenceSimilarityMatrix
                                    object.
          gapInitiationPenalty : The penalty to initiate a gap.
          gapExtensionPenalty  : The penalty to extend a gap.
          threshold            : The maximum score for which to report
                                  low scoring blocks (DEFAULT: 1).
          columns              : A collection of low scoring columns
                                  with start/end position:
                                  startPosition =
                                    $columns->[ 0-numColumns ]->[ 0 ]
                                  endPosition =
                                    $columns->[ 0-numColumns ]->[ 1 ]
          valArray             : A reference to an array of block scores
                                  across the reference sequence.

        The low scoring subsequences are found by generating the average
        substitution score per column using a given scoring system, inverting
        the scores and using the Ruzzo & Tompa algorithm: "A Linear Time
        Algorithm For Finding All Maximal Scoring Subsequences" to identify
        the blocks.

        For instance, given the multiple alignment:

              AACT-TTGACC---CCacta
            GAAAGT-TTCTCCAGTCCacta
            G-AACTATTC-CCA--CCacta
            GAAACT-TTG-CC-----acta

        And the reference:

            GAAACT-TTN-CCA--CCACTA

        The first column has three aligned sequences.  If G to G is scored
        as +10 the average score for column one is 30/3 = +10.  The seventh
        column however has four sequences aligning (with deletions) and
        therefore would be scored as 0 + 0 + -40 + 0 = -40 assuming the gap
        open penalty is -40.  The average for the columne is -40/4 = -10.

        For this alignment the column score array would be:

            10, -7.33, 9, 9, 3.75, 9, -10, 9, 9 -1, -20, 10, 10, -15.5, -10,
            -3.75, 3.75, 3.75, 9, 10, 9, 9

        Inverting the scores and and applying the Ruzzo Tompa algorithm will
        give:

            [[ 7.33 ], [ 10 ], [ 1, 20, -10, -10, 15.5, 10, 3.75 ]]

        With the following block ranges:

            [[ 1, 1 ], [ 6, 6 ], [ 9, 15 ]]

        This routine will report any block above the supplied threshold.  Using
        the default threshold of 1.0 all the above blocks will be returned in effect
        flagging the following columns of the alignment as low scoring:

        low: *    *  *******
              AACT-TTGACC---CCacta
            GAAAGT-TTCTCCAGTCCacta
            G-AACTATTC-CCA--CCacta
            GAAACT-TTG-CC-----acta

        """
        gapChars = ("-", ".", " ")
        profile = [0] * (len(self.reference))
        posCounts = [0] * (len(self.reference))
        matrix_r = {}

        if isinstance(matrix, SubstitutionMatrix):
            alphaArray = matrix.matrix
            alphabet_r = matrix.alphabet_r
            # TODO: Test that this works with matrix, Substitution matrix args
        else:
            # A    R    G   C    Y    T    K     M   S     W    N   V   H   D   B
            alphaArray = [
                [
                    9,
                    1,
                    -6,
                    -15,
                    -16,
                    -17,
                    -12,
                    -2,
                    -10,
                    -4,
                    -1,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    1,
                    1,
                    1,
                    -15,
                    -15,
                    -16,
                    -6,
                    -6,
                    -6,
                    -7,
                    -1,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -6,
                    1,
                    10,
                    -15,
                    -15,
                    -15,
                    -2,
                    -10,
                    -2,
                    -10,
                    -1,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -15,
                    -15,
                    -15,
                    10,
                    2,
                    -6,
                    -9,
                    -2,
                    -2,
                    -9,
                    -1,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -16,
                    -15,
                    -15,
                    1,
                    1,
                    1,
                    -6,
                    -7,
                    -7,
                    -7,
                    -1,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -17,
                    -16,
                    -15,
                    -6,
                    1,
                    9,
                    -2,
                    -12,
                    -11,
                    -4,
                    -1,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -12,
                    -6,
                    -2,
                    -11,
                    -6,
                    -2,
                    -2,
                    -11,
                    -7,
                    -7,
                    -1,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -2,
                    -6,
                    -10,
                    -2,
                    -7,
                    -12,
                    -11,
                    -2,
                    -7,
                    -7,
                    -1,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -10,
                    -6,
                    -2,
                    -2,
                    -7,
                    -11,
                    -7,
                    -7,
                    -2,
                    -10,
                    -1,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -4,
                    -7,
                    -10,
                    -11,
                    -7,
                    -4,
                    -7,
                    -7,
                    -10,
                    -4,
                    -1,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
                [
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                    -2,
                ],
            ]

            alphabet_r = [
                "A",
                "R",
                "G",
                "C",
                "Y",
                "T",
                "K",
                "M",
                "S",
                "W",
                "N",
                "V",
                "H",
                "D",
                "B",
            ]

        for i in range(len(alphabet_r)):
            for j in range(len(alphabet_r)):
                matrix_r[(alphabet_r[i], alphabet_r[j])] = alphaArray[i][j]

        for seqNum, querySeq in enumerate(self.aligned_sequences):
            inGap = False
            for baseNum, queryBase in enumerate(querySeq):
                if baseNum >= self.getAlignedStart(
                    seqNum
                ) and baseNum <= self.getAlignedEnd(seqNum):
                    refBase = self.reference[baseNum]
                    if (refBase in gapChars and queryBase not in gapChars) or (
                        queryBase in gapChars and refBase not in gapChars
                    ):
                        if inGap:
                            profile[baseNum] += gapExtensionPenalty
                            posCounts[baseNum] += 1
                        else:
                            profile[baseNum] += gapInitiationPenalty
                            posCounts[baseNum] += 1
                            inGap = True

                    elif refBase in gapChars and queryBase in gapChars:
                        posCounts[baseNum] += 1

                    # don't execute if we have not reached the aligned start or have passed the aligned end
                    # This statement was not needed in the previous version due to the
                    # OG data structure not storing leading / trailing sequence gaps
                    else:
                        # print("Profile / pos counts index: " + str(baseNum ))
                        profile[baseNum] += matrix_r[
                            (refBase.upper(), queryBase.upper())
                        ]
                        posCounts[baseNum] += 1
                        inGap = False

                # Now exiting big double for loop

        for j in range(len(posCounts)):
            if posCounts[j] > 0:
                profile[j] = profile[j] / posCounts[j]

        # I have added this so we
        # 1. Don't have to add seperate function with repeated code if said funciton
        #   won't be used
        # 2. Can get scores to test the 'showScores' portion of the printAlignment
        #   function (until we have a way of getting other scores)
        # Fixme: keep? (if you delete, fix printAlignments() too)
        if justScores:
            return profile

        # to find worst scoring regions, must calculate inverse of profile signs
        for i in range(len(profile)):
            profile[i] = profile[i] * -1

        # TODO: For debugging -- remove when done
        # print("Profile: " + str(profile))
        # expectedProfile = [6, 23.33, 6, 17, 12.75, 6, 30, 15, -9, 10.5, 20, -10, 15, 20, 17.5, 11.25, 8.25, 15, 0, 0, 0, 0]
        # for i, val in enumerate(expectedProfile):
        #    dif = val - profile[i]
        #    if dif > 0.1:
        #        print("Index of failure: " + str(i) + " Expected: " + str(val) + "  \tGot: " + str(profile[i]))
        """
        bases = "ACTG"
        for char in bases:
            for char2 in bases:
                 print(char + char2 + " = " + str(matrix_r[(char, char2)]))
        """
        # Deprecated function replaced with a validated one
        # myValArray = self._ruzzoTompaBasedScoreSequences(profile)
        RTA, IA, myValArray = MultAlign._ruzzoTompaFindAllMaximalScoringSubsequences(
            profile
        )

        seqPosAvg = myValArray.copy()
        columns = []

        # find low quality columns:
        inCol = 0
        colStart = -1
        for i, val in enumerate(
            seqPosAvg
        ):  # Find start and end positions where score is greater than threshold
            if val >= threshold:
                if colStart == -1:
                    inCol = 1
                    colStart = i

            elif inCol == 1:
                inCol = 0
                columns.append(
                    [colStart, i - 1]
                )  # FIXD: what is even going on in OG code here?
                colStart = -1  # Fixd: append, don't extend

            # end for loop
        if inCol == 1:
            columns.append([colStart, seqPosAvg])

        return (columns, myValArray)

    # ===========================================================================

    def _ruzzoTompaFindAllMaximalScoringSubsequences(scores):
        """
        An implementation of an algorithm by Walter Ruzzo and Martin
        Tompa ("A Linear Time Algorithm for Finding All Maximal Scoring
        Subsequences" - 7th Intl Conf. Intelligent Systems for Mol
        Biology). This was based on the python implementation found
        on the wikipedia page for this algorithm.

        Given the score array example: 4,-5,3,-3,1,2,-2,2,-2,1,5

        Identify the maximal scoring subsequences.  This is traditionally
        returned as a list of lists. Each list contains the position ordered
        scores for a maximal subsequence.  For the example above it would
        be:
               [[4], [3], [1, 2, -2, 2, -2, 1, 5]]

        Additionally this implementation provides the interval as a list
        of [start, end] lists.  For the example above:

               [[0, 1], [2, 3], [4, 11]]

        Finally this also provides a score mask.  This is list the size
        of the input score array where each position is either 0 or it
        contains the score of the maximal scoring range it is a member of.
        For instance in the above example:

               [4, 0, 3, 0, 7, 7, 7, 7, 7, 7, 7]

        Arg: Scores: accepts a single 1d list of scores to be processed
        """
        k = 0
        total = 0
        # Allocating arrays of size n
        I, L, R, Lidx = [[0] * len(scores) for _ in range(4)]
        for i, s in enumerate(scores):
            total += s
            if s > 0:
                # store I[k] by (start,end) indices of scores
                I[k] = (i, i + 1)
                Lidx[k] = i
                L[k] = total - s
                R[k] = total
                while True:
                    maxj = None
                    for j in range(k - 1, -1, -1):
                        if L[j] < L[k]:
                            maxj = j
                            break
                    if maxj is not None and R[maxj] < R[k]:
                        I[maxj] = (Lidx[maxj], i + 1)
                        R[maxj] = total
                        k = maxj
                    else:
                        k += 1
                        break
        # There are now k-1 valid ranges in I
        score_mask_array = [0] * len(scores)
        interval_array = []
        for i in range(k):
            interval_array.append([I[i][0], I[i][1]])
            for l in range(I[i][0], I[i][1]):
                score_mask_array[l] = R[i] - L[i]
        return (
            [scores[I[l][0] : I[l][1]] for l in range(k)],
            interval_array,
            score_mask_array,
        )

    ##
    ## DEPRECATED see above
    ##
    def _ruzzoTompaBasedScoreSequences(self, scores):
        """
        This algorithm is the same as the
        '_ruzzoTompaFindAllMaximalScoringSubsequences' from the previous implementation,
        which uses an algorithm from a similiarly named paper.

        The name of the function has been changed to reflect the fact that neither
        implementation is actually a direct implementation of Ruzzo and Tompa's
        algorithm, but rather an algorithm based on said algorithm.

        The orginal algorithm published by R&T would return a list of the maximal
        scoring subsequences of the input scores, so input scores of :
        [58, 85, 24, 68, 51, 24, 40, 60, -36, 42, 80, -40, 60, 80, 15, 15, 58, 60, -36, -40, -36, -36, 0, 0]
        would have an output of:
        [[58, 85, 24, 68, 51, 24, 40, 60, -36, 42, 80, -40, 60, 80, 15, 15, 58, 60]]

        On the other hand, this algorithm (with same input) would output:
        [744, 744, 744, 744, 744, 744, 744, 744, 744, 744, 744, 744, 744, 744, 744, 744, 744, 60, 0, 0, 0, 0, 0, 0]
        which is a list the length of the scores list, where each index holds the
        score of the maximally scoring subsequence to which it belongs

        The computation of this 'ScoreSequence' cannot be removed from this function
        because it requires intermediary data from R&T's original algorithm

        Arg: Scores: accepts a single 1d list of scores to be processed

        #TODO: Validate output (especially after R&T is done)!!!!!!!!
        #FIXME: ^^^^^^^^^^^^^^^^
        """
        k = 0
        total = 0
        # Allocating arrays of size n
        I, L, R, Lidx = [[0] * len(scores) for _ in range(4)]
        for i, s in enumerate(scores):
            total += s
            if s > 0:
                # store I[k] by (start,end) indices of scores
                I[k] = (i, i + 1)
                Lidx[k] = i
                L[k] = total - s
                R[k] = total
                while True:
                    maxj = None
                    for j in range(k - 1, -1, -1):
                        if L[j] < L[k]:
                            maxj = j
                            break
                    if maxj is not None and R[maxj] < R[k]:
                        I[maxj] = (Lidx[maxj], i + 1)
                        R[maxj] = total
                        k = maxj
                    else:
                        k += 1
                        break
        # Getting maximal subsequences using stored indices
        ruzzoTompaResults = [scores[I[l][0] : I[l][1]] for l in range(k)]
        # This concludes RuzzoTompa algorithm

        output = [0] * len(scores)
        for i in range(len(I)):
            if I[i] != 0:
                for j in range(I[i][0], I[i][1]):
                    output[j] = R[i] - L[i]

        """
        print("\n\nHere is Ruzzo Tompa Results")
        print(ruzzoTompaResults)
        print("\nHere is Scores")
        print(scores)
        print("\nHere is output")
        print(output)
        """

        return output

    # ===========================================================================

    def instRefStart(self, seqNum):
        """
        Returns the reference base pair coordinate (no gaps) where the specified aligned
        Sequence begins. E.g.:
        Ref     ACGT---ACGTACGT
        Seq1           ACGTACGT

        The alignedStart of seq1 is 7 (Its 1st bp is at the gapped reference's
        7th position) while the instRefStart is 5 (the 7th index of the ref and
        seq1 is the 5th bp of Ref)

        Arg: seqNum: accepts the index of the seq in aligned_sequences that the
                        instRefStart will be computed for.
        """
        gappedStart = self.getAlignedStart(seqNum)
        relevantRefPortion = (
            self.reference[0 : (gappedStart + 1)]
            .strip(internalGapChar)
            .strip(externalGapChar)
        )
        ungappedStart = len(relevantRefPortion)

    # ===========================================================================

    def _pickHistogramPeaks(
        self,
        windowSize=11,
        threshold=None,
        perPosCount=None,
        useHighestInWindow=False,
        histogram=None,
    ):
        """
        This function is untested and unreviewed, and is based on a (potentially
        bugged) implementation of the same name in the original that was also
        experimental, and was only called by a function that was also experimental
        and not used (which is not yet implemented here at this time) - so its
        testing and implementation of debug options found in the original has been
        halted...
        """

        score = 0
        peakHighScore = 0
        peakPos = -1
        inPeak = 0
        peakWindowStart = -1
        peakList = []
        windowFlankingSize = (windowSize - 1) / 2

        windowNameHash = {}
        for i in range(windowFlankingSize):
            score += histogram[i]
            for j in range(len(self.aligned_seuqences)):
                if self.instRefStart(j) == i:
                    if self.seqMetaData[i]["name"] in windowNameHash:
                        windowNameHash[self.seqMetaData[i]["name"]] += 1
                    else:
                        windowNameHash[self.seqMetaData[i]["name"]] = 0

        i = -1
        numUniqSequences = -1
        sig = -1
        distSinceLastCall = windowFlankingSize

        for i in range(len(histogram)):
            windowStart = i - windowFlankingSize
            windowEnd = i + windowFlankingSize
            distSinceLastCall += 1

            if windowEnd < len(histogram):
                score += histogram[windowEnd]

            if (windowStart - 1) >= 0:
                score -= histogram[windowStart - 1]

            for j in range(len(self.aligned_sequences)):
                if self.instRefStart(j) == windowEnd:
                    windowNameHash[self.seqMetaData[j]["name"]] += 1

                elif self.instRefEnd(j) == (windowStart - 1):
                    windowNameHash[self.seqMetaData[j]["name"]] -= 1
                    if windowNameHash[self.seqMetaData[j]["name"]] == 0:
                        windowNameHash.pop(self.seqMetaData[j]["name"], None)

            numUniqSequences = len(windowNameHash)
            sig = 0
            if numUniqSequences >= 3:
                sig = round((score / numUniqSequences), 2)

            # TODO: add that debug option here (and everywhere below that says 'add debug')

            if inPeak == 0 and sig > 0.5 and distSinceLastCall > windowFlankingSize:
                # add debug
                peakWindowStart = 0
                if windowStart >= 0:
                    peakWindowStart = windowStart
                inPeak = 1

            if inPeak == 1 and sig <= 0.5:
                # add debug
                max = 0
                maxPos = -1

                for j in range(peakWindowStart, i):
                    if max < histogram[j]:
                        max = histogram[j]
                        maxPos = j

                peakDict = {}
                peakDict["pos"] = maxPos
                peakDict["score"] = peakHighScore
                peakList.append(peakDict)

                inPeak = 0
                distSinceLastCall = 0
                peakHighScore = -1

            if inPeak == 1 and peakHighScore < score:
                peakHighScore = score

        # trailing case
        if inPeak == 1:
            # add debug
            max = 0
            maxPos = -1
            for j in range(peakWindowStart, i):
                if max < histogram[j]:
                    max = histogram[j]
                    maxPos = j

            peakDict = {}
            peakDict["pos"] = maxPos
            peakDict["score"] = peakHighScore
            peakList.append(peakDict)

        return peakList

    # ===========================================================================

    def filter(
        self,
        minGC=None,
        maxGC=None,
        minDiv=None,
        maxDiv=None,
        minSrcDiv=None,
        maxSrcDiv=None,
        removeSequencesLackingData=False,
    ):
        """
        All previous filter functions rolled into 1 function. Each argument implies
        the maximum or minimum allowable value for which it's named, e.g.:
        minGC = 2 will filter any sequences with a GC value lower than 2 (those
        equal to 2 will be spared)
        Leaving a value as None will not filter based on that variable, so maxCG
        of None will result in no upper bound for the cg value

        This function alters the data structure, removing all sequences (with the
        exception of the reference) that do not meet the specified filter ranges.
        This may alter the lengths of all sequences (ref included) if a sequence
        is removed that had the only non-gap character for that column of the
        alignmnent.

        By default, sequences that do not have metadata that is being filtered on
        will NOT be removed, though they can be by setting the arg
        'removeSequencesLackingData' to True. Sequences that lack metadata that is
        NOT being filtered on (the min and max remain None) will not be affected
        by this.

        TODO:
        --  FIX metadata
        """

        indicesToRemove = []
        for i, metaData in enumerate(self.seqMetaData):

            # Starting with CG
            if minGC != None:
                if metaData["CG"] != None:
                    if metaData["CG"] < minGC:
                        indicesToRemove.append(i)
                        continue  # don't bother checking others if we're cutting
                        # it anyway...
                elif removeSequencesLackingData:
                    indicesToRemove.append(i)
                    continue

            if maxGC != None:
                if metaData["GC"] != None:
                    if metaData["GC"] > maxGC:
                        indicesToRemove.append(i)
                        continue
                elif removeSequencesLackingData:
                    indicesToRemove.append(i)
                    continue

            # Now src_div
            if minSrcDiv != None:
                if metaData["srcDiv"] != None:
                    if metaData["srcDiv"] < minSrcDiv:
                        indicesToRemove.append(i)
                        continue
                elif removeSequencesLackingData:
                    indicesToRemove.append(i)
                    continue

            if maxSrcDiv != None:
                if metaData["srcDiv"] != None:
                    if metaData["srcDiv"] > maxSrcDiv:
                        indicesToRemove.append(i)
                        continue
                elif removeSequencesLackingData:
                    indicesToRemove.append(i)
                    continue

            # Now div
            if minDiv != None:
                if metaData["Div"] != None:
                    if metaData["Div"] < minDiv:
                        indicesToRemove.append(i)
                        continue
                elif removeSequencesLackingData:
                    indicesToRemove.append(i)
                    continue

            if maxDiv != None:
                if metaData["Div"] != None:
                    if metaData["Div"] > maxDiv:
                        indicesToRemove.append(i)
                elif removeSequencesLackingData:
                    indicesToRemove.append(i)

        # This removes duplicates so we don't screw up the datastructe by repeatedly
        # removing the nth index if it failed multiple tests or lacked multiple
        # types of metadata
        indicesToRemove = list(set(indicesToRemove))

        # This removes the undesirable alignments
        self.aligned_sequences = np.delete(self.aligned_sequences, indicesToRemove, 0)

        # Here we remove the metadata at each index (ordered largest to smallest to
        # prevent index mixups)
        for index in sorted(indicesToRemove, reverse=True):
            del self.seqMetaData[index]

        # Finally, we must check if any columns in the 2d array are now just gaps,
        # and remove them
        colsToRemove = []
        for i in range(len(self.aligned_sequences[0])):
            unique_chars = np.unique(self.aligned_sequences[:, i])
            # print("\n")
            # print(uniq_chars)
            # print(counts)
            if len(unique_chars) == 1:
                if unique_chars[0] in (self.internalGapChar, self.externalGapChar):
                    colsToRemove.append(i)

        self.aligned_sequences = np.delete(self.aligned_sequences, colsToRemove, 1)
        self.reference = np.delete(self.aligned_sequences, colsToRemove, 1)

        return (indicesToRemove, colsToRemove)

    # ===========================================================================

    # TODO: Make a CpG 1,1/10 transition version...weighted average?
    def kimura_divergence(self, only_ref_sites=False):
        """ """
        # TODO: Move this up
        # Mutation Types:
        #      Purines
        #      A--i--G
        #      | \ / |
        #      v  v  v
        #      | / \ |
        #      C--i--T
        #    Pyrimidines
        #  i = Transitions ( more frequent )
        #  v = Transversions ( rarer )
        #
        #  This lookup structure encodes
        #  transitions as "1" and transversions
        #  as "2".
        mut_types = {
            "CT": 1,
            "TC": 1,
            "AG": 1,
            "GA": 1,
            "GT": 2,
            "TG": 2,
            "GC": 2,
            "CG": 2,
            "CA": 2,
            "AC": 2,
            "AT": 2,
            "TA": 2,
        }
        unambiguous_bases = {"A": 1, "C": 1, "G": 1, "T": 1}

        if only_ref_sites:
            if self.reference is not None:
                # A reference sequence, in this context, is either a single sequence
                # from which all sequences have been aligned and padded to correspond
                # column-for-column with the MSA, or a special sequence which denotes
                # the "true" columns vs the insertion columns.  A good example of the
                # latter is the ref line in a HMMER stockholm file which uses "X"/"x" or
                # "-"/"." to characterize the match or insertion columns within the MSA.
                # The match states may also be denoted using the consensus base.
                only_positions = [
                    i for (i, n) in enumerate(self.reference[:]) if n not in "-. "
                ]
            else:
                raise Exception(
                    "Asked to calculate divergence with only_ref_sites, but no reference was provided"
                )
        else:
            only_positions = None

        # TODO: take parameters for the consensus calling step
        parameter_set = "Linup"
        consensus, cons_positions = ConsensusCaller._call_consensus_ndarray(
            self.aligned_sequences,
            parameter_set,
            cpg_adjustment=True,
            only_positions=only_positions,
        )

        total_divergence = 0
        for seq in self.aligned_sequences:
            total_pairs = 0
            transitions = 0
            transversions = 0

            for (cons_base, pos) in zip(consensus, cons_positions):
                seq_base = seq[pos]
                if seq_base == " ":
                    continue

                if cons_base in unambiguous_bases and seq_base in unambiguous_bases:
                    total_pairs += 1

                    if mut_types.get(cons_base + seq_base) == 1:
                        transitions += 1
                    elif mut_types.get(cons_base + seq_base) == 2:
                        transversions += 1

            kimura_div = 1.0
            if total_pairs >= 1:
                p = transitions / total_pairs
                q = transversions / total_pairs
                # Problem case 1: when (1-2q)<0, its square root is no longer in the reals.
                if q <= 0.5:
                    log_operand = (1 - (2 * p) - q) * (1 - (2 * q)) ** 0.5
                    # Problem case 2: when p and q are large enough, the log_operand can end up negative
                    if log_operand > 0:
                        kimura_div = -0.5 * math.log(log_operand)

            total_divergence += kimura_div

        return total_divergence / len(self.aligned_sequences)

    # ===========================================================================
    def __repr__(self):
        """
        __repr__() - Generic representation of an object in JSON format
        """
        return json.dumps(self.__dict__, indent=4)

    # ===========================================================================
    """
    MultAlign object is subscriptalbe - using multAlign[i] will return a string
    version of the ith aligned sequence, while multAlign['ref'] (or 'reference')
    will return a string version of the reference
    """

    def __getitem__(self, item):
        if item == "ref" or item == "reference":
            return "".join(self.reference)
        elif type(item) == int:
            return "".join(self.aligned_sequences[item])
        else:
            raise ValueError(
                "\nMultAlign Subscript Error: MultALign object is "
                + "subscriptable, but subscript must either be 'reference'"
                + ", 'ref', or index of desired aligned sequence"
            )

    # ===========================================================================

    # TODO: I wonder if 'length' shouldn't be the number of aligned bp, and
    # something else return the number of sequences? Just a thought
    # TODO: perhaps some of these could be normal properties - tbd
    @property
    def length(self):
        return len(self.aligned_sequences)

    @property
    def colLength(self):
        return len(self.aligned_sequences[0])

    @property
    def sequence(self):
        return self.aligned_sequences

    @property
    def reference_name(self):
        return self.refMetaData["name"]

    @property
    def reference_seq(self):
        return str("".join(self.reference))

    @property
    def refStart(self):
        return self.refMetaData["startPos"]

    @property
    def refEnd(self):
        return self.refMetaData["endPos"]

    # ===========================================================================

    def serializeOuptut(self, newFileName, relFolderPath=None, newFolderName=None):
        """
        This function supports several options for its output file.
        By default, the file is created in the working directory.
        It can also be stored:

        --in a new folder in the working direcotry (just give new folder name),
        --directory of choice (just specify a relFolderPath)
        --in a new directory within a specified directory (give folder name and folder path)


        TODO: Sean, this function could be further optimized - if Robert wants to keep
        all this functionality, fix it up nicely - Sean
        """
        sequenceHolder = self.aligned_sequences
        refHolder = self.reference
        self.aligned_sequences = self.aligned_sequences.tolist()
        for i in range(len(self.aligned_sequences)):
            self.aligned_sequences[i] = "".join(self.aligned_sequences[i])
        self.reference = "".join(self.reference.tolist())
        if newFolderName == None and relFolderPath == None:  # store file in working dir
            newOutputFile = open(newFileName, "x")
            newOutputFile.write(self.__repr__())
            newOutputFile.close()
        elif (
            newFolderName != None and relFolderPath == None
        ):  # store file in new dir created in working dir
            currentDirectory = os.getcwd()
            path = os.path.join(currentDirectory, newFolderName)
            try:
                os.mkdir(path)
            except OSError as error:
                print(error)
            newFile = os.path.join(path, newFileName)
            newOutputFile = open(newFile, "x")
            newOutputFile.write(self.__repr__())
            newOutputFile.close()
        elif (
            newFolderName == None and relFolderPath != None
        ):  # store file in existing folder at path
            print(os.getcwd())
            currentDirectory = os.getcwd()
            totalPath = os.path.join(currentDirectory, relFolderPath)
            newFile = os.path.join(totalPath, newFileName)
            newOutputFile = open(newFile, "x")
            newOutputFile.write(self.__repr__())
            newOutputFile.close()
        elif newFolderName != None and relFolderPath != None:
            currentDirectory = os.getcwd()
            pathToNewDirectory = os.path.join(currentDirectory, relFolderPath)
            newDirectory = os.path.join(pathToNewDirectory, newFolderName)
            try:
                os.mkdir(newDirectory)
            except OSError as error:
                print(error)
            newFile = os.path.join(newDirectory, newFileName)
            newOutputFile = open(newFile, "x")
            newOutputFile.write(self.__repr__())
            newOutputFile.close()
        self.aligned_sequences = sequenceHolder
        self.reference = refHolder

    # ===========================================================================

    # Todo: add option to not check for illegal chars
    def serializeInput(self, fileToReadFrom, deleteFile=False):
        objectFile = open(fileToReadFrom, "r")
        data = json.load(objectFile)
        objectFile.close()
        allowed_keys = set([])
        self.__dict__.update((k, None) for k in allowed_keys)
        self.__dict__.update((k, v) for k, v in data.items() if k in allowed_keys)

        if "aligned_sequences" in data:
            # Sequences should already be padded to be the same length
            seqs = data.get("aligned_sequences", None)

            # It's faster to build out a python 2d array and then convert it to a
            # ndarray than to call ndarray.append()
            aligned_data = []
            for seq in seqs:
                # print("Adding sequence with length: " + str(len(seq)))
                if isinstance(seq, Sequence):
                    aligned_data.append(
                        [ch.upper() for ch in seq.sequence_as_unicode_array]
                    )
                elif isinstance(seq, str):
                    aligned_data.append(list(seq.upper()))
            self.aligned_sequences = np.asarray(aligned_data)

            # Toss out our garbage
            aligned_data = None
            seqs = None
            # print(self.aligned_sequences)
            # print("dtype is: ",self.aligned_sequences.dtype)
            # np.set_printoptions(threshold=np.nan)
        if "reference" in data:
            self.reference = np.array(list(data.get("reference", None).upper()))
        if deleteFile == True:
            os.remove(fileToReadFrom)

        if "seqMetaData" in data:
            self.seqMetaData = data.get("seqMetaData", None)

        if "refMetaData" in data:
            self.refMetaData = data.get("refMetaData", None)

    # ===========================================================================
