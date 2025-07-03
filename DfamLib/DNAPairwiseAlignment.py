# -*- coding: utf-8 -*-
"""
    DNAPairwiseAlignment : A class representing a generic DNA pairwise
                           sequence alignment.

SEE ALSO: Based on previous RepeatMasker::SearchResult.pm

AUTHOR(S):
    Sean Rice <sean.rice@isbscience.org>
    Robert Hubley <rhubley@isbscience.org>

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
import math
import re
import json
import logging

LOGGER = logging.getLogger(__name__)


class DNAPairwiseAlignment:
    """
    DNA pairwise sequence alignment base class.

    A DNA to DNA sequence alignment may involve the product of a query sequence
    ( consensus sequence, profile Hidden Markov query etc ) vs. a target sequence or
    it may simply be result of searching two regions of a genome against each other.  In
    either case it is often useful to designate one of the sequences as the "query/model"
    ( what is being searched for ) and the other as the "subject/target" ( what is being
    searched in ).  In this class the terms "query" and "target" are used for
    to make a distinction when necessary.

    This class is designed to hold alignments from a wide range of DNA/DNA sequence aligners
    and therefore has some special conventions for the designation of the sequences, a fixed
    ordering of sequence ranges/orientations, and a wide range of optional metadata
    properties.

    Conventions:

      * Sequence statistics that are directional ( raw scores from asymmetric
        scoring matrices, complexity_adjusted score etc ) are assumed to be
        generated query-to-target unless specified otherwise using the
        'reference' property.

      * The orientation of the "query" sequence is required to be on the forward strand.  The
        "target" sequence strand is either forward or reverse and is determined by the
        orientation property.

      * The sequence alphabet is restricted to a-z, A-Z and "-".  The case of
        aligned sequences is preserved.

      * The sequence coordinates are 1-based, fully closed coordinates.


    Examples:

        An alignment in it's simplest form may simply be a set of ranges and
        a score:

        score       : 200
        query_id    : 'MIR'
        query_start : 1
        query_end   : 42
        orientation : '+'
        target_id   : 'chr2'
        target_start: 32832
        target_end  : 32881

        In this case only the coordinates are validated to be > 0 and no
        further checking is done.  In addition in this example the reference
        sequence ( the sequence that is scored against ) is assumed to be the
        query sequence.

        An alignment of a consensus sequence query vs a single genomic locus
        might simply contain the following details:

          Score = 150
          MIR#SINE/MIR         102 ACTTAACCTCTCTGTGCCTCAGTTTCCTCATCTGTAAAAT 141
                                        --- i              ii i   v     i
          chr1               26356 ACTTA---TTTCTGTGCCTCAGTTCTCCCATATGTAAGAT 26392

        It is up to the DNAPairwiseAlignment creator to decide which sequence
        in the above alignment is more usefully defined as the "query" and which one
        is the "target", while aligners typically have a preferred order, there
        are times when it is more useful to use the opposite ordering when
        building this object.  Here is how this might be defined in a
        DNAPairwiseAlignment object:

        score             : 150
        query_id          : 'MIR#SINE/MIR'
        query_start       : 102
        query_end         : 141
        aligned_query_seq : "ACTTAACCTCTCTGTGCCTCAGTTTCCTCATCTGTAAAAT"
        target_id         : 'chr1'
        target_start      : 26356
        target_end        : 26392
        sub_pct           : 0.0
        query_gap_pct     : 7.5
        target_gap_pct    : 16.22
        aligned_target_seq: "ACTTA---TTTCTGTGCCTCAGTTCTCCCATATGTAAGAT"

        Three additional (optional) statistics fields may be set for this
        alignment:

            sub_pct  -  The percentage of substitutions relative to the
                        aligned (non-gap) bp in the alignment.
            query_gap_pct  - The percentage of gap characters in the query
                             relative to the aligned query sequence length.
            target_gap_pct - The percentage of gap characters in the target
                             relative to the aligned target sequence length.

    Attributes:
        align_id           : string  - May be used for a globally unique
                                       identifier or to link alignment
                                       fragments together in a collection.
                                       [optional]
        score              : int     - A default integer score field for
                                       an alignment.  This is primarly used
                                       by the crossmatch_encode method.
                                       [required]
        e_value            : float   - E-value for alignment if provided
                                       by alignment method. [optional]
        p_value            : float   - P-value for alignment if provided
                                       by alignment method. [optional]
        bias               : float
        sub_pct            : float
        query_gap_pct      : float
        target_gap_pct     : float
        kimura_div         : float
        query_id           : string
        query_start        : int
        query_end          : int
        query_len          : int
        query_name         : string
        env_start          : int
        env_end            : int
        aligned_query_seq  : string
        orientation        : "+" or "-"
        target_id          : string
        target_start       : int
        target_end         : int
        target_len         : int
        aligned_target_seq : string
        matrix_name        : string
        pp_string          : string
        rf_string          : string

    Derived Attributes:
        sub_pct = number_of_substititons / aligned_query_length
        query_gap_pct = number_of_gap_characters_in_query / aligned_query_length
        target_gap_pct = number_of_gap_characters_in_target / aligned_target_length

    """

    def __init__(self, *args, **kwargs):
        allowed_keys = set(
            [
                "align_id",
                "score",
                "bit_score",
                "e_value",
                "p_value",
                "bias",
                "perc_sub",
                "perc_ins",
                "perc_del",
                "kimura_div",
                "query_id",
                "query_start",
                "query_end",
                "query_len",
                "query_name",
                "env_start",
                "env_end",
                "aligned_query_seq",
                "orientation",
                "target_id",
                "target_start",
                "target_end",
                "target_len",
                "aligned_target_seq",
                "matrix_name",
                "pp_string",
                "rf_string",
                "reference",
                "query_name",
                "sub_pct",
                "query_gap_pct",
                "target_gap_pct",
                "target_name",
            ]
        )
        self.minimumArgumentEnforcer(kwargs)
        if "check_types" not in kwargs or kwargs["check_types"] == True:
            self.keyArgumentTypeEnforcer(kwargs)
        self.sequenceFilter(
            kwargs.get("aligned_query_seq"), kwargs.get("aligned_target_seq")
        )
        self.__dict__.update((k, None) for k in allowed_keys)
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)
        self.subInsDelEvaluator()
        if self.reference == None:
            self.reference = "query"

    # Method used to calculate percentages of insertions and substitutions when
    # no values are given
    def subInsDelEvaluator(self):
        if self.perc_sub == None and self.sub_pct != None:
            self.perc_sub = self.sub_pct
        if self.reference == None or self.reference == "query":
            if self.perc_del == None and self.target_gap_pct != None:
                self.perc_del = self.target_gap_pct
            if self.perc_ins == None and self.query_gap_pct != None:
                self.perc_ins = self.query_gap_pct
        elif self.reference == "target":
            if self.perc_ins == None and self.target_gap_pct != None:
                self.perc_ins = self.target_gap_pct
            if self.perc_del == None and self.query_gap_pct != None:
                self.perc_del = self.query_gap_pct
        else:
            raise ValueError(
                "Error: 'reference' kwarg must be 'query', 'target', or omitted (default is query)"
            )

    # Method used by initializer to make sure sequences are same length and all legal chars
    def sequenceFilter(self, querySeq, targetSeq):
        if querySeq == targetSeq == None:
            return
        elif querySeq == None:
            raise ValueError(
                "Error: querySeq was given, but not Target - must give both or neither seq"
            )
        elif targetSeq == None:
            raise ValueError(
                "Error: targetSeq was given, but not query - must give both or neither seq"
            )
        if len(querySeq) != len(targetSeq):
            raise ValueError(
                "query Seq and Target Seq are different lengths - must be same lenght"
            )
        allowedChars = r"[^\-a-zA-Z.]"
        if re.search(allowedChars, querySeq):
            raise ValueError(
                "query Seq contains illegal characters (only alphabet and '-' allowed)"
            )
        if re.search(allowedChars, targetSeq):
            raise ValueError(
                "Target Seq contains illegal characters - (only alphabet and '-' allowed)"
            )

    # function to ensure the minimum number of kwargs are submitted
    def minimumArgumentEnforcer(self, keyArgs):
        minArgs = [
            "score",
            "query_id",
            "query_start",
            "query_end",
            "target_id",
            "target_start",
            "target_end",
            "orientation",
        ]
        for val in minArgs:
            if val not in keyArgs:
                raise ValueError("Error: " + val + " is a required argument")

    # Method to enforce type fidelity of constructor keyword args
    def keyArgumentTypeEnforcer(self, keyArgs):
        allowedTypes = {
            "align_id": int,
            "score": int,
            "e_value": float,
            "p_value": float,
            "bit_score": float,
            "bias": float,
            "perc_sub": float,
            "perc_ins": float,
            "perc_del": float,
            "kimura_div": float,
            "query_id": str,
            "query_start": int,
            "query_end": int,
            "query_len": int,
            "query_name": str,
            "env_start": int,
            "env_end": int,
            "aligned_query_seq": str,
            "target_id": str,
            "target_start": int,
            "target_end": int,
            "target_len": int,
            "aligned_target_seq": str,
            "matrix_name": str,
            "pp_string": str,
            "rf_string": str,
            "orientation": str,
            "query_name": str,
            "sub_pct": float,
            "query_gap_pct": float,
            "target_gap_pct": float,
            "reference": str,
            "target_name": str,
        }
        for key in keyArgs:
            if type(keyArgs[key]) != allowedTypes[key]:
                raise ValueError(
                    "Error: " + key + " must be of type " + str(allowedTypes[key])
                )
            elif allowedTypes[key] == float or allowedTypes[key] == int:
                if keyArgs[key] < 0:
                    raise ValueError(
                        "Error: "
                        + str(key)
                        + " was less than zero: "
                        + str(keyArgs[key])
                    )
        if ("orientation" in keyArgs) and (
            keyArgs["orientation"] != "+" and keyArgs["orientation"] != "-"
        ):
            raise ValueError("Error: Orientation must be either '+' or '-'")

    #
    # High Priority Methods
    #
    def __RLE_encode(self, input_string, markup_singletons=False):
        """
        __RLE_encode() - Private method to generate a Run Length Encoding of
                         an input string.

        Run length encoding replaces tandem occurrences of a symbol in a string
        with a count of tandem occurrences of the symbol followed by a single
        instance of the symbol. For instance the following input string
        "AAACTTAA" could strictly be encoded as "3A1C2T2A".  There are multiple
        conventions for RLE including variations that put the count after they
        symbol and only encoding larger than 1 or 2 tandem occurrences in the
        input string.  In this implementation the count appears before the
        symbol and the optionally singleton symbols may be left unalered.

        Args:

            input_string: String to encode

        Returns:
            input_string encoded using RLE.
        """
        if markup_singletons:  # borrowed from NHMMSearchResult.py
            return re.sub(
                r"(.)\1*", lambda m: str(len(m.group(0))) + m.group(1), input_string
            )
        else:
            return re.sub(
                r"(.)\1*",
                lambda m: str(len(m.group(0))) + m.group(1)
                if len(m.group(0)) > 1
                else m.group(1),
                input_string,
            )
        # NOTE: See NHMMSearchResult for RLE Decode method
        # originally put by Robert VVVV
        # self.aligned_query_seq = re.sub(r'(\d+)(\D)', lambda m: m.group(2) * int(m.group(1)), input_string)

    def TSV_encode(self, version=2):
        """
        TSV_encode() - Convert to a Tab Separated Value format

        Generate a Dfam tab-separated-value encoding of a subset of the alignment
        data.

        Dfam TSV format:

                target_id
                query_id
                query_name
                score
                e_value
                bias
                query_start
                query_end
                orientation
                target_start
                target_end
                env_start
                env_end
                target_len
                CIGAR             [version 2]
                kimura_divergence [version 2]

        Args: None
        Returns: Encoded string
        """
        if self.reference == "target":
            field_order = [
                "query_id",
                "target_id",
                "target_name",
                "bit_score",
                "e_value",
                "bias",
                "target_start",
                "target_end",
                "orientation",
                "query_start",
                "query_end",
                "env_start",
                "env_end",
                "query_len",
            ]
        else:
            field_order = [
                "target_id",
                "query_id",
                "query_name",
                "bit_score",
                "e_value",
                "bias",
                "query_start",
                "query_end",
                "orientation",
                "target_start",
                "target_end",
                "env_start",
                "env_end",
                "target_len",
            ]
        out_str = None
        out_str = "\t".join(
            "" if getattr(self, x) is None else str(getattr(self, x))
            for x in field_order
        )
        if version == 2:
            CIGAR = self.CIGAR_encode()
            kimuraOut = self.kimura_divergence(True)
            out_str += "\t" + str(CIGAR) + "\t" + str(kimuraOut)
        return out_str

    def CIGAR_encode(self, reference=None, style="singleMatchTag"):
        """
        Concise Idiosyncratic Gapped Alignment Report (CIGAR) Format:

            Code   Description
            ----   ------------------------------
             M      Alignment match or mismatch
             I      Insertion relative to the reference
             D      Deletion relative to the reference
             N      Skipped region of the reference
             S      Soft clipping
             H      Hard clipping
             P      Padding
             =      Match
             X      Mismatch

        Args:
          reference : The sequence that the edit codes are relative to.
                      May be either 'query' or 'target'.

          style     : The code group to use.
                        'singleMatchTag' : Encode using "M","I", and "D" only
                        'dualMatchTag'   : Encode using "=","X", "I" and "D" only

        Returns: CIGAR string
        """
        if self.aligned_query_seq == self.aligned_target_seq == None:
            return None

        if reference == None:  # FIXME: THIS MAY BE INCORRECT BY THE END
            if self.reference != None:
                reference = self.reference
            else:
                reference = "query"
        seq = None  # inspired by NHMMSearchResult.py
        refseq = None  # ask robert if this orientation is correct

        if reference == "query":
            seq = self.aligned_target_seq
            refseq = self.aligned_query_seq
        elif reference == "target":
            seq = self.aligned_query_seq
            refseq = self.aligned_target_seq
        else:
            raise ValueError(
                "Reference parameter must be either 'query' or 'target'.  '"
                + str(reference)
                + "' was proivided"
            )

        if style != "singleMatchTag":
            seq = seq.upper()
            refseq = refseq.upper()

        if seq is None:
            raise Exception("Search result does not contain alignment data")

        cigar = []
        for pair in zip(seq, refseq):
            if pair[1] == "-" or pair[1] == ".":
                # Insertion relative to the reference
                cigar.append("I")
            elif pair[0] == "-" or pair[0] == ".":
                # Deletion relative to the reference
                cigar.append("D")
            else:
                # Match/Mismatch
                if style == "singleMatchTag":
                    cigar.append("M")
                else:
                    if pair[0] == pair[1]:
                        cigar.append("=")
                    else:
                        cigar.append("X")

        seq = self.__RLE_encode("".join(cigar), True)
        return seq

    # TODO: Raise question with Robert about adjustments made to mismatches
    def CAF_encode(self, reference=None):
        """
        Compressed Align Format (CAF):

           Human             991050 GTAG--------TTGCCAAAAGCAGGATGGGGTGGGGGTGGTGGGAATGG 991091
                                        --------  i   i i  ----     ----     v
           L1MA7_3end#LI        711 GTAGAATGGTGGTTACCAGAGGC----TGGGG----GGTGGGGGGAATGG 752

           This format uses three symbols to compress merge the two sequences of an aligment
           into one.  The compression rate increases with the similarity of the sequences
           alinged.  The symbols are "/" for mismatch, "+" for insertions and "-" for deletions.
           If the sequences are identical a single copy of the base is emitted.  If the sequences
           differ at a position the bases are emitted as query base, forward slash, target base
           (E.g. the first transition in the above alignment would be denoted "G/A" ).  If an
           insertion occurs ( gaps in the query ) the inserted sequence is emitted flanked by "+".
           If a deletion occurs the deleted sequence is emitted flanked with "-".  Thus the
           example alignment above would be represented in CAF as:

                   GTAG+AATGGTGG+TTG/ACCAA/GAA/GGC-AGGA-TGGGG-TGGG-GGTGGT/GGGGAATGG

           This format is useful if the alignments are a mixture of various queries and
           targets.  If the same query is used for a large set of alignments it other formats
           (including A3M) offer better compression.  Also note that this is not a case
           sensitive format and may be used with derived sequences that are "soft masked"
        """
        if self.aligned_query_seq == self.aligned_target_seq == None:
            return None

        if reference == None:  # FIXME: THIS MAY BE INCORRECT BY THE END
            if self.reference != None:
                reference = self.reference
            else:
                reference = "query"
        seq = None
        refseq = None
        if reference == "query":
            refseq = self.aligned_query_seq
            seq = self.aligned_target_seq
        elif reference == "target":
            refseq = self.aligned_target_seq
            seq = self.aligned_query_seq
        else:
            raise ValueError(
                "Reference parameter must be either 'query' or 'target'.  '"
                + str(reference)
                + "' was proivided"
            )

        if seq is None:
            raise Exception("Search result does not contain alignment data")

        # seq = seq.upper() NOTE: we are now preserving case of non reference sequence
        refseq = refseq.upper()
        idx = 0
        caf = []
        inDel = 0
        inIns = 0
        for pair in zip(seq, refseq):
            if pair[1] == "-" or pair[1] == ".":
                # Insertion relative to the reference
                if not inIns:
                    caf.append("+")
                caf.append(pair[0])
                inIns = 1
            elif pair[0] == "-" or pair[0] == ".":
                # Deletion relative to the reference
                if not inDel:
                    caf.append("-")
                caf.append(pair[1])
                inDel = 1
            else:
                if inDel:
                    caf.append("-")
                    inDel = 0
                elif inIns:
                    caf.append("+")
                    inIns = 0
                if pair[0].upper() == pair[1]:
                    caf.append(pair[0])
                else:
                    caf.append(pair[1])  # TODO: changed these values - correct?
                    caf.append("/")
                    caf.append(pair[0])
        return "".join(caf)

    @staticmethod
    def CAF_decode(caf_string):
        """
        CAF_decode( caf_string ) - Populate object with values stored in CAF format.

        Args:
            caf_string  :  The CAF encoded alignment string

        Returns:
            [ reference_string, derived_string ]
        """
        inDel = False
        inInsert = False
        misMatch = False
        refSeq = []
        derivedSeq = []
        for i in range(len(caf_string)):
            if caf_string[i] == "-":
                inDel = not inDel
            elif caf_string[i] == "+":
                inInsert = not inInsert
            else:
                if inDel:
                    refSeq.append(caf_string[i].upper())
                    derivedSeq.append("-")
                elif inInsert:
                    refSeq.append("-")
                    derivedSeq.append(caf_string[i])
                elif caf_string[i] == "/":
                    misMatch = True
                else:
                    if misMatch == True:
                        derivedSeq[-1] = caf_string[i]
                        misMatch = False
                    else:
                        refSeq.append(caf_string[i].upper())
                        derivedSeq.append(caf_string[i])
        return ["".join(refSeq), "".join(derivedSeq)]

    def A2MA3M_encode(self, version=2, reference="query", use_rle=False):
        # FIXME: The use_rle will fail the test because in the example, two
        # consecutive same characters will be left (GG -> GG) while the current
        # implementation of RLE in this object will replace 2 consecutive
        # characters (i.e., GG -> 2G, and TTT -> 3T). Should we leave it as such,
        # implement an option in RLE to only replace 3 chars at a time, or make it
        # replace 2 chars by default?
        # TODO: in saying that in A3M the spacer characters "." or " " are omitted,
        # is it implied that they will be replaced with "-" or that spaces will simply
        # be smashed together? Either way, version 3 needs to be implemented
        # FIXME: There can currently be no sequences that have " " or "." because
        # the object is programmed not to allow sequences with chars other than
        # the alphabet and "-". Should we change this convention?
        """
        A2M/A3M/A3MRLE Align Formats:

           Multiple alignment representation:

           The A2M format encodes a set of sequences aligned to each other using character
           case and the "-" and "." symbols to represet gaps.  Upper case characters (matches)
           and "-" (deletions) are considered aligned columns and all sequences must have the
           same number of aligned columns.  Lowercase characters and "." (or " ") denote
           inserted sequence and spacers respectively.  A3M format differs from A2M in that the
           spacer characters "." or " " are omitted.

           Single alignment representation:

           Human             991050 GTAG--------TTGCCAAAAGCAGGATGGGGTGGGGGTGGTGGGAATGG 991091
                                        --------  i   i i  ----     ----     v
           L1MA7_3end#LI        711 GTAGAATGGTGGTTACCAGAGGC----TGGGG----GGTGGGGGGAATGG 752


           The A2M/A3M format is primarily a multiple alignment format, however the encoding
           may also benefit single alignments with additional conventions.  In A2M the
           sequence case is meaningful and prior information encoded using case will be
           discarded.  Uppercase characters and "-" characters denote aligned and deleted
           positions ( relative to the query ).  Lowercase characters encode insertion bases.
           Using this scheme the above alignment would become:

                        GTAGaatggtggTTACCAGAGGC----TGGGG----GGTGGGGGGAATGG

           This format is only meaningful if two additional pieces of information accompany
           this string: A reference to the query sequence, and the query start position.  It
           is also worth noting that the A2M format could be run-length-encoded to further
           compress gap and mononucleotide repeat sequences.  E.g the above would become:

                        GTAGaatggtggTTACCAGAGGC4-T4G4-GGT6GAATGG

           We will denote this run-length-encoded extension A3MRLE
        """
        if self.aligned_query_seq == self.aligned_target_seq == None:
            return None

        seq = None
        refseq = None
        if reference == "query":
            seq = self.aligned_target_seq
            refseq = self.aligned_query_seq
        elif reference == "target":
            seq = self.aligned_query_seq
            refseq = self.aligned_target_seq
        else:
            raise ValueError(
                "Reference parameter must be either 'query' or 'target'.  '"
                + str(reference)
                + "' was proivided"
            )

        if seq is None:
            raise Exception("Search result does not contain alignment data")

        seq = seq.upper()

        idx = 0
        for pair in zip(seq, refseq):
            if pair[1] == "-" or pair[1] == ".":
                seq = seq[:idx] + pair[0].lower() + seq[idx + 1 :]
            idx += 1

        if use_rle == True:
            seq = self.__RLE_encode(seq)
            seq += "    H e r e   I   a m  "
        return seq

    # Default (query, target)...see below
    # TODO: Tested by hand, not formally tested
    # TODO: while this appears to work, it could use some simplification / beautification
    # FIXME: This has been written to preserve case - keep as such?
    # TODO: has not been tested with show_alignment = true and reference = query

    def crossmatch_encode(
        self, show_alignment=False, out_file_format=False, reference=None
    ):

        # This section is come clerical setup
        if (
            self.aligned_query_seq == self.aligned_target_seq == None
            and show_alignment == True
        ):
            return None

        # default reference is "query" by convention if none given
        if reference == None:
            if self.reference == None:
                reference = "query"
            else:
                reference = self.reference

        ########################################################################

        # This section creates and populates dictionary of values that appear in
        # crossmatch alignment opening or "top" line summary
        topFieldVals = {
            "score": None,
            "sub_pct": None,
            "target_gap_pct": None,
            "query_gap_pct": None,
            "query_id": None,
            "query_start": None,
            "query_end": None,
            "query_left": None,
            "orientation": None,
            "target_id": None,
            "target_left": None,
            "target_start": None,
            "target_end": None,
            "align_id": None,
        }

        # Load most data into topFieldVals
        for key in topFieldVals:
            if key in self.__dict__:
                if isinstance(self.__dict__[key], float):
                    topFieldVals[key] = format(self.__dict__[key], ".2f")
                else:
                    topFieldVals[key] = self.__dict__[key]
                    # We must properly format the floats, may as well do it here!

        topFieldVals["query_left"] = "(" + str(self.query_len - self.query_end) + ")"
        topFieldVals["target_left"] = "(" + str(self.target_len - self.target_end) + ")"

        # orientation differs based on object orientation and out_file_format arg
        if topFieldVals["orientation"] == "-":
            topFieldVals["orientation"] = "C"
        else:
            if out_file_format == True:
                topFieldVals["orientation"] = "+"
            else:
                topFieldVals["orientation"] = ""

        ########################################################################
        # This section will be ordering the values based on factors like orientation,
        # object reference, and the reference argument
        # Some of these sections may need to be reversed (like TFSecondTwo Based
        # on reference, or TFBotStrandData depending on Orientation) - it is also
        # possible that the entire ordering may need to be swapped (like changing
        # which is top and bottom strand in some cases)
        TFFirstTwo = ["score", "sub_pct"]
        TFSecondTwo = ["target_gap_pct", "query_gap_pct"]
        TFTopStrand = ["query_id"]
        TFTopStrandData = ["query_start", "query_end", "query_left"]
        TFOrient = ["orientation"]
        TFBotStrand = ["target_id"]
        TFBotStrandData = ["target_start", "target_end", "target_left"]
        reverseComplement = False  # This is a handy tool that will help us later

        #           [0] / [5]      [1] / [6]      [2]           [3]            [4]
        topField = [
            TFFirstTwo,
            TFSecondTwo,
            TFTopStrand,
            TFTopStrandData,
            TFOrient,
            TFBotStrand,
            TFBotStrandData,
        ]

        # THIS IS THE ORIGINAL
        if reference == "target":
            topField[1].reverse()

        if self.reference != reference and topFieldVals["orientation"] == "C":
            reverseComplement = True
            # This is just a neat way to switch list element positions without holder
            topField[2], topField[5] = topField[5], topField[2]
            topField[3], topField[6] = topField[6], topField[3]

        if topFieldVals["orientation"] == "C":
            topField[6].reverse()

        ########################################################################

        # This section will deal with adding data for out_file_format
        # This MAY need to be adjusted away from just checking self.target_id
        # but it's working for now!

        if out_file_format == True:
            topField[5].append("query_class")
            topField[6].append("align_id")
            if ("#") in self.target_id:
                topFieldVals["target_id"] = self.target_id.split("#")[0]
                topFieldVals["query_class"] = self.target_id.split("#")[1]
            else:
                topFieldVals["query_class"] = "None"

        ########################################################################
        # Here we will construct the top summary string - show_alignment is False,
        # This summary data will be all that we return - otherwise, we process
        # alignmnet data below and return everything later
        topString = ""
        for valueSet in topField:
            for val in valueSet:
                if topFieldVals[val] != None:
                    topString += str(topFieldVals[val])
                else:
                    topString += "0"
                topString += "\t"

        if show_alignment == False:
            return topString

        ########################################################################
        # Here we will evaluate which strand is on "top" in the alignment, and which
        # is on bottom - this will affect printing and base pair coordinates
        outputStr = topString + "\n\n"

        topStrandName = topFieldVals[topField[2][0]]
        botStrandName = topFieldVals[topField[5][0]]
        topStrand = ""
        botStrand = ""
        if topStrandName == self.query_id:
            topStrand = self.aligned_query_seq
            botStrand = self.aligned_target_seq
            topCount = self.query_start
            if self.orientation == "+":
                bottomCount = self.target_start
            else:
                bottomCount = self.target_end
        else:
            topStrand = self.aligned_target_seq
            botStrand = self.aligned_query_seq
            topCount = self.target_start
            if self.orientation == "+":
                bottomCount = self.target_start
            else:
                bottomCount = self.query_end

        if reverseComplement == True:
            botStrandName = "C " + botStrandName
            topStrandName = "  " + topStrandName
            topStrand, botStrand = botStrand, topStrand
        elif self.orientation == "-":
            # FIXME: THIS DOESN'T SEEM RIGHT BUT IT'S WHAT THE TEST WANTS...
            # topStrandName = "C " + topStrandName
            botStrandName = "C " + botStrandName
            topStrandName = "  " + topStrandName
        else:
            botStrandName = "  " + botStrandName
            topStrandName = "  " + topStrandName

        topLine = topStrandName.ljust(30) + str(topCount).rjust(10) + " "
        bottomLine = botStrandName.ljust(30) + str(bottomCount).rjust(10) + " "
        annotation = ""
        midLineSpace = "".ljust(41)  # Note: writing this to preserve case - keep?

        ########################################################################
        # Here we reverse the sequences if we are dealing with reverse complements -
        # we will complement them right below
        complementBases = None
        if reverseComplement == True:
            topHolder = ""
            botHolder = ""
            complementBases = {"A": "T", "C": "G", "G": "C", "T": "A"}
            for i in range(len(topStrand)):
                if topStrand[i].upper() in complementBases:
                    topHolder = complementBases[topStrand[i].upper()] + topHolder
                else:
                    topHolder = topStrand[i] + topHolder

                if botStrand[i].upper() in complementBases:
                    botHolder = complementBases[botStrand[i].upper()] + botHolder
                else:
                    botHolder = botStrand[i] + botHolder
            topStrand = botHolder
            botStrand = topHolder

        #######################################################################
        # Here we will set up some variables we will need to get the transitions
        # and transitions, kimura, etc.
        # must keep track of these for final report and indexing
        transitions = 0
        transversions = 0
        gapTop = 0
        gapBottom = 0
        finalIndex = 0  # this value is used to know how much more to print
        gapInitCount = 0  # (since we only print every 50 chars and will )
        # (have to print a few stragglers at the end)

        # List of mutation types to help detect transitions and transversions
        mutType = {
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

        ########################################################################
        # Here we parse through the sequences (1 pass) and add to transi, transv,
        # gaptop, gapInit, and the middle (annotation) line appropriately
        # FIXME: THIS DOES NOT PRESERVE CASE YET!!!!!!!
        for i in range(len(topStrand)):
            # Here every 50 chars we construct the top, middle (annotation line)
            # and bottom of 1 'chunk' of the alignment
            # TODO: neaten this up Sean, it looks pretty horrific - Sean
            if i % 50 == 0 and i > 0:
                # This is where we get our indexes to flank the sequences with
                topRightIndex = topCount + i - gapTop - 1
                if self.orientation == "+":  # again, direction we count depends on
                    bottomRightIndex = bottomCount + i - gapBottom - 1  # orientation
                    bottomLeftIndex = bottomRightIndex + 1
                else:
                    bottomRightIndex = bottomCount - i + gapBottom + 1
                    bottomLeftIndex = bottomRightIndex - 1
                topLeftIndex = topRightIndex + 1

                # Here we construct the lines and append them to the outputStr
                # Note the string slicing with the index to get the right 50 bp
                topLine += topStrand[(i - 50) : i] + " " + str(topRightIndex)
                bottomLine += botStrand[(i - 50) : i] + " " + str(bottomRightIndex)
                outputStr = (
                    outputStr
                    + topLine
                    + "\n"
                    + midLineSpace
                    + annotation[(i - 50) : i]
                    + "\n"
                    + bottomLine
                    + "\n\n"
                )

                # Reset everything for the next pass...
                topLine = topStrandName.ljust(30) + str(topLeftIndex).rjust(10) + " "
                bottomLine = (
                    botStrandName.ljust(30) + str(bottomLeftIndex).rjust(10) + " "
                )
                midLine = "".ljust(41)
                # set this vv value so we know where to slice for the final stragglers
                finalIndex = i

            # Here is where we are actually evaluating the sequences and how they
            # match - this part must come /after/ the code to construct the output
            # or the indexing for the annotation line will get screwed up (it will
            # evaluate the 51st base pair before calculating coordinates, which may
            # ruin coordinate calculation if the base is a gap "-")
            if topStrand[i] == botStrand[i]:
                annotation += " "
            elif topStrand[i] == "-":
                annotation += "-"
                gapTop += 1
                if i == 0:
                    gapInitCount += 1
                elif i > 0 and topStrand[i - 1] != "-":
                    gapInitCount += 1
            elif botStrand[i] == "-":
                annotation += "-"
                gapBottom += 1
                if i == 0:
                    gapInitCount += 1
                elif i > 0 and botStrand[i - 1] != "-":
                    gapInitCount += 1
            else:
                mutPair = (topStrand[i] + botStrand[i]).upper()
                type = mutType[mutPair]
                if type == 1:
                    annotation += "i"
                    transitions += 1
                elif type == 2:
                    annotation += "v"
                    transversions += 1

                # Here is where we get the stragglers, or all bp if there are less
                # than 50 - Note: this is not superfluous, there will ALWAYS be
                # stragglers even if the seqence length is divisible by 50
            if i == len(topStrand) - 1:
                topRightIndex = str(topCount + i - gapTop)
                if self.orientation == "+":  # see above
                    bottomRightIndex = str(bottomCount + i - gapBottom)
                else:
                    bottomRightIndex = str(bottomCount - i + gapBottom)
                topLine += topStrand[finalIndex:] + " " + topRightIndex
                bottomLine += botStrand[finalIndex:] + " " + bottomRightIndex
                outputStr = (
                    outputStr
                    + topLine
                    + "\n"
                    + midLineSpace
                    + annotation[finalIndex:]
                    + "\n"
                    + bottomLine
                    + "\n\n"
                )

        ########################################################################
        # After the alignment data is complete, here we construct the final summary
        # using data collected during the pass
        if self.matrix_name is not None:
            matrixLine = "\nMatrix = " + self.matrix_name
        else:
            matrixLine = "\nMatrix = Unknown"
        kimuraVal = None
        if self.kimura_div != None:
            kimuraVal = self.kimura_div
        else:
            kimuraVal = self.kimura_divergence(True, reference=reference)
        kimuraLine = "\nkimura (with divCpGmod) = " + str(kimuraVal)
        if transversions != 0:
            trans = transitions / transversions
        else:
            trans = 0
        transLine = "\nTransitions / Transversion = {:.2f} ({} / {})".format(
            trans, transitions, transversions
        )
        nonGapChars = len(topStrand) - gapBottom - gapTop - 1
        gapChars = gapBottom + gapTop
        if gapChars != 0:
            gapInitRate = round((gapInitCount / gapChars), 2)
        else:
            gapInitRate = 0
        if gapInitCount != 0:
            avgGap = round((gapChars / gapInitCount), 2)
        else:
            avgGap = 0
        gapLine = "\ngap_init rate: {} ({} / {}), avg gap size = {} ({} / {})".format(
            gapInitRate, gapInitCount, nonGapChars, avgGap, gapChars, gapInitCount
        )

        # append the final summary data and ship it off!
        outputStr += matrixLine + transLine + kimuraLine + gapLine + "\n"
        # print("\n" + outputStr + "\n")
        return outputStr

    # TODO: properly test div_cpg_mod = True functionality
    def kimura_divergence(
        self, div_cpg_mod=False, reference="query"
    ):  # TODO: should we move 'wellCharacterizedBases'?
        """
         kimura_divergence() - Calculate Kimura sequence divergence

         Calculate a standard or CpG aware Kimura divergence from
         the alignment strings ( if present ).

         Args:

             divCpGMod:   Treat CpG sites specially.  At a
                      CpG site single transitions will be
                      counted as 1/10 of a transition and
                      two transitions will be counted as
                      one.  Transversions are unaffected
                      by this option and will be counted
                      normally.

        Returns:

             kimura:
             transitions:
             transversions:
             wellCharacterizedBases:
             CpGSites:
        """
        if self.aligned_query_seq == self.aligned_target_seq == None:
            return None

        wellCharacterizedBases = {
            "CG": 1,
            "CA": 1,
            "CC": 1,
            "CT": 1,
            "TA": 1,
            "TC": 1,
            "TG": 1,
            "TT": 1,
            "AA": 1,
            "AC": 1,
            "AG": 1,
            "AT": 1,
            "GA": 1,
            "GC": 1,
            "GG": 1,
            "GT": 1,
        }
        mutType = {
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
        # if (div_cpg_mod == False):
        #    raise NotImplementedError( "missing required implementation: see SearchResult.pm::calcKimuraDiv()" )
        if reference == "query":
            subjData = self.aligned_target_seq
            queryData = self.aligned_query_seq
        else:
            subjData = self.aligned_query_seq
            queryData = self.aligned_target_seq
        if div_cpg_mod == True:
            queryData = queryData.lower().replace("cg", "CG")
        transI = 0
        transV = 0
        wellChar = 0
        cpgStartMut = 0

        for i in range(0, len(queryData)):
            qBase = queryData[i]
            sBase = subjData[i]

            mutPair = (qBase + sBase).upper()
            if mutPair in wellCharacterizedBases:
                wellChar = wellChar + 1
            else:
                continue

            if qBase.upper() != sBase.upper():
                mtype = mutType[mutPair]
                if mtype == 1:
                    if qBase == "C":
                        transI = transI + 1 / 10
                        cpgStartMut = 1
                    elif qBase == "G":
                        if cpgStartMut == 1:
                            transI = transI + 9 / 10
                        else:
                            transI = transI + 1 / 10
                        cpgStartMut = 0
                    else:
                        transI = transI + 1
                else:
                    transV = transV + 1
            else:
                cpgStartMut = 0
        kimura = 100.00
        if wellChar >= 1:
            p = transI / wellChar
            q = transV / wellChar
            logOperand = 0
            try:
                logOperand = (1 - (2 * p) - q) * math.pow(1 - (2 * q), 0.5)
            except ValueError:
                LOGGER.warning(
                    "kimuraDivergence: p/q are not valid numbers: p = %d, q = %d. Sequences:",
                    p,
                    q,
                )
                LOGGER.warning("%s", queryData)
                LOGGER.warning("%s", subjData)
            if logOperand > 0:
                kimura = math.fabs((-0.5 * math.log(logOperand))) * 100
        return round(kimura, 2)

    def __repr__(self):
        """
        __repr__() - Generic representation of an object in JSON format
        """
        return json.dumps(self.__dict__, indent=4)

    #
    # LOW PRIORITY METHODS
    #
    def seed_pattern_count(self, bit_pattern):
        """
        seed_pattern_count() - Assess seed pattern matches in an existing
                               alignment.

        Run the seed pattern over the alignment and calculate the number of
        seed matches using given a seed bit pattern. This can be used to
        determine the sensitivity of a seed pattern on a given alignment.

        NOTE: The un-interrupted seed may not be visible in the final alignment
              due to gaps being inserted.  Use the ungapped sequences to
              find the potential seed matches.

        Args:

            bit_pattern    : The seed pattern to run over the alignment.
                             ie. "10101011110111" or "111111111"

        Returns:
            Count (int) of seed pattern matches

        """
        raise NotImplementedError("low-priority unimplemented method")

    def fragment_alignment(self, region_list):
        """
        fragment_alignment() - Fragment alignment into individual objects...
        """
        raise NotImplementedError("low-priority unimplemented method")

    def raw_to_bit_score(self, lam, mu):
        """
        Calculate the bitscore given using the rawscore stored in the
        searchresult and the scoring system's lambda, and mu parameters.
            = ( (score * lam) - log(mu) ) / log(2)
        """
        raise NotImplementedError("low-priority unimplemented method")

    def rescore_alignment(
        self,
        matrix,
        gap_init,
        gap_ext,
        ins_gap_ext,
        del_gap_ext,
        score_cpg_mod,
        div_cpg_mod,
        x_drop,
        complexity_adjust,
    ):
        """
          rescore_alignment() -

          Use the provided scoring system to rescore the alignment data
          stored in the object.  Does not alter objects values.  This
          routine will report the score ( complexity_adjusted if specified ),
          the Kimura divergence ( NOTE: not the same as CrossMatch reports ),
          the number of CpG sites in alignment, the percent insertions, and
          percent deletions.

          TODO: At this time, only the score is recalculated for
                complexity adjustment and returned.

          Well characterized bases is a metric of the number of positions where
          there are match/mismatch aligned bases [ACGT].  IUB codes, and gap
          portions of the alignment are not included.

          In addition to a overall score an array of per-position accumulated
          score is returned ( positionScores ).

          This routine assumes that the query sequence is the consensus and
          the target sequence is the genomic sequence.

          Args:
                   matrix:   An instance of Matrix.pm
                 gap_init:   A number (typically negative) used
                             to penalize the initiation of a gap
                             and includes the first position of
                             the gap.
                  gap_ext:   The penalty to apply to extensions
                             of a gap beyond the first position
                             and applies equally to insertions
                             and deletions.  Optionally the
                             insertion/deletion penalty may be
                             separated using the next two options.
              ins_gap_ext:   The penalty to apply to extensions of
                             an insertion beyond the first
                             position.
              del_gap_ext:   The penalty to apply to extensions of
                             an deletion beyond the first
                             position.
        complexity_adjust:   Use Phil Green's complexity adjusted
                             scoring cross_match function.
                   x_drop:   Calculate an xDrop function accross
                             the alignment and determine where the
                             high scoring subalignments are found.
              div_cpg_mod:   Treat CpG sites specially in the
                             Kimura divergence calculation.  At a
                             CpG site single transitions will be
                             counted as 1/10 of a transition and
                             two transitions will be counted as
                             one.  Transversions are unaffected
                             by this option and will be counted
                             normally.  In human transitions are
                             15x more likely in CpG sites in
                             rodents they are less likely.
            score_cpg_mod:   Only score transversions at CpG
                             sites. Use the standard matrix
                             C->C and G->G scores in place
                             of transitions.


        """

        # TODO: This code is very incompletely translated from SearchResult.pm.
        # At this time only complexity-adjustment is implemented, and only
        # returns correct results if the original search was done with the
        # basic scoring mode.

        # TODO: actually re-score the alignment first
        score = self.score

        query_frequencies = {}
        for s_base, q_base in zip(self.aligned_target_seq, self.aligned_query_seq):
            if s_base == "-":
                continue

            if q_base == "-":
                continue

            if q_base not in query_frequencies:
                query_frequencies[q_base] = 0
            query_frequencies[q_base] += 1

        if complexity_adjust:
            t_factor = 0
            t_sum = 0
            t_counts = 0
            for (ch, freq) in matrix.background_freqs.items():
                count = query_frequencies.get(ch, 0)
                if count > 0 and freq > 0 and math.log(freq) != 0:
                    t_factor += count * math.log(count)
                    t_sum += count * math.log(freq)
                    t_counts += count

            if t_counts != 0:
                t_factor -= t_counts * math.log(t_counts)
            t_sum -= t_factor

            lmbda = matrix.getLambda()

            # To more closely match the method employed in cross_match we use their
            # rounding scheme rather than rounding or string-formatting

            new_score = int(score + t_sum / lmbda + 0.999)
            new_score = max(0, new_score)
        else:
            raise NotImplementedError(
                "Only complexity adjustment is implemented so far!"
            )

        return new_score

    def nishimaki_sato_divergence(self, div_cpg_mod=False):
        """
        nishimaki_sato_divergence() - Calculate Nishimaki/Sato divergence

        A gap-aware extension to the Kimura divergence metric ( Nishimaki,
        Sato 2019 ).  NOTE: This treats each gap site as independent. The
        net effect will be that a 2 - 3bp gaps have the same divergence
        as 1 - 6bp gap.

        Args:

            div_cpg_mod:   Treat CpG sites specially.  At a
                           CpG site single transitions will be
                           counted as 1/10 of a transition and
                           two transitions will be counted as
                           one.  Transversions are unaffected
                           by this option and will be counted
                           normally.
        Returns:

            K2PGap:
            transitions:
            transversions:
            wellCharacterizedBases:
            CpGSites:
            gapLen:
        """
        raise NotImplementedError("low-priority unimplemented method")

    def __RLE_decode(self, input_str):
        """
        RLE_decode() - Generate alignment from Run Length Encoding string
        """
        raise NotImplementedError(
            "missing required implementation: see NHMMSearchResult.py"
        )
