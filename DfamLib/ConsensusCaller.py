# -*- coding: utf-8 -*-
"""
SEE ALSO: Dfam: http://www.dfam.org

AUTHOR(S):
    Robert Hubley <rhubley@systemsbiology.org>
    Jeb Rosen <jrosen@systemsbiology.org>

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

import numpy as np

# NB: this function is defined outside of the MultAlign class so that it is easier
# to test and add features (such as scoring matrices) separately from changes
# elsewhere in MultAlign.
def _call_consensus_ndarray(
    array, parameter_set, cpg_adjustment=True, only_positions=None
):
    """
    Call the consensus of a 2d ndarray of sequence data.

    This consensus algorithm was originally developed by Dr. Arian Smit for use
    in the development of Transposable Element families in predominantly
    mammalian genomes. The caller differs from a majority-rule consensus in two
    important ways:

      - The calls are made by scoring each column of the multiple alignment
        against a matrix, choosing the highest scoring base.  The matrix is
        a neutral-evolving DNA matrix developed from mammalian genomes which
        typically have a strong A/T bias.

      - After the consensus is called an attempt is made to restore CpG sites
        where there is evidence of C-deamination and its byproducts.

    Parameters:
        array           : Numpy 2D sequence array.

        parameter_set   : A string describing the set of parameters to use.
                          Allowed values:
                            * TODO: add support for SubstitutionMatrix instances
                            * "Linup" - the scoring matrix from Linup. This matrix is tuned for
                              the A/T bias in mammals and has support for the CpG adjustment.
                            * TODO: add more predefined options

        cpg_adjustment  : A flag controlling use of the CpG adjustment described above.
                          If the chosen parameter_set does not support the CpG adjustment,
                          an error will be raised.
                          Default: True

        only_positions  : An array of indices that should be called.
                          Default: None (all positions).

    Returns:
        A 2-element tuple (consensus, consensus_positions). Columns can be
        "missing" from consensus and consensus_positions if they would have
        been called as a gap.

    Examples:
        _call_consensus_ndarray(sequences, "Linup", cpg_adjustment=True)
        _call_consensus_ndarray(sequences, "Linup", cpg_adjustment=True, only_positions=[1,2,7])
    """

    if parameter_set == "Linup":

        # For mammals where there is a strong A/T bias
        #      A    R    G    C    Y    T    K    M    S   W   N   X   Z   -
        matrix = [
            [9, 0, -8, -15, -16, -17, -13, -3, -11, -4, -2, -7, -3, -6],
            [2, 1, 1, -15, -15, -16, -7, -6, -6, -7, -2, -7, -3, -6],
            [-4, 3, 10, -14, -14, -15, -2, -9, -2, -9, -2, -7, -3, -6],
            [-15, -14, -14, 10, 3, -4, -9, -2, -2, -9, -2, -7, -3, -6],
            [-16, -15, -15, 1, 1, 2, -6, -7, -6, -7, -2, -7, -3, -6],
            [-17, -16, -15, -8, 0, 9, -3, -13, -11, -4, -2, -7, -3, -6],
            [-11, -6, -2, -11, -7, -3, -2, -11, -6, -7, -2, -7, -3, -6],
            [-3, -7, -11, -2, -6, -11, -11, -2, -6, -7, -2, -7, -3, -6],
            [-9, -5, -2, -2, -5, -9, -5, -5, -2, -9, -2, -7, -3, -6],
            [-4, -8, -11, -11, -8, -4, -8, -8, -11, -4, -2, -7, -3, -6],
            [-2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -1, -7, -3, -6],
            [-7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -3, -6],
            [-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -6],
            [-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 3],
        ]
        alph_r = ["A", "R", "G", "C", "Y", "T", "K", "M", "S", "W", "N", "X", "Z", "-"]
        alph_h = {
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

        treat_as_n = ["B", "D", "H", "V"]

        # Consensus calling parameters
        #
        # This is what boosts the value of the CG consensus comparison if
        # the instance dinucs are CA or TG.
        ##
        ## Method better at younger families.  This table is based
        ## on the fixed matrix above.
        ##
        ##  CpG
        ##  Hypothesis  Obs   Score   Description
        ##  -----------------------------------------------------------
        ##  1           TG      +12    Direct result of current strand
        ##                               CpG deaminating the C and converting
        ##                               to a TG.
        ##  2           CA      +12    Indirect result of CpG on opp strand
        ##                               converting to TG and an incorrect
        ##                               repair of the current strand.
        ##  3           TA      -5     Two step result.  CpG -> TG followed
        ##                               by a common transition of either
        ##                               forward strand TG->TA or reverse
        ##                               strand CA->TA.
        ##  4           TT      -13    Normal transition C->T followed by
        ##                               transversion.  Scored as +2 for
        ##                               transition and matrix value for
        ##                               transversion: Matrix G->T = -15
        ##                               [+2] = -13. Could have started as
        ##                               a CG->TG->TT
        ##  5           TC      -13    (dito) Matrix G->C = -15 = -13
        ##  6           AA      -13    (dito)Matrix C->A = -15 = -13
        ##  7           GA      -13    (dito)Matrix C->G = -15 = -13
        ##  8           GT      -30    Normal matrix score C->G = -15,
        ##                               G->T = -15, Total = -30
        ##  9           AG      -5     Normal matrix score C->A = -15,
        ##                               G->G = 10, Total = -5

        # TG or CA match is 19, TG <-> CA mismatch -12. Previously set at
        # 14. Seems to overestimate in very old elements.
        cg_param = 12

        # TG or CA to TA is -4, so slightly worse than that that could have
        # arisen from CpG site.
        ta_param = -5

        # Adjust scores of transitions/transversion pairs that could have
        # arisen from CpG site.
        cg_trans_param = 2

    else:
        raise NotImplementedError(
            "The only supported parameter_set at this time is: 'Linup' (got: {})".format(
                parameter_set
            )
        )

        # TODO (for implementing custom parameter_sets): check that CpG score parameters
        # are actually available, if cpg_adjustmnet=True was requested

    array = np.char.upper(array)
    seq_cnt = array.shape[0]

    if only_positions:
        positions = list(only_positions)
    else:
        positions = list(range(array.shape[1]))

    # Build an initial consensus. Each column's value in the consensus is
    # set to whichever letter has the best sum-of-scores against all other
    # letters in the column
    consensus = ""

    # Also track which position in the alignment generated each consensus base
    cons_positions = []

    for col in positions:
        max_base = None
        max_score = None
        n_score = None

        # Tally the counts in this column
        bases, counts = np.unique(array[:, col], return_counts=True)

        # We do not want to call an "N" if there is are not any aligned
        # sequences in the column.  For instance:
        #     '    AGGA'
        #     '   AATGA'
        #     '   AAGGA'  The first three columns should not be called 'N's
        if not only_positions and bases[0] == " " and counts[0] == seq_cnt:
            continue

        # Try each possible base for this consensus column...
        for (cons_base_idx, cons_base) in enumerate(alph_r):
            score = 0

            for obs_base, obs_cnt in zip(bases, counts):
                # Don't consider ' ' in the score. ' ' represents data that is
                # "outside" the alignment at this position, and should not
                # participate positively or negatively in the score calculation
                if obs_base == " ":
                    continue

                if obs_base in treat_as_n:
                    obs_base = "N"

                obs_idx = alph_h[obs_base]

                # TODO: Track and account for different penalty scores for gap
                # insertions vs deletions and opens vs extensions (as included
                # in the SubstitutionMatrix class)
                score += obs_cnt * matrix[cons_base_idx][obs_idx]

            if cons_base == "N":
                n_score = score

            # Keep the best non-N score, breaking ties in favor of A/C/G/T
            if (
                max_score is None
                or score > max_score
                or (score == max_score and cons_base in "ACGT")
            ):
                max_base = cons_base
                max_score = score

        # Break ties in favor of N
        if n_score is not None and n_score == max_score:
            max_base = "N"

        # TODO: If a column is entirely unoccupied (" "), it gets called as N

        # By only calling non-gaps, we "delete" rare insertions from the consensus
        if max_base is not None and max_base != "-":
            consensus += max_base
            cons_positions.append(col)

    # Now, go through the consensus and consider changing each dinucleotide to a 'CG'
    # if the cpg_adjustment was requested
    if cpg_adjustment:
        adjusted_cons = list(consensus)

        for i in range(0, len(consensus) - 1):
            # Gather di-nucleotide pair and set cons_lft/cons_rgt accordingly.
            cons_lft = consensus[i]

            # Don't try to change -X to CG
            if cons_lft == "-":
                continue

            # Gaps between pair in consensus are ok: CG, C---G, C-G are all considered
            # i2 is the "next consensus position that's not a gap"
            i2 = i
            while True:
                i2 += 1
                if i2 < len(consensus):
                    cons_rgt = consensus[i2]
                    if cons_rgt != "-":
                        allgaps = False
                        break
                else:
                    allgaps = True
                    break

            if allgaps:
                break

            cons_lft_idx = alph_h[cons_lft]
            cons_rgt_idx = alph_h[cons_rgt]

            # Get the positions in the original array. These are used to find
            # the corresponding dinucleotide in each of the other sequences.
            lft_col_idx = cons_positions[i]
            rgt_col_idx = cons_positions[i2]

            dn_score = 0
            cg_score = 0
            for row_idx in range(array.shape[0]):
                seq_lft = array[row_idx, lft_col_idx]
                if seq_lft in treat_as_n:
                    seq_lft = "N"
                seq_rgt = array[row_idx, rgt_col_idx]
                if seq_rgt in treat_as_n:
                    seq_rgt = "N"

                if seq_lft in "- " or seq_rgt in "- ":
                    continue

                seq_lft_idx = alph_h[seq_lft]
                seq_rgt_idx = alph_h[seq_rgt]

                seq_dinucl = seq_lft + seq_rgt
                dn_score += matrix[cons_lft_idx][seq_lft_idx]
                dn_score += matrix[cons_rgt_idx][seq_rgt_idx]

                if seq_dinucl in ["CA", "TG"]:
                    # CpG -> TG (on either strand)
                    cg_score += cg_param
                elif seq_dinucl == "TA":
                    # CpG -> TG, repaired incorrectly to TA
                    cg_score += ta_param
                elif seq_dinucl in ["TC", "TT"]:
                    # CpG -> TG transition, followed by G->[C or T] transversion
                    cg_score += cg_trans_param + matrix[alph_h["G"]][seq_rgt_idx]
                elif seq_dinucl in ["AA", "GA"]:
                    # CpG -> TG transition on opposite strand (CA), followed by C->[A or G] transversion
                    cg_score += cg_trans_param + matrix[alph_h["C"]][seq_lft_idx]
                else:
                    cg_score += matrix[alph_h["C"]][seq_lft_idx]
                    cg_score += matrix[alph_h["G"]][seq_rgt_idx]

            if cg_score > dn_score:
                adjusted_cons[i] = "C"
                adjusted_cons[i2] = "G"

        consensus = "".join(adjusted_cons)

    return consensus, cons_positions
