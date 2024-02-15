__all__ = [
    "sequence_identity",
    "weighted_sequence_identity",
    "sequence_identity_matrix",
    "gap_selector",
]
import numpy as np


def sequence_identity(a, b, gaps="y"):
    """Compute the sequence identity between two sequences.

    The definition of sequence_identity is ambyguous as it depends on how gaps are treated,
    here defined by the *gaps* argument. For details and examples, see
    `this page <https://pyaln.readthedocs.io/en/latest/tutorial.html#sequence-identity>`_

    Parameters
    ----------
    a : str
        first sequence, with gaps encoded as "-"

    b : str
        second sequence, with gaps encoded as "-"

    gaps : str
       defines how to take into account gaps when comparing sequences pairwise. Possible values:
       - 'y' : gaps are considered and considered mismatches. Positions that are gaps in both sequences are ignored.
       - 'n' : gaps are not considered. Positions that are gaps in either sequences compared are ignored.
       - 't' : terminal gaps are trimmed. Terminal gap positions in either sequences are ignored, others are considered as in 'y'.
       - 'a' : gaps are considered as any other character; even gap-to-gap matches are scored as identities.

    Returns
    -------
    float
        sequence identity between the two sequences

    Examples
    --------
    >>> sequence_identity('ATGCA',
    ...                   'ATGCC')
    0.8

    >>> sequence_identity('--ATC-GGG-',
                          'AAATCGGGGC',
                          gaps='y')
    0.6

    Note
    ----
    To compute sequence identity efficiently among many sequences, use :func:`~pyaln.Alignment.score_similarity` instead.

    See also
    --------
    pyaln.Alignment.score_similarity,  weighted_sequence_identity

    """

    if len(a) != len(b):
        raise IndexError(
            "sequence_identity ERROR sequences do not have the same length"
        )

    if gaps == "y":
        pos_to_remove = [i for i in range(len(a)) if a[i] == "-" and b[i] == "-"]
    elif gaps == "n":
        pos_to_remove = [i for i in range(len(a)) if a[i] == "-" or b[i] == "-"]
    elif gaps == "t":
        pos_to_remove = [i for i in range(len(a)) if a[i] == "-" and b[i] == "-"]
        for s in [a, b]:
            for i, c in enumerate(s):
                if c == "-":
                    pos_to_remove.append(i)
                else:
                    break
            for i, c in reversed(list(enumerate(s))):
                if c == "-":
                    pos_to_remove.append(i)
                else:
                    break
    elif gaps == "a":
        count_identical = sum([int(ca == b[i]) for i, ca in enumerate(a)])
        return count_identical / len(a) if len(a) else 0.0
    else:
        raise Exception(
            "sequence_identity ERROR gaps argument must be one of {a, y, n, t}"
        )

    exclude_pos = set(pos_to_remove)
    count_identical = sum(
        [
            int(ca == b[i] and ca != "-")
            for i, ca in enumerate(a)
            if not i in exclude_pos
        ]
    )
    denominator = len(a) - len(exclude_pos)
    return count_identical / denominator if denominator else 0.0


def weighted_sequence_identity(a, b, weights, gaps="y"):
    """Compute the sequence identity between two sequences, different positions differently


    The definition of sequence_identity is ambyguous as it depends on how gaps are treated,
    here defined by the *gaps* argument. For details and examples, see
    `this page <https://pyaln.readthedocs.io/en/latest/tutorial.html#sequence-identity>`_

    Parameters
    ----------
    a : str
        first sequence, with gaps encoded as "-"

    b : str
        second sequence, with gaps encoded as "-"

    weights : list of float
        list of weights. Any iterable with the same length as the two input sequences
        (including gaps) is accepted. The final score is divided by their sum
        (except for positions not considered, as defined by the gaps argument).

    gaps : str
       defines how to take into account gaps when comparing sequences pairwise. Possible values:
       - 'y' : gaps are considered and considered mismatches. Positions that are gaps in both sequences are ignored.
       - 'n' : gaps are not considered. Positions that are gaps in either sequences compared are ignored.
       - 't' : terminal gaps are trimmed. Terminal gap positions in either sequences are ignored, others are considered as in 'y'.
       - 'a' : gaps are considered as any other character; even gap-to-gap matches are scored as identities.

    Returns
    -------
    float
        sequence identity between the two sequences

    Examples
    --------
    >>> weighted_sequence_identity('ATGCA',
    ...                            'ATGCC', weights=[1, 1, 1, 1, 6])
    0.4

    >>> weighted_sequence_identity('ATGCA',
    ...                            'ATGCC', weights=[1, 1, 1, 1, 1])
    0.8

    Note
    ----
    To compute sequence identity efficiently among many sequences, use :func:`~pyaln.Alignment.score_similarity` instead.

    See also
    --------
    pyaln.Alignment.score_similarity,  weighted_sequence_identity

    """

    if len(a) != len(b):
        raise IndexError(
            "sequence_identity ERROR sequences do not have the same length"
        )

    if len(weights) != len(a):
        raise IndexError(
            "sequence_identity ERROR weights must be the same length as sequences"
        )

    if gaps == "y":
        pos_to_remove = [i for i in range(len(a)) if a[i] == "-" and b[i] == "-"]
    elif gaps == "n":
        pos_to_remove = [i for i in range(len(a)) if a[i] == "-" or b[i] == "-"]
    elif gaps == "t":
        pos_to_remove = [i for i in range(len(a)) if a[i] == "-" and b[i] == "-"]
        for s in [a, b]:
            for i, c in enumerate(s):
                if c == "-":
                    pos_to_remove.append(i)
                else:
                    break
            for i, c in reversed(list(enumerate(s))):
                if c == "-":
                    pos_to_remove.append(i)
                else:
                    break
    elif gaps == "a":
        total_weight = sum(weights)
        count_identical = sum([int(ca == b[i]) * weights[i] for i, ca in enumerate(a)])
        return count_identical / total_weight if total_weight else 0.0
    else:
        raise Exception(
            "sequence_identity ERROR gaps argument must be one of {a, y, n, t}"
        )

    exclude_pos = set(pos_to_remove)
    actual_weights = [w for i, w in enumerate(weights) if not i in exclude_pos]
    total_weight = sum(actual_weights)

    count_identical = sum(
        [
            int(ca == b[i] and ca != "-") * weights[i]
            for i, ca in enumerate(a)
            if not i in exclude_pos
        ]
    )
    return count_identical / (total_weight) if total_weight else 0.0


def gap_selector(npt, nps, gaps):
    if gaps == "a":
        selector = np.full((npt.shape[0], *nps.shape), True, dtype=bool)
    elif gaps == "n":
        ## which positions are taken into account:  those in which neither target and selfseq is gap
        selector = (npt != "-")[:, np.newaxis, :] & (nps != "-")[np.newaxis, :, :]
    elif gaps == "y":
        selector = (npt != "-")[:, np.newaxis, :] | (nps != "-")[np.newaxis, :, :]
    elif gaps == "t":
        # code from Alignment.terminal_gap_mask
        terminal_gaps_t = np.apply_along_axis(
            lambda x: (
                np.logical_and.accumulate(x) | np.logical_and.accumulate(x[::-1])[::-1]
            ),
            1,
            npt == "-",
        )
        terminal_gaps_s = np.apply_along_axis(
            lambda x: (
                np.logical_and.accumulate(x) | np.logical_and.accumulate(x[::-1])[::-1]
            ),
            1,
            nps == "-",
        )
        selector = (
            ~(terminal_gaps_t)[:, np.newaxis, :]
            & ~(terminal_gaps_s)[np.newaxis, :, :]
            & ((npt != "-")[:, np.newaxis, :] | (nps != "-")[np.newaxis, :, :])
        )
    return selector


def sequence_identity_matrix(npt, nps, selector=None, gaps=None, eq_matrix=None):
    if selector is None:
        if gaps is None:
            raise ValueError(
                "sequence_identity_matrix ERROR you must provide selector or gaps arguments"
            )
        selector = gap_selector(npt, nps, gaps)
    if eq_matrix is None:
        eq_matrix = np.char.equal(nps, npt[:, np.newaxis])

    sums = selector.sum(axis=2)
    sums[sums == 0] = 1  # avoid division 0/0
    return (
        (eq_matrix & selector).sum(axis=2)
        /
        # below: matrix of length of alignments, i.e. the positions actually used for each comparison
        sums
    )
