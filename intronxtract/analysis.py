"""Analysis of the collected introns."""

import logging
import random
from collections import Counter

logger = logging.getLogger("intronxtract")

def _intron_key(intron):
    info = intron["intron_info"]
    return (info["reference_name"], info["intron_start_ref"], info["intron_end_ref"])

def count_transcriptomic_support(introns):
    """
    Annotate each intron with the number of alignments supporting the same
    intron (same contig, start and end).
    """
    support_count = Counter(_intron_key(intron) for intron in introns)
    for intron in introns:
        intron["transcriptomic_support"] = support_count[_intron_key(intron)]
    return introns

def _choose_criterion(introns):
    """
    Decide how to pick the representative of redundant introns:
    'identity' if any read has an identity, otherwise 'alignment_score'
    if any read has an AS tag, otherwise 'random'.
    """
    if any(i["intron_info"]["identity"] is not None for i in introns):
        return "identity"
    if any(i["intron_info"]["alignment_score"] is not None for i in introns):
        logger.warning("Read identity not available (no identity tag or NM tag): redundant introns "
                       "will be resolved keeping the read with the best alignment score (AS).")
        return "alignment_score"
    logger.warning("Neither read identity (identity tag / NM) nor alignment score (AS) available: "
                   "redundant introns will be resolved by choosing one read at random.")
    return "random"

def remove_redundant_introns(introns, seed=0):
    """
    Keep one record per intron (same contig, start and end).
    The kept intron is the one from the read with the highest identity. If no
    read has an identity, the highest alignment score (AS) is used; if there is
    no AS either, one read is chosen at random.
    """
    criterion = _choose_criterion(introns) if introns else "identity"
    rng = random.Random(seed)

    def score(intron):
        value = intron["intron_info"][criterion]
        return float("-inf") if value is None else value

    unique_introns = {}
    seen = Counter()
    removed_introns = []
    for intron in introns:
        key = _intron_key(intron)
        seen[key] += 1
        if key not in unique_introns:
            unique_introns[key] = intron
            continue

        if criterion == "random":
            # Reservoir sampling: each read has the same chance of being kept
            replace = rng.random() < 1 / seen[key]
        else:
            replace = score(intron) > score(unique_introns[key])

        if replace:
            removed_introns.append(unique_introns[key])
            unique_introns[key] = intron
        else:
            removed_introns.append(intron)

    return list(unique_introns.values()), removed_introns