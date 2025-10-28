def count_transcriptomic_support(introns):
    """
    Calculates the transcriptomic support for each intron and adds it 
    to the intron dictionary.
    """
    support_count = {}

    for intron in introns:
        key = (intron['intron_info']['reference_name'],
               intron['intron_info']['intron_start_ref'],
               intron['intron_info']['intron_end_ref'])

        if key not in support_count:
            support_count[key] = 1
        else:
            support_count[key] += 1

    # assign intron count 
    for intron in introns:
        key = (intron['intron_info']['reference_name'],
               intron['intron_info']['intron_start_ref'],
               intron['intron_info']['intron_end_ref'])
        
        intron['transcriptomic_support'] = support_count[key]

    return introns

def remove_redundant_introns(introns):
    """Removes redundant introns, keeping only the one with the highest identity."""
    unique_introns = {}
    removed_introns = []

    for intron in introns:
        key = (intron['intron_info']['reference_name'],
               intron['intron_info']['intron_start_ref'],
               intron['intron_info']['intron_end_ref'])

        if key in unique_introns:
            if intron['intron_info']['identity'] > unique_introns[key]['intron_info']['identity']:
                removed_introns.append(unique_introns[key])
                unique_introns[key] = intron
        else:
            unique_introns[key] = intron

    return list(unique_introns.values()), removed_introns
