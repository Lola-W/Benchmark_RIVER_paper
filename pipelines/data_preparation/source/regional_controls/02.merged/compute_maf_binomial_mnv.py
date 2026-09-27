import sys
import re
import pandas as pd
import numpy as np
from scipy.stats import beta
from scipy import stats


def translate_bases(ref, depth, match):
    ref_lc = ref.lower()
   
   #04/22/22: base N(n) in the read pileup column

    bases = ['A','a','T','t','C','c','G','g','N','n']
    pos_n = {}
    count = {}
    for base in bases:
        count[base] = 0
    for base in bases:
        pos_n[base] = []
    I = 0; D = 0
    in_base= 0; del_base = 0; #mis_base = 0
    i = 0
    pos_q = 0
    end_reads = {}
    while(i < len(match)):
        if match[i] == '^':
            i += 2
        elif match[i] == "$":
            i += 1
        elif match[i] == ".":
            count[ref] += 1
            pos_n.setdefault(ref,[]).append(pos_q)
            pos_q += 1
            i += 1
        elif match[i] == ",":
            count[ref_lc] += 1
            pos_n.setdefault(ref_lc,[]).append(pos_q)
            pos_q += 1
            i += 1
        elif match[i] in bases:
            count[match[i]] += 1
            pos_n.setdefault(match[i],[]).append(pos_q)
            pos_q += 1
            i += 1
        elif match[i] == '+':
            I += 1
            if re.match(r"\d",match[i+2]):
                c = 10*int(match[i+1]) + int(match[i+2])
                i += c+3
                pos_q += c+3
            else:
                c = int(match[i+1])
                i += c+2
                pos_q += c+2
            in_base += c
        elif match[i] == '-':
            D += 1
            if re.match(r'\d',match[i+2]):
                c = 10*int(match[i+1])+ int(match[i+2])
                i += c+3
            else:
                c = int(match[i+1])
                i += c+2
            del_base += c
        elif match[i] == "*":
            D += 1
            i += 1
        else:
            sys.stderr.write("ERROR: " + str(i) + " " + str(match[i]) + "\n")
            sys.exit(0)

    return count, pos_n, I, D



def translate_qualities(ref, alt, pos_n, qualities): 
    #print(pos_n)
    quality_ref = ""
    quality_alt = ""
    for position in pos_n[ref]+pos_n[ref.lower()]:
        quality_ref += qualities[position]
    for position in pos_n[alt]+pos_n[alt.lower()]:
        quality_alt += qualities[position]
    quality_ref = re.sub(r'!','',quality_ref)
    quality_ref = re.sub(r'\\','',quality_ref)
    quality_alt = re.sub(r'!','',quality_alt)
    quality_alt = re.sub(r'\\','',quality_alt)
    return quality_ref, quality_alt


def compute_table(ref, alt, ref_qualities, alt_qualities):
    qualities = ref_qualities + alt_qualities
    ref_length = len(ref_qualities)
    alt_length = len(alt_qualities)
    length = ref_length + alt_length
    #print(length, ref_length)
    matrix = np.zeros((length + 1, length + 1))
    #compute matrix
    matrix[0,0] = 1
    for i in range(length):
        prob = 1- 10**-(ord(qualities[i])-33)
        if i >= ref_length:
            prob = 1-prob
        for j in range(i+1):
            matrix[j, i+1] += matrix[j, i] * prob
            matrix[j+1, i+1] += matrix[j, i] * (1-prob)
    scores = matrix[:,-1]
    return scores

def compute_MAF_and_CI(scores):
    length = len(scores)
    MAF = np.argmax(scores)/length
    scores_index = np.nonzero(scores)
    nonzero_scores = list(scores[scores_index])
    start_index = scores_index[0][0]
    if scores_index[0][0] != 0:
        start_index = start_index -1
        nonzero_scores = [0] + nonzero_scores
    if scores_index[0][-1] != len(scores) -1:
        nonzero_scores = nonzero_scores + [0]
    #expanding scores
    #print(sum(nonzero_scores))
    new_scores = []
    for i in range(len(nonzero_scores)-1):
        first = nonzero_scores[i]
        second = nonzero_scores[i+1]
        new_scores.append(first)
        new_scores.append(first + (second-first)/5)
        new_scores.append(first + (second-first)/5*2)
        new_scores.append(first + (second-first)/5*3)
        new_scores.append(first + (second-first)/5*4)
    new_scores.append(second)
    sum_scores = sum(new_scores)
    #print(sum_scores)
    #compute ci by adding up the scores from left and from right
    end_index = start_index + len(new_scores) - 1
    agg = 0
    for i in range(len(new_scores)):
        agg += new_scores[i]
        if agg >= 0.025*sum_scores:
            lower_CI = (start_index + i/5)/length
            break
    agg = 0
    for i in range(len(new_scores)):
        agg += new_scores[-(i+1)]
        if agg >= 0.025*sum_scores:
            upper_CI = (start_index+(end_index-i)/5)/length
            break
    return MAF, lower_CI, upper_CI


def clopper_binom_interval(success, total, confint=0.95):
    quantile = (1 - confint) / 2.
    lower = beta.ppf(quantile, success, total-success+1)
    upper = beta.ppf(1 - quantile, success+1, total-success)
    if np.isnan(lower):
        lower = 0
    return lower, upper

def wilson_binom_interval(success, total, alpha = 0.05):
    q_ = success / total
    crit = stats.norm.isf( alpha / 2.)
    crit2 = crit**2
    denom = 1 + crit2 / total
    center = (q_ + crit2 / (2 * total)) / denom
    dist = crit * np.sqrt(q_ * (1. - q_) / total + crit2 / (4. * total**2))
    dist /= denom
    ci_low = center - dist
    ci_upp = center + dist
    return ci_low, ci_upp
    


def parse_samtools(gt_position, gt_ref, gt_alt, file):
    gt_maf = None
    gt_lower_ci = None
    gt_upper_ci = None
    pre_upper_cis = []
    pos_upper_cis = []
    f = open(file, "r")
    gt_pos = int(gt_position.split(":")[-1])

    #04/12/25: handle the MNVs
    is_mnv = (len(gt_ref) > 1 and len(gt_alt) > 1)
    if is_mnv:
        ref_counts = []
        alt_counts = []
        for line in f:
            items = line.strip().split("\t")
            if len(items) < 6:
                continue
            chrom = items[0]
            pos = int(items[1])
            # Process only positions within the MNV span.
            if pos < gt_pos or pos >= gt_pos + len(gt_ref):
                continue
            ref_base_in_line = items[2].upper()
            if ref_base_in_line == 'N':
                continue
            depth = int(items[3])
            if depth <= 0:
                continue
            match = items[4]
            quality = items[5]
            count, pos_n, in_base, del_base = translate_bases(ref_base_in_line, depth, match)
            # Determine the position offset for the MNV
            idx = pos - gt_pos
            expected_ref = gt_ref[idx]
            expected_alt = gt_alt[idx]
            # Count both upper and lower cases.
            r_count = count.get(expected_ref, 0) + count.get(expected_ref.lower(), 0)
            a_count = count.get(expected_alt, 0) + count.get(expected_alt.lower(), 0)
            ref_counts.append(r_count)
            alt_counts.append(a_count)
        f.close()
        if len(ref_counts) == 0 or (sum(ref_counts) + sum(alt_counts)) == 0:
            return
        # Use the minimum across positions to approximate reads spanning the full haplotype.
        num_ref = min(ref_counts)
        num_alt = min(alt_counts)
        total = num_ref + num_alt
        gt_maf = num_alt / total if total != 0 else None
        gt_lower, gt_upper = wilson_binom_interval(num_alt, total, alpha=0.05)
        gt_ref_count = num_ref
        gt_alt_count = num_alt
        return gt_ref_count, gt_alt_count, gt_maf, gt_lower, gt_upper, pre_upper_cis, pos_upper_cis
    
    # Original SNP/indel handling.
    for line in f:
        #print(line.split("\t"))
        if len(line.split("\t")) > 0:
            items = line.strip().split("\t")
            chrom = items[0]
            pos = int(items[1])
            #print(pos) 
            if pos >  gt_pos and pos < gt_pos + len(gt_ref):
                continue
            position = str(chrom) + ":" + str(pos)
            #print(position)
            ref = items[2].upper()
            
            #04/22/22: skip the line when ref base is N
            
            if ref == 'N':
                continue
            depth = int(items[3])
            if depth > 0:
                match = items[4]
                quality = items[5]

                #04/22/22: handling ref base N condition
                #solution: change the ref base N to the major base from the read
                #plieup column and modify the read pileup column accordingly
                #however, need to also consider when ref base N and base N(n) in
                #the read pileup column simutaenously
                #for now, try skip when ref base N

                count, pos_n, in_base, del_base = translate_bases(ref, depth, match)
                #+- 5bp or this bp
                if position == gt_position:
                    if len(gt_ref) > len(gt_alt):
                        num_alt = del_base
                        num_ref = count[ref] + count[ref.lower()]
                    elif len(gt_ref) < len(gt_alt):
                        num_alt = in_base
                        num_ref = count[ref] + count[ref.lower()]
                    else:
                        alt = gt_alt
                        num_ref = count[ref] + count[ref.lower()]
                        num_alt = count[alt] + count[alt.lower()]
                    gt_ref_count = num_ref
                    gt_alt_count = num_alt
                    if gt_ref_count + gt_alt_count == 0:
                        return
                    else:
                        gt_maf = gt_alt_count / (gt_ref_count + gt_alt_count)
                        gt_lower, gt_upper = wilson_binom_interval(num_alt, num_alt + num_ref, alpha = 0.05)
                else:
                    allbases = {}
                    for base in ['A','T','C','G']:
                        count.setdefault(base, 0)
                        count.setdefault(base.lower(),0)
                        allbases[base] = count[base]+count[base.lower()] 
                    allbases["in"] = in_base
                    allbases["del"] = del_base
                    keys = sorted(allbases,key=lambda k:allbases[k],reverse=True)
                    if keys[0] == ref:
                        alt1 = keys[1]
                    else:
                        alt1 = keys[0]
                    num_ref = allbases[ref]
                    num_alt = allbases[alt1]                
                    ci_lower, ci_upper = wilson_binom_interval(num_alt, num_alt + num_ref, alpha = 0.05)
                    if pos < gt_pos:
                        pre_upper_cis.append(ci_upper)
                    else:
                        pos_upper_cis.append(ci_upper)
    if gt_maf is None:
        return
    f.close()
    return gt_ref_count, gt_alt_count, gt_maf, gt_lower, gt_upper, pre_upper_cis, pos_upper_cis


def main(argv):
    if len(argv) != 5:
        sys.stderr.write("usage: " + argv[0] + "<samtools_output> <chrom:pos> <REF> <ALT>\n")
        sys.exit(2)
    file = argv[1]
    gt_position = argv[2]
    gt_ref = argv[3]
    gt_alt = argv[4]
    output =  parse_samtools(gt_position, gt_ref, gt_alt, file)
    if output is None:
        gt_ref_count = np.nan
        gt_alt_count = np.nan
        gt_maf = np.nan
        gt_lower = np.nan
        gt_upper = np.nan
        greater = np.nan
    else: 
        gt_ref_count, gt_alt_count, gt_maf, gt_lower, gt_upper, pre_upper_cis, pos_upper_cis = output
        greater = ""
        for value in pre_upper_cis:
            if gt_lower <= value:
                greater += "F"
            else:
                greater += "P"
        greater += "_"
        for value in pos_upper_cis:
            if gt_lower <= value:
                greater += "F"
            else:
                greater += "P"
    #writing  
    print(gt_ref_count, gt_alt_count, gt_maf, gt_lower, gt_upper, greater)


main(sys.argv)
