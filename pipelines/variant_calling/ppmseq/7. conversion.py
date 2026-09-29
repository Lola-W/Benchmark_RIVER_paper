#!/usr/bin/env python3

import gzip
import os
import sys
from argparse import ArgumentParser
import pandas as pd


def rev_triN(triN: str) -> str:
    """Reverse complement of a trinucleotide context."""
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return comp[triN[6]] + '[' + comp[triN[2]] + '>' + comp[triN[4]] + ']' + comp[triN[0]]


def vcf_to_df(vcf_path: str) -> pd.DataFrame:
    """Convert ppmSeq outMap featuremap VCF to DataFrame."""
    features = [
        'X_EDIST','X_FC1','X_FC2','X_FILTERED_COUNT','X_FLAGS',
        'X_INDEX','X_LENGTH','X_MAPQ','X_READ_COUNT',
        'X_SCORE','X_SMQ_LEFT','X_SMQ_RIGHT',
        'a3','rq','st','tm','et',
        'max_softclip_length','hmer_context_ref','hmer_context_alt',
        'next_1','next_2','next_3','prev_1','prev_2','prev_3','ML_QUAL'
    ]
    flag_features = ['is_cycle_skip', 'is_forward'] #존재 여부로 0/1 판단

    records = {k: [] for k in features + flag_features +
               ['chrom','pos','ref','alt','qual','filt','triN','prev_3bp','next_3bp']}

    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, pos, _, ref, alt, qual, filt, info = fields[:8]

            # Parse INFO into dict; flag features separately
            info_dict, row_flags = {}, {k: 0 for k in flag_features}
            for item in info.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v
                else:
                    row_flags[item] = 1

            # Build trinucleotide context
            v = info_dict.get('trinuc_context_with_alt', '')
            triN = f"{v[0]}[{v[1]}>{v[3]}]{v[2]}" if len(v) >= 4 else None
            if info_dict.get('X_FLAGS') in ('16', '1040') and triN:
                triN = rev_triN(triN)

            # prev/next 3bp context
            p = [info_dict.get(f'prev_{i}') for i in (3, 2, 1)]
            n = [info_dict.get(f'next_{i}') for i in (1, 2, 3)]

            records['chrom'].append(chrom)
            records['pos'].append(int(pos))
            records['ref'].append(ref)
            records['alt'].append(alt)
            records['qual'].append(float(qual))
            records['filt'].append(filt)
            records['triN'].append(triN)
            records['prev_3bp'].append(''.join(p) if all(p) else None)
            records['next_3bp'].append(''.join(n) if all(n) else None)
            for k in features:
                records[k].append(info_dict.get(k))
            for k in flag_features:
                records[k].append(row_flags.get(k, 0))

    df = pd.DataFrame(records)

    int_cols = ['X_EDIST','X_FC1','X_FC2','X_FILTERED_COUNT','X_FLAGS','X_INDEX',
                'X_LENGTH','X_MAPQ','X_READ_COUNT','a3','max_softclip_length',
                'hmer_context_ref','hmer_context_alt','is_cycle_skip','is_forward']
    float_cols = ['X_SCORE','X_SMQ_LEFT','X_SMQ_RIGHT','rq','ML_QUAL']

    for c in int_cols:
        df[c] = pd.to_numeric(df[c], errors='coerce').astype('Int32')
    for c in float_cols:
        df[c] = pd.to_numeric(df[c], errors='coerce').astype('float32')

    return df


def main():
    ap = ArgumentParser(description="Convert ppmSeq outMap featuremap VCF to TSV DataFrame.")
    ap.add_argument("--vcf", required=True, help="Input outMap .vcf.gz")
    ap.add_argument("--out", required=True, help="Output TSV path")
    args = ap.parse_args()

    sys.stderr.write(f"[INFO] Converting VCF to TSV: {args.vcf}\n")
    df = vcf_to_df(args.vcf)
    df.to_csv(args.out, sep="\t", index=False)
    sys.stderr.write(f"[INFO] Wrote: {args.out}\n")


if __name__ == "__main__":
    main()
