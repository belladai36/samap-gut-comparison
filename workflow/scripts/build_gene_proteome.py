"""Collapse Ensembl peptide isoforms to the longest protein per gene."""
import argparse, gzip, re
from pathlib import Path

def records(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    header=None; seq=[]
    with opener(path,'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None: yield header,''.join(seq)
                header=line[1:].strip(); seq=[]
            else: seq.append(line.strip())
    if header is not None: yield header,''.join(seq)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--input',required=True); parser.add_argument('--output',required=True); args=parser.parse_args()
    longest={}
    for header,seq in records(args.input):
        match=re.search(r'(?:^|\s)gene:([^\s]+)',header)
        if not match: continue
        gene=match.group(1).split('.')[0]
        if gene not in longest or len(seq)>len(longest[gene]): longest[gene]=seq
    output=Path(args.output); output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w') as handle:
        for gene in sorted(longest):
            handle.write(f'>{gene}\n')
            seq=longest[gene]
            for start in range(0,len(seq),80): handle.write(seq[start:start+80]+'\n')
    print(f'Wrote {len(longest)} longest-per-gene proteins to {output}')

if __name__=='__main__': main()

