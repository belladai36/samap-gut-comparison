"""Build an annotated uninfected adult-zebrafish AnnData object."""
import argparse, json
from pathlib import Path
import pandas as pd
import scanpy as sc

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--matrix-dir',required=True); parser.add_argument('--cluster-json',required=True); parser.add_argument('--output',required=True); args=parser.parse_args()
    matrix_dir=Path(args.matrix_dir)
    adata=sc.read_10x_mtx(matrix_dir, var_names='gene_ids', make_unique=True, prefix='GSM7184372_PBS_')
    cluster=json.loads(Path(args.cluster_json).read_text())['data']
    annotations=pd.Series(cluster['annotations'], index=pd.Index(cluster['cells'], name='barcode'), name='cell_type')
    retained=adata.obs_names.intersection(annotations.index)
    if len(retained)!=len(annotations): raise ValueError(f'Only {len(retained)} of {len(annotations)} published annotated barcodes matched the count matrix')
    adata=adata[retained].copy(); adata.obs['cell_type']=annotations.loc[adata.obs_names].astype('category'); adata.obs['sample_id']='GSM7184372_PBS'; adata.obs['condition']='uninfected'; adata.obs['species']='Danio rerio'; adata.obs['dataset_accession']='GSE230044'
    adata.var['gene_symbol']=adata.var['gene_symbols'].astype(str); adata.var_names.name='ensembl_gene_id'; adata.layers['counts']=adata.X.copy()
    output=Path(args.output); output.parent.mkdir(parents=True,exist_ok=True); adata.write_h5ad(output,compression='gzip')
    print(adata); print(adata.obs['cell_type'].value_counts().to_string())

if __name__=='__main__': main()

