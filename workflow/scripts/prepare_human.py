"""Build a healthy-adult human gut AnnData object from the atlas H5AD."""
import argparse
from pathlib import Path
import scanpy as sc

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--input',required=True); parser.add_argument('--output',required=True); args=parser.parse_args()
    source=sc.read_h5ad(args.input, backed='r')
    mask=(source.obs['Diagnosis'].astype(str).to_numpy()=='Healthy adult') & (source.obs['Age_group'].astype(str).to_numpy()=='Adult')
    view=source[mask].to_memory(); source.file.close()
    view.obs['cell_type']=view.obs['annotation'].astype(str).astype('category'); view.obs['sample_id']=view.obs['Sample name'].astype(str).astype('category'); view.obs['region']=view.obs['Region'].astype(str).astype('category'); view.obs['species']='Homo sapiens'; view.obs['dataset_accession']='Space-Time Gut Cell Atlas'; view.obs['condition']='healthy_adult'
    view.var['gene_symbol']=view.var_names.astype(str); view.var_names=view.var['gene_ids'].astype(str); view.var_names.name='ensembl_gene_id'; view.var_names_make_unique(); view.layers['counts']=view.X.copy()
    output=Path(args.output); output.parent.mkdir(parents=True,exist_ok=True); view.write_h5ad(output,compression='gzip'); print(view); print(view.obs['cell_type'].value_counts().head(30).to_string())

if __name__=='__main__': main()

