from google.cloud import storage
from pathlib import Path
from typing import List, Dict
import scanpy as sc
import anndata as ad
from anndata import AnnData
import pandas as pd
import requests
import tarfile


def download_and_extract_file(my_file, base_dir="data"):
    """
    Downloads and extracts a .tar.gz matrix file from IGVF for a given FileAccession.
    """
    url = f"https://api.data.igvf.org/matrix-files/{my_file}/@@download/{my_file}.tar.gz"
    output_dir = Path(base_dir) / my_file
    output_dir.mkdir(parents=True, exist_ok=True)
    tar_path = output_dir / f"{my_file}.tar.gz"

    # Download
    response = requests.get(url)
    with open(tar_path, 'wb') as f:
        f.write(response.content)

    # Extract
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=output_dir)

    return output_dir


def download_gs_file(gs_path: str, output_folder: str = "../data") -> Path:
    """
    Download a files from a gs link using Google Cloud Storage library.
    """
    client = storage.Client()
    bucket_name, blob_name = gs_path.split('/', 1)
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    output_path = Path(output_folder) / f"{blob_name.split('/')[-1]}"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    blob.download_to_filename(output_path)

    return output_path


def download_h5ads(h5ads_gs: List, output_folder: str = "../data") -> Dict:
    """
    Download h5ad files from a list of gs paths.
    Args:
        h5ads_gs (List): List of gs paths to h5ad files.
        output_folder (str, optional): Directory to save the downloaded files. Defaults to "../data".

    Returns:
        Dict: sample_id -> Path
    """
    h5ad_map = {}
    for h5ad_gs in h5ads_gs:
        h5ad_path = Path(output_folder) / h5ad_gs.split('/')[-1]
        if not h5ad_path.exists():
            h5ad_path = download_gs_file(h5ad_gs, output_folder)
        h5ad_map[h5ad_path.stem] = h5ad_path

    return h5ad_map


def build_adata(h5ads_map: Dict) -> AnnData:
    """
    Builds an AnnData object from a previously downloaded and extracted IGVF tar.gz file.
    Assumes files live under data/{my_file}/ and uses global paths for metadata CSVs.
    """

    adatas = {}

    for sample_id, filename in h5ads_map.items():
        try:
            sample_adata = sc.read(filename)
        except FileNotFoundError:
            print(f"File not found: {filename}")
            continue
        sample_adata.var_names_make_unique()
        adatas[sample_id] = sample_adata

    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()

    return adata


def add_metadata_to_adata(adata: AnnData, assay_metadata: str, variable_to_map: str = "FileAccession", assay="snRNA-seq", barcode_map_file: str = None) -> AnnData:
    # Merge common assay metadata
    meta_df = pd.read_csv(assay_metadata)

    # Ensure variable_to_map exists in adata.obs
    if variable_to_map not in adata.obs.columns:
        raise ValueError(f"'{variable_to_map}' column not found in adata.obs. Available columns: {list(adata.obs.columns)}")

    # Ensure variable_to_map exists in meta_df
    if variable_to_map not in meta_df.columns:
        raise ValueError(f"'{variable_to_map}' column not found in assay metadata. Available columns: {list(meta_df.columns)}")

    barcode_map = None
    if barcode_map_file:
        barcode_map = pd.read_csv(barcode_map_file)

    common_fields = [
        'FileAccession', 'Institution', 'Assay', 'Preferred_Assay_Title',
        'AnalysisSet', 'MeasurementSet', 'MultiplexedSample', 'Subpool'
    ]
    meta_df = meta_df[meta_df['Assay'] == assay]

    adata.obs = adata.obs.merge(meta_df[common_fields], on=variable_to_map, how='left')

    # Determine assay & institution
    assay_title = adata.obs['Preferred_Assay_Title'].values[0]
    institution = adata.obs['Institution'].values[0]

    # Sample-level metadata for 10x/Broad SHARE
    sample_fields = [
        variable_to_map, 'Donor', 'Sample', 'lab_sample_id', 'Mouse_ID', 'Tissue',
        'Genotype', 'SampleType', 'Sex', 'Age_weeks', 'DOB', 'Age_days', 'Body_weight_g',
        'Estrus_cycle', 'Dissection_date', 'Dissection_time', 'ZT', 'Dissector', 'Tissue_weight_mg'
    ]
    if assay_title == "10x Multiome" or (assay_title == "SHARE-seq" and institution == "Broad"):
        adata.obs = adata.obs.merge(meta_df[sample_fields], on=variable_to_map, how='left')
        # Splitting the barcode rounds.
        adata.obs['bc1'] = adata.obs["cell_barcode"].str[:8]
        adata.obs['bc2'] = adata.obs["cell_barcode"].str[8:16]
        adata.obs['bc3'] = adata.obs["cell_barcode"].str[16:24]
        # Map barcodes to plate.
        if barcode_map is not None:
            adata.obs = adata.obs.merge(barcode_map, left_on='bc1', right_on='barcode', how='left')

    # Barcode-based metadata for Parse/UCI SHARE
    if assay_title == "Parse Split-seq" or (assay_title == "SHARE-seq" and institution == "UCI"):
        adata.obs['bc1'] = adata.obs['cell_barcode'].str[-8:]
        adata.obs['bc2'] = adata.obs['cell_barcode'].str[8:16]
        adata.obs['bc3'] = adata.obs['cell_barcode'].str[:8]

        if assay_title == "Parse Split-seq" and barcode_map is not None:
            adata.obs = adata.obs.merge(barcode_map, left_on='bc1', right_on='barcode', how='left')

        elif assay_title == "SHARE-seq" and institution == "UCI" and barcode_map is not None:
            adata.obs = adata.obs.merge(barcode_map, left_on='bc1', right_on='barcode', how='left')
            adata = adata[~adata.obs['well'].isna()].copy()

    return adata


def add_gene_names_to_adata(adata: AnnData, gene_metadata: str) -> AnnData:
    """
    Add gene names to an AnnData object from a gene metadata TSV.
    """
    gene_df = pd.read_csv(gene_metadata, header=None, sep="\t", names=['transcript_id', 'gene_id', 'gene_name', 'gene_name_unique', 'chromosome', 'start', 'end', 'strand'])
    # Keep only gene_id and gene_name and remove duplicates
    gene_df = gene_df[['gene_id', 'gene_name']].drop_duplicates()
    gene_df = gene_df.set_index('gene_id')
    adata.var['gene_name'] = adata.var.index.map(gene_df['gene_name'])
    return adata


def run_scrublet(adata):
    """
    Runs Scrublet doublet detection on an AnnData object.
    For multiplexed assays (Parse Split-seq or UCI SHARE-seq), it runs per well.
    For other assays, it runs on the entire dataset.
    Wells with fewer than 30 cells are skipped and retain default values.
    """

    # Initialize scrublet columns
    adata.obs['doublet_score'] = 0.0
    adata.obs['predicted_doublet'] = False

    # Multiplexed assays: run per well
    try:
        for well, subset in adata.obs.groupby(['well', 'AnalysisSet']):
            if len(subset) <= 30:
                print(f"Skipping well '{well}' (only {len(subset)} cells)")
                continue

            well_subset = adata[subset.index].copy()
            try:
                sc.pp.scrublet(well_subset, n_prin_comps=30)
            except Exception:
                print(f"Error running Scrublet on well '{well}' with {len(well_subset)} cells.")
                continue

            adata.obs.loc[subset.index, 'doublet_score'] = well_subset.obs['doublet_score'].astype(float)
            adata.obs.loc[subset.index, 'predicted_doublet'] = well_subset.obs['predicted_doublet'].astype(bool)
    except KeyError:
        print("Missing 'well' column in obs for multiplexed assay.\n Running on entire dataset.")
        sc.pp.scrublet(adata, n_prin_comps=30)

    return adata
