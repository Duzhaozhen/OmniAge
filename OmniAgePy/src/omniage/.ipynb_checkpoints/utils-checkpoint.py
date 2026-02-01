import os
import pandas as pd
import numpy as np
import warnings
from typing import List, Union
import requests
import tarfile
from io import BytesIO
from pathlib import Path
import shutil

def get_data_path(filename: str) -> str:
    """
    Retrieves the absolute path for a file within the package's 'data' directory.
    
    This utility ensures consistent file access regardless of the execution context.
    
    Args:
        filename (str): Name of the file (e.g., 'Horvath2013.csv').
        
    Returns:
        str: Normalized absolute path to the file.
    """
    # Get the directory where utils.py resides (.../src/omniage)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(current_dir, 'data', filename)
    return os.path.normpath(data_path)

def load_clock_coefs(clock_name: str) -> pd.DataFrame:
    """
    Loads coefficient data for a specific clock from the internal data store.WW
    
    Args:
        clock_name (str): The name of the clock (file is expected to be named '{clock_name}.csv').
        
    Returns:
        pd.DataFrame: Loaded coefficients.
        
    Raises:
        FileNotFoundError: If the coefficient file does not exist.
    """
    filename = f"{clock_name}.csv"
    path = get_data_path(filename)
    
    if not os.path.exists(path):
        raise FileNotFoundError(f"Coefficient file not found: {path}")
        
    return pd.read_csv(path)



def check_and_download_zenodo(target_dir: Path, download_url: str):
    """
    Validates the existence of data in `target_dir`. If missing, automates the download 
    and extraction from Zenodo. Uses streaming to handle large datasets efficiently.

    Args:
        target_dir (Path): Local directory where the data should be located.
        download_url (str): Direct download URL for the Zenodo archive.
    """
    # 1. Check if the directory exists and contains files to avoid redundant downloads
    if target_dir.exists() and any(target_dir.iterdir()):
        print(f"✅ [{target_dir.name}] Data already exists locally. Skipping download.")
        return

    print(f"📌 [{target_dir.name}] Local data not found.")
    print(f"Attempting automatic download from Zenodo...")
    print(f"Source: {download_url}")
    
    # Define a temporary path for the downloaded archive
    temp_archive = target_dir.parent / f"{target_dir.name}_temp.tar.gz"
    
    try:
        # Ensure the parent directory structure exists
        target_dir.parent.mkdir(parents=True, exist_ok=True)

        # 2. Stream the download to disk to maintain low memory footprint
        # This prevents loading several gigabytes (e.g., 4.82 GB) into RAM
        response = requests.get(download_url, stream=True, timeout=60)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        print(f"Downloading archive (Size: {total_size / (1024**3):.2f} GB)...")
        
        with open(temp_archive, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024 * 1024):  # 1MB chunks
                if chunk:
                    f.write(chunk)
        
        # 3. Extract the downloaded archive
        print("Extracting files... (This may take a few minutes for large datasets)")
        with tarfile.open(temp_archive, mode="r:gz") as tar:
            # Extraction is performed in the parent directory, assuming the 
            # archive contains a top-level folder matching target_dir.name.
            tar.extractall(path=target_dir.parent)
            
        print(f"✅ Successfully installed data to: {target_dir}")
        
    except Exception as e:
        # --- Error Handling and Cleanup ---
        # Remove the target directory if it was created but remains empty/corrupted
        if target_dir.exists() and not any(target_dir.iterdir()):
            try:
                os.rmdir(target_dir)
            except OSError:
                pass

        abs_path = target_dir.resolve()
        
        # Determine the specific file to check for manual verification based on the dataset type
        verify_file = (
            abs_path / 'PCHorvath2013' / 'model.parquet' 
            if 'PCClocks' in str(target_dir) 
            else abs_path / 'model.parquet'
        )

        error_msg = (
            f"\n\n{'='*60}\n"
            f"❌ AUTOMATIC DOWNLOAD FAILED\n"
            f"{'='*60}\n"
            f"Reason: {str(e)}\n\n"
            f"Please install the data manually following these steps:\n"
            f"1. Download the file from this URL:\n"
            f"   {download_url}\n\n"
            f"2. Extract the downloaded .tar.gz file.\n\n"
            f"3. Place the extracted folder ('{target_dir.name}') exactly at:\n"
            f"   {abs_path}\n\n"
            f"4. Verify the existence of this file to confirm installation:\n"
            f"   {verify_file}\n"
            f"{'='*60}\n"
        )
        raise RuntimeError(error_msg) from e

    finally:
        # Clean up the temporary archive to free up disk space
        if temp_archive.exists():
            temp_archive.unlink()
# --- Example Usage ---
# check_and_download_zenodo(
#     target_dir=Path("./data/SystemsAge"), 
#     download_url="https://zenodo.org/records/18287372/files/SystemsAge.tar.gz?download=1"
# )



def PASTA_create_pasta_pseudobulks(
    adata, 
    pool_by: List[str] = ["cell_type", "age"], 
    chunk_size: int = 1000, 
    min_cells: int = 10,
    layer: str = None
) -> "AnnData":
    """
    Python implementation of the pseudobulk generation logic (akin to 'making_pseudobulks_from_seurat' in R).
    
    This function aggregates single-cell expression profiles into 'pseudobulk' samples.
    It groups cells based on specified metadata columns and randomly assigns them to 
    chunks of a target size to create robust aggregated profiles.

    Parameters
    ----------
    adata : AnnData
        Input single-cell object containing expression data and metadata.
    pool_by : list of str
        List of column names in `adata.obs` to use for grouping cells (e.g., ['cell_type', 'age']).
    chunk_size : int, default=1000
        Target number of cells to aggregate per pseudobulk sample.
    min_cells : int, default=10
        Minimum number of cells required to form a group. Groups with fewer cells are discarded.
    layer : str, optional
        Key in `adata.layers` to use for aggregation. If None, `adata.X` is used.
        Note: Pseudobulking is typically performed on raw counts, not log-normalized data.

    Returns
    -------
    AnnData
        A new AnnData object where observations (obs) represent pseudobulk samples
        and variables (var) correspond to genes from the input.
    """
    try:
        import anndata as ad
        from scipy import sparse
    except ImportError:
        raise ImportError("The 'anndata' and 'scipy' packages are required for pseudobulking.")

    # 1. Validate Grouping Columns
    for col in pool_by:
        if col not in adata.obs.columns:
            raise ValueError(f"Column '{col}' not found in adata.obs")

    # 2. Prepare Grouping Key
    obs = adata.obs.copy()
    # Create a composite key (e.g., "Astrocytes-50yr")
    obs['group_key'] = obs[pool_by].astype(str).agg('-'.join, axis=1)
    
    # 3. Assign Chunk IDs
    # Simulates the logic of splitting groups into random chunks
    obs['chunk_id'] = None 
    obs['chunk_id'] = obs['chunk_id'].astype(object)
    
    for group_name, group_df in obs.groupby('group_key'):
        n_cells = len(group_df)
        if n_cells < min_cells:
            continue
            
        # Determine number of chunks required
        n_chunks = int(np.ceil(n_cells / chunk_size))
        
        # Generate chunk indices (e.g., [0, 0, ..., 1, 1, ...])
        # Truncate to match exact cell count
        chunk_indices = np.repeat(np.arange(n_chunks), chunk_size)[:n_cells]
        
        # Shuffle to ensure random assignment within the group
        np.random.shuffle(chunk_indices)
        
        # Assign unique chunk IDs: "GroupKey-ChunkIdx" (e.g., "Astrocytes-50yr-0")
        full_chunk_ids = [f"{group_name}-{i}" for i in chunk_indices]
        obs.loc[group_df.index, 'chunk_id'] = full_chunk_ids

    # Drop cells that were filtered out (e.g., small groups)
    obs = obs.dropna(subset=['chunk_id'])
    
    if len(obs) == 0:
        warnings.warn("No cells remained after grouping. Check 'pool_by' criteria or 'min_cells' threshold.")
        return None

    # 4. Aggregation (Matrix Multiplication)
    # Uses One-Hot Encoding via get_dummies for efficient summation
    # Operation: (Chunks x Cells) * (Cells x Genes) = (Chunks x Genes)
    
    print(f"Aggregating {len(obs)} cells into pseudobulks...")
    
    # Subset AnnData to valid cells
    subset_adata = adata[obs.index]
    
    # Create aggregation matrix
    dummies = pd.get_dummies(obs['chunk_id'], dtype=float)
    chunks = dummies.columns.tolist() # Sorted chunk IDs
    agg_mat = dummies.values.T # Shape: (n_chunks, n_cells)
    
    # Convert to sparse matrix for performance
    agg_mat_sparse = sparse.csr_matrix(agg_mat)
    
    # Select source data layer
    if layer is not None:
        X_source = subset_adata.layers[layer]
    else:
        X_source = subset_adata.X
        
    # Ensure source is sparse CSR for fast dot product
    if not sparse.issparse(X_source):
        X_source = sparse.csr_matrix(X_source)
    
    # Perform Aggregation (Sum Counts)
    X_agg = agg_mat_sparse.dot(X_source)
    
    # 5. Construct Result AnnData
    new_obs = pd.DataFrame(index=chunks)
    
    # Restore metadata
    # Map chunk_id back to original pool_by attributes (constant within a chunk)
    meta_mapping = obs.drop_duplicates(subset='chunk_id').set_index('chunk_id')
    new_obs = new_obs.join(meta_mapping[pool_by])
    
    # Record cell count per pseudobulk
    new_obs['cell_count'] = dummies.sum(axis=0).values

    pb_adata = ad.AnnData(X=X_agg, obs=new_obs, var=adata.var)
    
    # Note: The output X contains summed raw counts.
    # Downstream clocks (like PASTA) typically require log-normalized data.
    # Users should run `sc.pp.normalize_total` and `sc.pp.log1p` on pb_adata if needed.
    
    print(f"Created {pb_adata.n_obs} pseudobulk samples.")
    return pb_adata

