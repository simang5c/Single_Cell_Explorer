import logging
import re
import tempfile
import os
from shutil import copyfileobj
from typing import Dict, Any, List, Optional
import scanpy as sc
import pandas as pd
from .state import solara_state, progress, cluster_id, AppState

# Set up logging
logger = logging.getLogger(__name__)

def fix_gene_name(genes: List[str]) -> List[str]:
    """Standardize gene names."""
    return [
        re.sub(r'(r|R)(i|I)(k|K)$', 'Rik', g.strip())[0].upper() + g.strip()[1:]
        if g else g for g in genes
    ]

def load_file(file: Dict[str, Any]) -> None:
    """Load and preprocess h5ad file with robust error handling."""
    logger.info("Starting file upload: %s", file["name"])
    solara_state.set(AppState(
        genes=solara_state.value.genes,
        organism=solara_state.value.organism,
        is_loading=True,
        loading_error=None
    ))
    progress.value = 0

    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as temp_file:
            copyfileobj(file["file_obj"], temp_file)
            temp_file_path = temp_file.name

        logger.info("Reading h5ad file: %s", temp_file_path)
        adata = sc.read_h5ad(temp_file_path)
        logger.info("AnnData object loaded with %d cells and %d genes", adata.n_obs, adata.n_vars)

        if "X_umap" not in adata.obsm:
            logger.error("No UMAP embeddings found in AnnData object")
            raise ValueError("No UMAP embeddings found in the uploaded file.")

        for col in adata.obs.columns:
            if col.startswith("RNA_snn_res."):
                adata.obs[col] = adata.obs[col].astype(int)
            elif col == "orig.ident":
                adata.obs[col] = adata.obs[col].astype(str)

        adata.raw = adata
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        clustering_ids = [
            col for col in adata.obs.columns
            if col in {"orig.ident", "seurat_clusters"} or col.startswith("RNA_snn_res.")
        ]
        clustering_ids = sorted(
            clustering_ids,
            key=lambda x: (
                0 if x == "orig.ident" else
                1 if x == "seurat_clusters" else
                2 + float(x.split("RNA_snn_res.")[-1]) if x.startswith("RNA_snn_res.") else 3
            )
        )
        logger.info("Available clustering columns: %s", clustering_ids)

        new_cluster_id = clustering_ids[0] if clustering_ids else "None"
        logger.info("Set cluster_id to: %s", new_cluster_id)

        df = pd.DataFrame(
            adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index
        ).assign(Cluster=adata.obs[new_cluster_id].astype(str))

        solara_state.set(AppState(
            genes=solara_state.value.genes,
            organism=solara_state.value.organism,
            adata=adata,
            df=df,
            clustering_ids=clustering_ids,
            is_loading=False,
            loading_error=None
        ))
        cluster_id.value = new_cluster_id
        progress.value = 100
        logger.info("Loaded .h5ad file successfully")
    except Exception as e:
        logger.error("Error loading file: %s", str(e))
        solara_state.set(AppState(
            genes=solara_state.value.genes,
            organism=solara_state.value.organism,
            loading_error=f"Error loading file: {str(e)}",
            is_loading=False
        ))
        progress.value = 100
    finally:
        if 'temp_file_path' in locals():
            try:
                os.unlink(temp_file_path)
                logger.info("Temporary file deleted: %s", temp_file_path)
            except Exception as e:
                logger.error("Error deleting temp file: %s", str(e))

def update_progress(value: Optional[float]) -> None:
    """Update progress state."""
    if value is not None:
        new_progress = min(max(float(value), 0), 100)
        logger.debug("Progress callback triggered: %.1f%%", new_progress)
        progress.value = new_progress

