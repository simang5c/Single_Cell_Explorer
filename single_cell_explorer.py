import logging
import solara
import re
import pandas as pd
from utils.state import solara_state, cluster_id, progress, AppState, ORGANISMS
from utils.utils import load_file, update_progress, fix_gene_name
from utils.plotting import create_umap_plot, create_violin_plot, create_dot_plot, create_feature_plot

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Debugging print to confirm script loading
print("Loading single_cell_explorer.py with modular structure and utils folder")

@solara.component
def Page():
    """Main Solara component for the Single Cell Explorer."""
    text = solara.use_reactive("")
    submitted = solara.use_reactive("")
    gene_error = solara.use_reactive("")

    def update_umap():
        adata = solara_state.value.adata
        cluster_col = cluster_id.value
        if adata and cluster_col in adata.obs.columns:
            logger.info("Updating UMAP with cluster: %s", cluster_col)
            solara_state.set(AppState(
                genes=solara_state.value.genes,
                organism=solara_state.value.organism,
                adata=adata,
                df=pd.DataFrame(
                    adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index
                ).assign(Cluster=adata.obs[cluster_col].astype(str)),
                clustering_ids=solara_state.value.clustering_ids,
                is_loading=solara_state.value.is_loading,
                loading_error=solara_state.value.loading_error
            ))
        else:
            logger.warning("No valid AnnData or cluster column for UMAP update")
            if solara_state.value.df is not None:
                solara_state.set(AppState(
                    genes=solara_state.value.genes,
                    organism=solara_state.value.organism,
                    adata=solara_state.value.adata,
                    df=None,
                    clustering_ids=solara_state.value.clustering_ids,
                    is_loading=solara_state.value.is_loading,
                    loading_error=solara_state.value.loading_error
                ))

    solara.use_effect(update_umap, dependencies=[solara_state.value.adata, cluster_id.value])

    def monitor_progress():
        logger.debug("Progress state changed: %.1f%%", progress.value)
    solara.use_effect(monitor_progress, dependencies=[progress.value])

    def handle_value_change(new_value: str):
        if not new_value:
            return
        logger.info("Handling gene input: %s", new_value)
        submitted.value = new_value
        text.value = ""
        gene_error.value = ""
        words = re.split(r"[,\s]+", new_value.strip())
        adata = solara_state.value.adata
        if adata:
            var_names_map = {str(name).lower(): name for name in adata.raw.var_names}
            valid_words, invalid_words = [], []
            for w in words:
                w_norm = str(w).lower()
                if w_norm in var_names_map:
                    valid_words.append(var_names_map[w_norm])
                else:
                    invalid_words.append(w)
            if not valid_words:
                gene_error.value = f"No valid genes: {words}. Available: {list(adata.raw.var_names[:5])}"
                logger.error("No valid genes entered: %s", words)
            elif invalid_words:
                gene_error.value = (
                    f"<span style='color:red'>Invalid genes: {invalid_words}</span><br>"
                    f"<span style='color:green'>Valid: {valid_words}</span>"
                )
                logger.warning("Some invalid genes: %s", invalid_words)
            solara_state.set(AppState(
                genes=[w[0].upper() + w[1:].lower() for w in valid_words] if solara_state.value.organism == "Mouse" else
                      [w.upper() for w in valid_words],
                organism=solara_state.value.organism,
                adata=solara_state.value.adata,
                df=solara_state.value.df,
                clustering_ids=solara_state.value.clustering_ids,
                is_loading=solara_state.value.is_loading,
                loading_error=solara_state.value.loading_error
            ))
        else:
            solara_state.set(AppState(
                genes=words,
                organism=solara_state.value.organism,
                adata=None,
                df=None,
                clustering_ids=[],
                loading_error="No data loaded to validate genes"
            ))

    with solara.AppBarTitle():
        solara.Text("Single Cell Explorer")

    with solara.Sidebar():
        with solara.Card(title="Upload Data", subtitle="Seurat data in h5ad format"):
            solara.FileDrop(
                label="Drop .h5ad file here",
                on_file=load_file,
                on_total_progress=update_progress
            )
            solara.ProgressLinear(value=progress.value, color="red")
            if solara_state.value.loading_error:
                solara.Error(solara_state.value.loading_error)
            if solara_state.value.is_loading:
                solara.Text("Loading file...", style={"color": "blue"})    
            solara.Text(f"Progress: {progress.value:.1f}%")
            solara.Markdown(" ")
            solara.Markdown(" ")
            solara.Text("Select Organism", style={"font-size": "18px"})
            solara.Markdown(" ")
            solara.ToggleButtonsSingle(value=solara_state.value.organism, values=ORGANISMS)
            solara.Markdown(" ")
            solara.Markdown(" ")
            solara.Text("Clustering", style={"font-size": "18px"})
            if solara_state.value.adata and solara_state.value.clustering_ids:
                solara.Select(
                    label="Select_cluster",
                    value=cluster_id,
                    values=solara_state.value.clustering_ids,
                    disabled=False
                )
            else:
                solara.Text("Please upload data to select clusters", style={"color": "gray"})
            solara.InputText(
                label="Enter genes",
                value=text,
                on_value=handle_value_change,
                update_events=["keyup.enter"],
                continuous_update=False
            )
            if submitted.value:
                genes = fix_gene_name(re.split(r'[,\s]+', submitted.value))
                solara.Text(f"Gene Markers: {genes}")
                if gene_error.value:
                    solara.Markdown(gene_error.value)
                    
    #creates a row with two columns of equal width
    with solara.Columns([1, 1], style={"gap": "10px"}):
        if solara_state.value.df is not None and not solara_state.value.df.empty:
            fig = solara.use_memo(
                lambda: create_umap_plot(solara_state.value.df, cluster_id.value, submitted.value),
                dependencies=[solara_state.value.df, cluster_id.value, submitted.value]
            )
            if fig:
                solara.FigurePlotly(fig)
            else:
                solara.Error("Failed to create UMAP plot. Check logs for details.")
        else:
            solara.Text("No data loaded for UMAP plot", style={"font-size": "18px"})

        if solara_state.value.adata:
            fig = create_violin_plot(solara_state.value.adata, solara_state.value.genes, cluster_id.value)
            if fig:
                solara.FigurePlotly(fig)
            else:
                solara.Error(f"None of the requested genes found: {solara_state.value.genes}")
        else:
            solara.Text("No data loaded for violin plot", style={"font-size": "18px"})
    
    # creates another row with two columns of equal width
    with solara.Columns([1, 1], style={"gap": "10px"}):
        if solara_state.value.adata:
            fig, error = create_dot_plot(
                id(solara_state.value.adata), tuple(solara_state.value.genes), groupby=cluster_id.value
            )
            if fig:
                solara.FigurePlotly(fig)
            else:
                solara.Error(f"Failed to create Dot plot: {error}")
        else:
            solara.Text("No data loaded for Dot plot", style={"font-size": "18px"})

        if solara_state.value.adata:
            fig = create_feature_plot(solara_state.value.adata, solara_state.value.genes, cluster_id.value)
            if fig:
                solara.FigurePlotly(fig)
            else:
                solara.Error(f"None of the requested genes found: {solara_state.value.genes}")
        else:
            solara.Text("No data loaded for Feature plot", style={"font-size": "18px"})
