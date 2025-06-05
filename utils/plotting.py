import logging
from typing import Optional, Tuple, List
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from functools import lru_cache
from .state import solara_state, cluster_id

# Set up logging
logger = logging.getLogger(__name__)

#####################################
## Function to create a UMAP plot
#####################################
def create_umap_plot(df: pd.DataFrame, cluster_id: str, submitted: str) -> Optional[go.Figure]:
    """Create UMAP scatter plot."""
    logger.info("Creating UMAP plot with cluster: %s", cluster_id)
    if df[["UMAP1", "UMAP2"]].isna().any().any():
        logger.error("UMAP1 or UMAP2 contains NaN values")
        return None
    try:
        color_column = submitted if submitted in df.columns else "Cluster"
        ordered_clusters = sorted(
            df["Cluster"].unique(),
            key=lambda x: float(x) if x.replace(".", "").isdigit() else str(x)
        )
        plot_df = df.assign(
            Cluster=pd.Categorical(df["Cluster"], categories=ordered_clusters, ordered=True)
        ).sort_values("Cluster")
        fig = px.scatter(
            plot_df, x="UMAP1", y="UMAP2", color="Cluster",
            labels={"UMAP1": "UMAP 1", "UMAP2": "UMAP 2"},
            category_orders={"Cluster": ordered_clusters}
        )
        fig.update_layout(
            autosize=True, margin=dict(l=10, r=10, t=30, b=10),
            legend=dict(title="Cluster", itemsizing="constant", traceorder="normal")
        )
        # for centering the cluster label
        for cluster in ordered_clusters:
            cluster_data = plot_df[plot_df["Cluster"] == cluster]
            if not cluster_data.empty:
                centroid = cluster_data[["UMAP1", "UMAP2"]].mean()
                if not centroid.isna().any():
                    fig.add_annotation(
                        x=centroid["UMAP1"], y=centroid["UMAP2"], text=str(cluster),
                        showarrow=False, font=dict(size=10, weight="bold"), yshift=10
                    )
        logger.info("UMAP plot created successfully")
        return fig
    except Exception as e:
        logger.error("Error creating UMAP plot: %s", str(e))
        return None

###############################
#Function to create Violin Plot
###############################
def create_violin_plot(adata, genes: List[str], cluster_id: str) -> Optional[go.Figure]:
    """Create violin plot for gene expression."""
    adata_source = adata.raw if adata.raw else adata
    var_names_map = {str(name).lower(): name for name in adata_source.var_names}
    available_genes = [var_names_map[str(g).lower()] for g in genes if str(g).lower() in var_names_map]
    sanitized_genes = [f"Gene_{g}" if g[0].isdigit() else g for g in available_genes]
    gene_mapping = dict(zip(sanitized_genes, available_genes))

    if not available_genes:
        logger.error("None of the requested genes found: %s", genes)
        return None

    logger.info("Creating violin plot for genes: %s", available_genes)
    expr = adata_source[:, available_genes].X
    expr = expr.toarray() if not isinstance(expr, np.ndarray) else expr
    expr_data = pd.DataFrame(expr, columns=sanitized_genes, index=adata.obs_names)
    expr_data[cluster_id] = adata.obs[cluster_id].astype(str)

    clusters = sorted(
        expr_data[cluster_id].unique(),
        key=lambda x: int(float(x)) if x.replace(".", "").isdigit() else str(x)
    )
    color_palette = px.colors.qualitative.Plotly
    cluster_colors = {c: color_palette[i % len(color_palette)] for i, c in enumerate(clusters)}

    n_cols = min(len(available_genes), 3)
    n_rows = (len(available_genes) + n_cols - 1) // n_cols
    fig = make_subplots(
        rows=n_rows, cols=n_cols, subplot_titles=[gene_mapping[sg] for sg in sanitized_genes],
        shared_yaxes=True
    )

    # Loop through genes and create the plot
    for idx, sg in enumerate(sanitized_genes):
        row, col = idx // n_cols + 1, idx % n_cols + 1
        for cluster in clusters:
            cluster_data = expr_data[expr_data[cluster_id] == cluster][sg]
            fig.add_trace(
                go.Violin(
                    y=cluster_data, name=cluster, x0=cluster, bandwidth=0.4, width=1.0,
                    showlegend=(idx == 0), box_visible=True,
                    box=dict(line=dict(color="black", width=2), fillcolor="white"),
                    points=False, scalemode="width",
                    line_color=cluster_colors[cluster], fillcolor=cluster_colors[cluster], opacity=0.6
                ),
                row=row, col=col
            )
            fig.add_trace(
                go.Bar(
                    x=[cluster], y=[cluster_data.mean()], width=0.1,
                    marker=dict(color="black", opacity=0.8), showlegend=False
                ),
                row=row, col=col
            )

    fig.update_layout(
        autosize=True, margin=dict(l=10, r=10, t=30, b=10),
        legend_title_text=cluster_id, plot_bgcolor="white", paper_bgcolor="white",
        violingap=0.05, violingroupgap=0.05
    )
    for idx in range(len(available_genes)):
        row, col = idx // n_cols + 1, idx % n_cols + 1
        fig.update_yaxes(range=[0, None], showgrid=True, gridcolor="lightgray", row=row, col=col)
        fig.update_xaxes(tickangle=45, showgrid=False, tickvals=clusters, ticktext=clusters, row=row, col=col)
    logger.info("Violin plot created successfully")
    return fig
    
#############################
#Function to create dotplot
#############################
@lru_cache(maxsize=32)
def create_dot_plot(
    adata_id: int,
    genes: tuple,
    groupby: str,
    expression_cutoff: float = 0.0,
    mean_only_expressed: bool = False,
    log: bool = False,
    layer: Optional[str] = None,
    cmap: str = "Reds"
) -> Tuple[Optional[go.Figure], Optional[str]]:
    """Generate a dot plot with optimized data processing."""
    logger.info("Generating dot plot for genes: %s, groupby: %s", genes, groupby)
    adata = solara_state.value.adata
    if not adata:
        logger.error("No AnnData object loaded")
        return None, "No data loaded."

    var_names_map = {str(name).lower(): name for name in adata.raw.var_names}
    genes_normalized = [str(g).lower() for g in genes]
    available_genes = [var_names_map[gn] for gn in genes_normalized if gn in var_names_map]
    if not available_genes:
        logger.error("No valid genes found: %s", genes)
        return None, f"No valid genes found: {genes}"

    logger.info("Available genes for dot plot: %s", available_genes)
    sanitized_genes = [f"Gene_{g}" if g[0].isdigit() else g for g in available_genes]
    gene_mapping = dict(zip(sanitized_genes, available_genes))

    data = (adata.raw[:, available_genes].X if adata.raw and not layer else
            adata[:, available_genes].layers[layer] if layer else
            adata[:, available_genes].X)
    data = data.toarray() if not isinstance(data, np.ndarray) else data

    clusters = adata.obs[groupby].astype(str).values
    try:
        cluster_order = sorted(np.unique(clusters), key=lambda x: int(float(x)))
    except (ValueError, TypeError):
        cluster_order = sorted(np.unique(clusters))

    plot_data = []
    for gene_idx, (sanitized, orig) in enumerate(zip(sanitized_genes, available_genes)):
        expr = data[:, gene_idx]
        if log:
            expr = np.log1p(expr)

        for cluster in cluster_order:
            mask = clusters == cluster
            cluster_expr = expr[mask]
            expr_cells = cluster_expr > expression_cutoff
            pct_expr = np.sum(expr_cells) / np.sum(mask)
            mean_expr = (
                cluster_expr[expr_cells].mean() if pct_expr > 0 and mean_only_expressed
                else cluster_expr.mean()
            )
            plot_data.append({
                "Gene": orig, "SanitizedGene": sanitized, "Cluster": cluster,
                "PercentExpressing": pct_expr * 100, "AverageExpression": mean_expr
            })

    df = pd.DataFrame(plot_data)
    df["Cluster"] = pd.Categorical(df["Cluster"], categories=cluster_order, ordered=True)
    df["DotSize"] = df["PercentExpressing"].map(
        lambda x: 20 if x >= 90 else 14 if x >= 75 else 9 if x >= 50 else 5 if x >= 25 else 0
    )

    x_categories = list(dict.fromkeys(df["Gene"]))
    x_map = {g: i for i, g in enumerate(x_categories)}
    df["x"] = df["Gene"].map(x_map)
    df["y"] = df["Cluster"].cat.codes

    #Lowest and highest
    vmin = df["AverageExpression"].min()
    vmax = df["AverageExpression"].max()
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["x"], y=df["y"], mode="markers",
        marker=dict(
            size=df["DotSize"], color=df["AverageExpression"], colorscale=cmap,
            colorbar=dict(title="Mean-Expression", x=1.3, tickvals=[vmin, vmax], ticktext=[f"{vmin:.2f}", f"{vmax:.2f}"]),
            showscale=True, line=dict(width=0)
        ),
        text=[f"Gene: {g}<br>Cluster: {c}<br>%: {p:.1f}<br>Mean: {m:.2f}"
              for g, c, p, m in zip(df["Gene"], df["Cluster"], df["PercentExpressing"], df["AverageExpression"])],
        hoverinfo="text"
    ))

    legend_sizes = {"> 75%": 14, "50%": 9, "25%": 5}
    shapes, annotations = [], []
    for i, (label, size) in enumerate(legend_sizes.items()):
        y_pos = 0.9 - i * 0.1
        radius = size / 1000
        shapes.append(dict(
            type="circle", xref="paper", yref="paper",
            x0=1.05 - radius, x1=1.05 + radius, y0=y_pos - radius, y1=y_pos + radius,
            fillcolor="gray", line=dict(width=0), opacity=1
        ))
        annotations.append(dict(
            x=1.10, y=y_pos, xref="paper", yref="paper", text=label,
            showarrow=False, font=dict(size=12), xanchor="left", yanchor="middle"
        ))

    fig.update_layout(
        shapes=shapes, annotations=annotations,
        yaxis=dict(
            title=groupby, tickmode="array",
            tickvals=df["Cluster"].cat.codes, ticktext=df["Cluster"].cat.categories
        ),
        xaxis=dict(
            title="Gene", tickmode="array", tickvals=list(x_map.values()), ticktext=x_categories,
            range=[-0.5, max(x_map.values()) + 0.5]
        ),
        autosize=True, margin=dict(l=10, r=10, t=30, b=10),
        plot_bgcolor="white", paper_bgcolor="white", font=dict(size=14)
    )
    fig.update_yaxes(showgrid=True, gridcolor="lightgray")
    logger.info("Dot plot generated successfully")
    return fig, None

#################################
# Function to create a feature plot
#################################
def create_feature_plot(adata, genes: List[str], cluster_id: str) -> Optional[go.Figure]:
    """Create feature plot for gene expression."""
    adata_source = adata.raw if adata.raw else adata
    var_names_map = {str(name).lower(): name for name in adata_source.var_names}
    available_genes = [var_names_map[str(g).lower()] for g in genes if str(g).lower() in var_names_map]
    sanitized_genes = [f"Gene_{g}" if g[0].isdigit() else g for g in available_genes]
    gene_mapping = dict(zip(sanitized_genes, available_genes))

    if not available_genes:
        logger.error("None of the requested genes found: %s", genes)
        return None

    logger.info("Creating feature plot for genes: %s", available_genes)
    expr = adata_source[:, available_genes].X
    expr = expr.toarray() if not isinstance(expr, np.ndarray) else expr
    expr_data = pd.DataFrame(expr, columns=sanitized_genes, index=adata.obs_names)
    expr_data[cluster_id] = adata.obs[cluster_id].astype(str)
    umap = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    umap_expr = umap.join(expr_data)

    long_df = umap_expr[["UMAP1", "UMAP2", cluster_id] + sanitized_genes].melt(
        id_vars=["UMAP1", "UMAP2", cluster_id],
        value_vars=sanitized_genes, var_name="Gene", value_name="Expression"
    )
    long_df["Gene"] = long_df["Gene"].map(gene_mapping)
    fig = px.scatter(
        long_df, x="UMAP1", y="UMAP2", color="Expression", facet_col="Gene", facet_col_wrap=2,
        color_continuous_scale="Reds", custom_data=[cluster_id],
        hover_data={"UMAP1": ":.2f", "UMAP2": ":.2f", "Expression": ":.2f", cluster_id: True}
    )
    fig.update_traces(
        hovertemplate=(
            "<b>Gene:</b> %{facet_col}<br><b>Cluster:</b> %{customdata[0]}<br>"
            "<b>Expression:</b> %{customdata[1]:.2f}<br><b>UMAP1:</b> %{x:.2f}<br>"
            "<b>UMAP2:</b> %{y:.2f}<extra></extra>"
        ),
        customdata=long_df[[cluster_id, "Expression"]]
    )
    fig.update_layout(autosize=True, margin=dict(l=10, r=10, t=30, b=10))
    fig.update_yaxes(showgrid=True, gridcolor="lightgray")
    fig.update_xaxes(showgrid=False)
    logger.info("Feature plot created successfully")
    return fig
