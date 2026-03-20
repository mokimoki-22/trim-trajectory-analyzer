import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import io

def plot_pseudotime_umap(adata):
    """
    Returns a Plotly figure showing UMAP colored by dpt_pseudotime.
    """
    if 'dpt_pseudotime' not in adata.obs:
        return None
    
    df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
    df['Pseudotime'] = adata.obs['dpt_pseudotime'].values
    df['Index'] = adata.obs_names
    
    fig = px.scatter(df, x='UMAP1', y='UMAP2', color='Pseudotime', 
                     hover_data=['Index'], color_continuous_scale='Viridis')
    fig.update_layout(title="Pseudotime UMAP", 
                      plot_bgcolor='white', 
                      xaxis=dict(showgrid=False), 
                      yaxis=dict(showgrid=False))
    return fig

def plot_gene_dynamics_heatmap(adata, genes):
    """
    Returns a Matplotlib figure for a heatmap of gene expression vs pseudotime.
    """
    # Sort by pseudotime
    adata_sorted = adata[adata.obs['dpt_pseudotime'].argsort(), :]
    
    import seaborn as sns
    expression = adata_sorted[:, genes].X
    if hasattr(expression, "toarray"):
        expression = expression.toarray()
    
    df = pd.DataFrame(expression, columns=genes)
    df.index = np.arange(len(df)) # Pseudotime order
    
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(df.T, cmap='viridis', ax=ax, xticklabels=False)
    ax.set_xlabel("Pseudotime (cells sorted)")
    ax.set_ylabel("Genes")
    plt.tight_layout()
    return fig

def plot_pseudotime_bar(adata, group_by='leiden'):
    """
    Boxplot/Bar chart of pseudotime per group.
    """
    df = pd.DataFrame({
        'Pseudotime': adata.obs['dpt_pseudotime'],
        'Group': adata.obs[group_by]
    })
    
    fig = px.box(df, x='Group', y='Pseudotime', color='Group', title=f"Pseudotime by {group_by}")
    fig.update_layout(showlegend=False, 
                      plot_bgcolor='white', 
                      xaxis=dict(showgrid=False), 
                      yaxis=dict(showgrid=False))
    return fig
