import streamlit as st
import scanpy as sc
import pandas as pd
import io
import os
import tempfile
import sys
from i18n import t, LANGUAGE_OPTIONS
from trim.utils import log_analysis, load_adata, export_h5ad
from trim.trajectory import compute_dpt
from trim.visualization import (
    plot_pseudotime_umap, plot_gene_dynamics_heatmap, plot_pseudotime_bar
)

# --- PREMIUM UI STYLE ---
st.set_page_config(page_title="TRIM", page_icon="📈", layout="wide")

st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Outfit:wght@400;600&display=swap');
    html, body, [class*="css"] {
        font-family: 'Outfit', sans-serif;
    }
    .main {
        background-color: #F8F9FA;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 24px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        white-space: pre-wrap;
        background-color: #FFFFFF;
        border-radius: 8px 8px 0px 0px;
        padding: 10px 20px;
        font-weight: 600;
        transition: 0.2s;
    }
    .stTabs [aria-selected="true"] {
        border-bottom: 2px solid #007AFF;
        color: #007AFF;
    }
    .stButton>button {
        border-radius: 12px;
        font-weight: 600;
        padding: 0.5rem 2rem;
        background: linear-gradient(135deg, #007AFF 0%, #0056b3 100%);
        color: white;
        border: none;
        transition: all 0.3s;
    }
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 15px rgba(0,122,255,0.3);
    }
</style>
""", unsafe_allow_html=True)

APP_NAME = "TRIM"
APP_SUBTITLE = "Trajectory Analysis in Minutes"

# --- SESSION STATE INITIALIZATION ---
defaults = {
    "adata": None,
    "analysis_log": [],
    "lang": "English",
    "root_cell": None,
    "dpt_done": False,
}
for k, v in defaults.items():
    if k not in st.session_state:
        st.session_state[k] = v

# --- SIDEBAR ---
with st.sidebar:
    st.title(f"📈 {APP_NAME}")
    st.caption(APP_SUBTITLE)
    st.session_state["lang"] = st.selectbox(t("lang_select"), LANGUAGE_OPTIONS, index=0)
    st.divider()
    if st.button(t("reset_btn"), type="secondary"):
        for k in defaults: st.session_state[k] = defaults[k]
        st.rerun()

# --- MAIN TAB NAVIGATION ---
tab_upload, tab_trajectory, tab_vis, tab_export = st.tabs([
    t("upload_tab"), t("trajectory_tab"), t("visualization_tab"), t("export_tab")
])

# --- TAB: UPLOAD ---
with tab_upload:
    st.header(t("upload_header"))
    col1, col2 = st.columns([2, 1])
    
    with col1:
        ftype = st.radio("File Type", ["h5ad", "csv"], horizontal=True)
        uploaded_file = st.file_uploader("Upload single-cell data", type=["h5ad", "csv"])
        
        if uploaded_file and st.button(t("upload_tab"), type="primary"):
            try:
                adata = load_adata(uploaded_file, ftype)
                # DPT requires neighbors and pca, let's pre-check
                if 'X_pca' not in adata.obsm:
                    sc.pp.normalize_total(adata, target_sum=1e4)
                    sc.pp.log1p(adata)
                    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
                    sc.pp.pca(adata, n_comps=50)
                if 'neighbors' not in adata.uns:
                    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
                
                st.session_state["adata"] = adata
                log_analysis(st.session_state["analysis_log"], "Data Loaded", f"Type: {ftype}, Shape: {adata.shape}")
                st.success("Data loaded and pre-processed for trajectory analysis!")
            except Exception as e:
                st.error(f"Error during upload: {e}")

    with col2:
        st.info("Sample Data (PBMC)")
        if st.button(t("sample_data")):
            with st.spinner("Downloading sample data..."):
                try:
                    adata = sc.datasets.pbmc3k_processed()
                    # Ensure UMAP and neighbors are there
                    # pbmc3k_processed already has them
                    st.session_state["adata"] = adata
                    log_analysis(st.session_state["analysis_log"], "Sample Data Loaded", "pbmc3k_processed")
                    st.success("Sample data loaded successfully!")
                except Exception as e:
                    st.error(f"Sample data download failed. Error: {e}")

    if st.session_state["adata"] is not None:
        adata = st.session_state["adata"]
        st.write("---")
        c1, c2, c3 = st.columns(3)
        c1.metric("Cells", adata.n_obs)
        c2.metric("Genes", adata.n_vars)
        if 'leiden' in adata.obs:
            c3.metric("Clusters", len(adata.obs['leiden'].unique()))

# --- TAB: TRAJECTORY ---
with tab_trajectory:
    if st.session_state["adata"] is None:
        st.warning("Please upload data first.")
    else:
        st.header(t("trajectory_header"))
        adata = st.session_state["adata"]
        
        st.subheader(t("root_select"))
        # Let user select root by index or cluster?
        # A simple way for a GUI is selecting one cell index or cluster starting point.
        # Often DPT users choose manually or from a known cluster.
        
        root_method = st.radio("Root Selection Method", ["Cluster starting point", "Specific Cell Index"], horizontal=True)
        
        if root_method == "Cluster starting point" and 'leiden' in adata.obs:
            cluster_id = st.selectbox("Select Cluster to find a root cell", sorted(adata.obs['leiden'].unique()))
            # Pick the cell in that cluster with smallest DPT in a local sense? 
            # Or just the first one?
            root_candidates = adata[adata.obs['leiden'] == cluster_id].obs_names.tolist()
            root_cell = root_candidates[0]
            st.session_state["root_cell"] = root_cell
            st.write(f"Automatically selected root cell index: **{root_cell}**")
        else:
            root_cell = st.text_input("Enter Cell Index (Cell ID)", value=st.session_state["root_cell"] or adata.obs_names[0])
            st.session_state["root_cell"] = root_cell
            
        if st.button(t("apply_dpt"), type="primary"):
            with st.spinner("Calculating Diffusion Pseudotime..."):
                try:
                    adata = compute_dpt(adata, st.session_state["root_cell"])
                    st.session_state["adata"] = adata
                    st.session_state["dpt_done"] = True
                    log_analysis(st.session_state["analysis_log"], "DPT Completed", f"Root: {st.session_state['root_cell']}")
                    st.success(t("success_msg"))
                except Exception as e:
                    st.error(f"DPT failed: {e}")

# --- TAB: VISUALIZATION ---
with tab_vis:
    if not st.session_state["dpt_done"]:
        st.warning("Please run DPT calculation first.")
    else:
        st.header(t("vis_header"))
        adata = st.session_state["adata"]
        
        vc1, vc2 = st.columns([2, 1])
        
        with vc1:
            st.subheader(t("umap_plot"))
            fig_umap = plot_pseudotime_umap(adata)
            st.plotly_chart(fig_umap, use_container_width=True)
            
            st.subheader(t("bar_chart"))
            group_col = st.selectbox("Group by for Pseudotime comparison", 
                                     [c for c in adata.obs.columns if adata.obs[c].dtype.name in ['category', 'object']])
            fig_bar = plot_pseudotime_bar(adata, group_col)
            st.plotly_chart(fig_bar, use_container_width=True)

        with vc2:
            st.subheader(t("heatmap_plot"))
            genes_text = st.text_area("Gene Markers (comma separated)", "CD3E, CD14, GNLY, CD8A, MS4A1")
            genes_to_plot = [g.strip() for g in genes_text.split(",") if g.strip() in adata.var_names]
            
            if genes_to_plot:
                if st.button("Generate Heatmap"):
                    import matplotlib.pyplot as plt
                    fig_heat = plot_gene_dynamics_heatmap(adata, genes_to_plot)
                    st.pyplot(fig_heat)
            else:
                st.caption("Enter valid gene names from your dataset.")

# --- TAB: EXPORT ---
with tab_export:
    st.header(t("export_header"))
    if st.session_state["adata"] is not None:
        adata_final = st.session_state["adata"]
        
        if st.session_state["analysis_log"]:
            st.subheader("Analysis Log")
            st.table(pd.DataFrame(st.session_state["analysis_log"]))
            
        c1, c2 = st.columns(2)
        with c1:
            if st.button("💾 Prepare h5ad Download"):
                with st.spinner("Preparing..."):
                    h5_data = export_h5ad(adata_final)
                    st.download_button("📦 Download Final h5ad", h5_data, "trim_results.h5ad")
        with c2:
            if 'dpt_pseudotime' in adata_final.obs:
                csv_data = adata_final.obs[['dpt_pseudotime']].to_csv().encode('utf-8')
                st.download_button("📊 Download Pseudotime (CSV)", csv_data, "pseudotime_results.csv")
