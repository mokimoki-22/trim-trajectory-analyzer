import scanpy as sc
import anndata
import os
import streamlit as st

def load_adata(uploaded_file, file_type):
    import tempfile
    with tempfile.NamedTemporaryFile(delete=False, suffix=f".{file_type}") as tmp:
        tmp.write(uploaded_file.getbuffer())
        tmp_name = tmp.name
    
    try:
        if file_type == "h5ad":
            adata = sc.read_h5ad(tmp_name)
        elif file_type == "csv":
            adata = sc.read_csv(tmp_name)
        else:
            raise ValueError(f"Unsupported file type: {file_type}")
    finally:
        os.remove(tmp_name)
    
    return adata

def log_analysis(log_list, action, details):
    import datetime
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_list.append({"Time": now, "Action": action, "Details": details})

def export_h5ad(adata):
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        adata.write_h5ad(tmp.name)
        with open(tmp.name, 'rb') as f:
            data = f.read()
        tmp_name = tmp.name
    os.remove(tmp_name)
    return data
