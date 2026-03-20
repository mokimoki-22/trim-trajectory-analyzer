LANGUAGE_OPTIONS = ["English", "日本語"]

I18N = {
    "English": {
        "upload_tab": "Upload",
        "trajectory_tab": "Trajectory",
        "visualization_tab": "Visualization",
        "export_tab": "Export",
        "upload_header": "Upload scRNA-seq Data",
        "sample_data": "Load Sample Data",
        "run_dpt": "Run Diffusion Pseudotime (DPT)",
        "root_select": "Select Root Cell/Node",
        "trajectory_header": "Trajectory Analysis",
        "vis_header": "Visualization",
        "export_header": "Export Results",
        "umap_plot": "UMAP Plot",
        "heatmap_plot": "Heatmap of Gene Dynamics",
        "bar_chart": "Bar Chart of Pseudotime",
        "apply_dpt": "Run DPT",
        "apply_root": "Set Root",
        "cells_after": "Cells Processed",
        "error_loading": "Failed to load data.",
        "success_msg": "Operation successful!",
        "reset_btn": "Reset All Data",
        "lang_select": "Language / 言語",
    },
    "日本語": {
        "upload_tab": "アップロード",
        "trajectory_tab": "軌跡解析",
        "visualization_tab": "可視化",
        "export_tab": "エクスポート",
        "upload_header": "scRNA-seq データのアップロード",
        "sample_data": "サンプルデータを取得",
        "run_dpt": "Diffusion Pseudotime (DPT) を実行",
        "root_select": "ルートセル/ノードの選択",
        "trajectory_header": "軌跡解析設定",
        "vis_header": "可視化設定",
        "export_header": "結果のエクスポート",
        "umap_plot": "UMAP プロット",
        "heatmap_plot": "遺伝子発現動態のヒートマップ",
        "bar_chart": "擬時間バーチャート",
        "apply_dpt": "DPTを実行",
        "apply_root": "ルートに設定",
        "cells_after": "処理済みセル数",
        "error_loading": "データのロードに失敗しました。",
        "success_msg": "完了しました！",
        "reset_btn": "全てのデータをリセット",
        "lang_select": "Language / 言語",
    }
}

def t(key):
    import streamlit as st
    lang = st.session_state.get("lang", "English")
    return I18N.get(lang, I18N["English"]).get(key, key)
