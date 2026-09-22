"""
run_figure3.py - End-to-end driver for the Figure 3 spatial niche analysis.

Runs the numbered scripts in order. Each step reads/writes the intermediate and
output paths defined in 00_config.py, so you can also run any single step on its
own after step 01 has produced the labelled objects.

    python run_figure3.py            # run all steps
    python run_figure3.py 03 06      # run only steps 03 and 06
"""
import os, sys, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))

STEPS = [
    "01_fig3A_senescence_gradient.py",
    "02_fig3BC_cluster_niche_maps.py",
    "03_fig3D_gradient_bin_composition.py",
    "04_fig3E_spatial_marker_maps.py",
    "05_fig3F_gradient_gene_heatmap.py",
    "06_fig3G_niche_deconvolution.py",
    "07_cross_species_niche_markers.py",
    "08_export_hotspot_objects.py",
]


def _load(step_file):
    spec = importlib.util.spec_from_file_location(
        step_file.removesuffix(".py"), os.path.join(HERE, step_file))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(requested=None):
    wanted = {f"{int(s):02d}" for s in requested} if requested else None
    for step_file in STEPS:
        if wanted and step_file[:2] not in wanted:
            continue
        print(f"\n=== {step_file} ===")
        _load(step_file).main()


if __name__ == "__main__":
    main(sys.argv[1:])
