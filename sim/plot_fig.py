# -*- coding: utf-8 -*-
"""plot_fig.py -- fig_curvature_sweeps for the Spectral Curvature paper.

Panels:
  (a) detuning schedules Delta(t)/Delta0 vs t/T (representative runs)
  (b) curvature profiles Lambda(t)/Lambda(t*): analytic (Lemma 1) lines
      + finite-difference ||[B, Bdot]||_op markers (self-check #1 visually)
  (c) leakage P_LZ vs normalized peak curvature Lambda(t*)/Omega^3 for all
      90 sweep runs, with BOTH closed-form overlays: the standard-LZ form
      exp(-pi Omega^3 / (4 Lambda*)) and the paper's printed form
      exp(-pi Omega^3 / (2 Lambda*)).

Colors: validated reference palette slots 1-3 (dataviz skill), fixed
order; theory lines wear ink/muted tones, never series colors.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SIM_DIR = os.path.dirname(os.path.abspath(__file__))
BUILD_DIR = os.path.dirname(SIM_DIR)
RESULTS_DIR = os.path.join(BUILD_DIR, "results")
FIG_DIR = os.path.join(BUILD_DIR, "figures")

# validated palette (light mode), fixed slot order
C_SCHED = {"linear": "#2a78d6", "tanh": "#1baf7a", "gauss_flat": "#eda100"}
LABELS = {"linear": "linear", "tanh": "tanh", "gauss_flat": "Gaussian-flattened"}
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#ffffff"          # print target: paper-white surface
MARKERS_OMEGA = {0.5: "o", 1.0: "s", 2.0: "^"}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 7.5,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.7,
    "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 6.8, "ytick.labelsize": 6.8,
    "grid.color": GRID, "grid.linewidth": 0.5,
    "legend.frameon": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


def style_axis(ax):
    ax.grid(True, which="major", zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    with open(os.path.join(RESULTS_DIR, "lzsm_sweeps.json"),
              encoding="utf-8") as f:
        res = json.load(f)
    with open(os.path.join(RESULTS_DIR, "rep_profiles.json"),
              encoding="utf-8") as f:
        prof = json.load(f)

    fig, axes = plt.subplots(1, 3, figsize=(7.08, 2.45),
                             constrained_layout=True)
    ax_a, ax_b, ax_c = axes

    # ------------------------------------------------ (a) schedules
    for kind in ("linear", "tanh", "gauss_flat"):
        p = prof[kind]
        ax_a.plot(p["t_over_T"], p["Delta_over_Delta0"],
                  color=C_SCHED[kind], lw=1.6, zorder=3)
    ax_a.set_xlabel(r"$t/T$")
    ax_a.set_ylabel(r"$\Delta(t)/\Delta_0$")
    ax_a.set_xlim(-0.5, 0.5)
    ax_a.set_ylim(-1.12, 1.12)
    # direct labels (relief rule for sub-3:1 slots)
    ax_a.text(0.335, 0.30, "linear", color=C_SCHED["linear"],
              fontsize=6.8, ha="left")
    ax_a.text(-0.13, 0.72, "tanh", color="#0d7a54",
              fontsize=6.8, ha="right")
    ax_a.text(0.02, -0.62, "Gaussian-\nflattened", color="#a06d00",
              fontsize=6.8, ha="left")
    style_axis(ax_a)
    ax_a.set_title("(a) detuning schedules", fontsize=7.8, color=INK)

    # ---------------------------------------- (b) curvature profiles
    for kind in ("linear", "tanh", "gauss_flat"):
        p = prof[kind]
        ax_b.plot(p["t_over_T"], p["lam_analytic_over_peak"],
                  color=C_SCHED[kind], lw=1.6, zorder=3)
        ax_b.plot(p["t_over_T_subsampled"],
                  p["lam_fd_over_peak_subsampled"],
                  ls="none", marker="o", ms=3.4, mfc="none",
                  mec=C_SCHED[kind], mew=0.8, zorder=4)
    ax_b.set_xlabel(r"$t/T$")
    ax_b.set_ylabel(r"$\Lambda(t)\,/\,\Lambda(t^*)$")
    ax_b.set_xlim(-0.5, 0.5)
    ax_b.set_ylim(0, 4.4)
    ax_b.text(0.28, 1.22, "linear", color=C_SCHED["linear"], fontsize=6.8)
    ax_b.text(-0.085, 0.42, "tanh", color="#0d7a54", fontsize=6.8,
              ha="right")
    ax_b.text(0.16, 3.55, "Gaussian-\nflattened", color="#a06d00",
              fontsize=6.8)
    ax_b.text(-0.47, 4.05,
              "lines: Lemma 1   markers: FD $\\|[B,\\dot B]\\|$",
              color=MUTED, fontsize=6.0)
    style_axis(ax_b)
    ax_b.set_title("(b) curvature profiles", fontsize=7.8, color=INK)

    # ------------------------------- (c) leakage vs peak curvature
    u = np.geomspace(0.085, 4.4, 300)
    ax_c.plot(u, np.exp(-np.pi / (4.0 * u)), color=INK, lw=1.4, zorder=2,
              label=r"standard LZ  $e^{-\pi\Omega^3/4\Lambda^*}$")
    ax_c.plot(u, np.exp(-np.pi / (2.0 * u)), color=MUTED, lw=1.2,
              ls=(0, (4, 2.5)), zorder=2,
              label=r"paper Eq. (5)  $e^{-\pi\Omega^3/2\Lambda^*}$")
    # All 90 runs are plotted. The three Omega branches (and the three
    # schedules) land on the SAME (Lambda*/Omega^3, P) points -- the
    # Omega-branch spread is < 6e-13 relative (exact scale invariance) --
    # so markers are nested by size to keep every schedule visible
    # (a bullseye = all three schedules coincide).
    SIZES = {"gauss_flat": 6.4, "tanh": 4.5, "linear": 2.6}
    from matplotlib.lines import Line2D
    for kind in ("gauss_flat", "tanh", "linear"):   # big first, small on top
        rows = [r for r in res["sweep_rows"] if r["schedule"] == kind]
        xs = [r["Lambda_star_over_Omega3"] for r in rows]
        ys = [r["P_sim"] for r in rows]
        ax_c.plot(xs, ys, ls="none", marker="o", ms=SIZES[kind],
                  mfc=C_SCHED[kind], mec="white", mew=0.5, zorder=4)
    sh = [Line2D([], [], ls="none", marker="o", ms=SIZES[k],
                 mfc=C_SCHED[k], mec="white", mew=0.5, label=LABELS[k])
          for k in ("linear", "tanh", "gauss_flat")]
    leg1 = ax_c.legend(handles=sh, loc="lower right", fontsize=6.0,
                       handletextpad=0.15, borderaxespad=0.2,
                       labelspacing=0.25)
    ax_c.add_artist(leg1)
    ax_c.legend(loc="upper left", fontsize=5.8, handlelength=1.5,
                handletextpad=0.35, borderaxespad=0.1, labelspacing=0.3)
    # (the 90-runs / Omega-collapse note lives in the figure caption)
    ax_c.set_xscale("log")
    ax_c.set_yscale("log")
    ax_c.set_xlim(0.085, 4.4)
    ax_c.set_ylim(1.2e-4, 1.6)
    ax_c.set_xlabel(r"$\Lambda(t^*)\,/\,\Omega^3$")
    ax_c.set_ylabel(r"$P_{\mathrm{LZ}}$  (TDSE)")
    style_axis(ax_c)
    ax_c.set_title("(c) leakage vs peak curvature", fontsize=7.8,
                   color=INK)

    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"fig_curvature_sweeps.{ext}"),
                    dpi=300)
    print("wrote", os.path.join(FIG_DIR, "fig_curvature_sweeps.{png,pdf}"))


if __name__ == "__main__":
    main()
