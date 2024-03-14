import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.axes import Axes
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from pymatgen.core.composition import Composition
import mpltern
from scipy.constants import pi
import numpy as np

def get_gs(df, Tl = 'Tliq', jez=False):
    df['K_W_Tc'] = (df.Tc - df.Tg) / df[Tl] # best one in the paper
    df['K_W_Tx'] = (df.Tx - df.Tg) / df[Tl]
    df['gamma_Tc'] = df.Tc / (df.Tg+df[Tl])
    df['H_prime_Tx'] = (df.Tx - df.Tg) / df.Tg
    df['K_H_Tc'] = (df.Tc - df.Tg) / (df[Tl] - df.Tc) # replaced Tmelt with Tliq
    df['H_prime_Tc'] = (df.Tc - df.Tg) / df.Tg
    df['K_H_Tx'] = (df.Tx - df.Tg) / (df[Tl] - df.Tx) # replaced Tmelt with Tliq
    df['deltaT_rg'] = (df.Tx - df.Tg) / (df[Tl] - df.Tg)
    df['K_cr'] = (df[Tl] - df.Tx) / (df[Tl] - df.Tg)
    if jez:
        df['Jezica'] = (df.ViscosityAtTl) - 2 * np.log10(df[Tl])
    return df
    
def get_eta_tl(df, Tl = 'Tliq'):
    return df['log10 (η∞)'] + (12-df['log10 (η∞)'])*(df.T12/df[Tl])*np.exp((df.m/(12-df['log10 (η∞)'])-1)*(df.T12/df[Tl] - 1))

def get_gfa(df, logXs = -2, logNs = 3, g=pi,  Tl = 'Tliq', **kw):
    Umax = 10 ** df.log_Umax
    
    tn = (10**logXs / (g * 10**logNs * Umax**2))**(1 / 2)

    df['GFA'] = -np.log10((df[Tl] - df.T_Umax) / tn)
    return df

# general plotting styles-related settings
plt.style.use("seaborn-v0_8-pastel")
rcParams.update(
    {
        "font.family": "Helvetica",
        "mathtext.fontset": "stixsans",
        "axes.labelsize": 20,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "xtick.major.size": 7,
        "ytick.major.size": 7,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "font.size": 20,
        "axes.linewidth": 1.0,
        "lines.dashed_pattern": (5, 2.5),
    }
)


def latexify(comp: Composition):
    c = Composition(comp)
    latex_str = c.to_latex_string().replace("$_{1}$", "")
    return latex_str


def plot_ternary(
    df: pd.DataFrame,
    column: str,
    ax: Axes,
    corner_columns: list = None,
    corner_labels: list = None,
    cmap: str = "coolwarm",
    vmax: float = None,
    vmin: float = None,
    title: str = None,
):
    if corner_labels is None:
        corner_labels = [latexify(col.split()[0]) for col in corner_columns]

    if vmax is None:
        #vmax = max(df[df["total_filter_score"] >= 0]["salt_cation_loading"])
        vmin = min(df[column])
    if vmin is None:
        vmax = max(df[column])
    norm = Normalize(vmin=vmin, vmax=vmax)

    # set ternary labels
    ax.set_tlabel(corner_labels[0])
    ax.set_llabel(corner_labels[1])
    ax.set_rlabel(corner_labels[2])
    if title is not None:
        ax.set_title(title, pad=20, bbox={"boxstyle": "round", "fc": "lightgrey"})

    # plot filter score < 0 points
    #df1 = df[df["total_filter_score"] < 0]
    #top = df1[corner_columns[0]]
    #left = df1[corner_columns[1]]
    #right = df1[corner_columns[2]]
    #ax.scatter(
    #    top,
    #    left,
    #    right,
    #    c="g",
    #    edgecolors="dimgrey",
    #    linewidths=0.25,
    #    s=12,
    #    zorder=100,
    #    clip_on=False,
    #)

    # plot filter score = 0 points, color per mass loading
    #df1 = df[df["total_filter_score"] >= 0]
    df1 = df.copy()
    top = df1[corner_columns[0]]
    left = df1[corner_columns[1]]
    right = df1[corner_columns[2]]
    pc = ax.scatter(
        top,
        left,
        right,
        s=15,
        #c=df1["salt_cation_loading"],
        c=df1[column],
        cmap=cmap,
        norm=norm,
        zorder=100,
        clip_on=False,
    )
   
    # add sciglass boundary
    if corner_columns == ["P2O5 at %", "Na2O at %", "Fe2O3 at %"]:
        top = [60, 60, 65, 65]
        left = [40, 0, 0, 35]
        right = [0, 40, 35, 0]
        ax.fill(top, left, right, '', alpha=0.25)
        top = [60, 40, 40]
        left = [0, 0, 20]
        right = [40, 60, 40]
        ax.fill(top, left, right, 'gray', alpha=0.25) 
        top = [60, 40, 60]
        left = [0, 20, 20]
        right = [40, 40, 20]
        ax.fill(top, left, right, 'green', alpha=0.25)

        top = [65, 65]
        left = [0, 35]
        right = [35, 0]
        ax.plot(top, left, right, '-k', linewidth=1.5)
        
        top = [60,60]
        left = [40,20]
        right = [0,20]
        ax.plot(top, left, right, '-k', linewidth=1.5)

        top = [60, 40]
        left = [20, 20]
        right = [20, 40]
        ax.plot(top, left, right, '-k', linewidth=1.5)

        top = [40, 40]
        left = [20, 0]
        right = [40,60]
        ax.plot(top, left, right, '-k', linewidth=1.5)

    if corner_columns == ["SiO2", "Na2O", "B2O3"]:
        top = [50, 0]
        left = [50, 42]
        right = [0, 68]
        ax.plot(top, left, right, '-k', linewidth=1.5 )

        top = [50, 0, 0, 100]
        left = [50, 42, 0, 0]
        right = [0, 68, 100, 0]
        ax.fill( top, left, right, 'gray', alpha=0.25)
    return pc

def make_ternaries(
    csv_path: str = "WS_Na2O_ADD_Fe2O3_P2O5.csv",
    column: str = "salt_cation_loading",
    column_label: str = 'Salt Cation Loading',
    corner_columns: list = None,
    corner_labels: list = None,
    anchor_column: str = None,
    anchor_values: list[float] = None,
    cmap: str = "coolwarm",
    save_to_disk: bool = True,
    output_dir: str = "figures",
    filename: str = None,
):
    df = pd.read_csv(csv_path)
    df = df[ (df[ corner_columns[0] ] != 0.0) & (df[ corner_columns[1] ] != 0.0) ] 
    df = df[ (df[ corner_columns[1] ] != 0.0) & (df[ corner_columns[2] ] != 0.0) ] 
    df = df[ (df[ corner_columns[0] ] != 0.0) & (df[ corner_columns[2] ] != 0.0) ] 

    if anchor_column is None:
        #vmax = max(df[df["total_filter_score"] >= 0]["salt_cation_loading"])
        vmin = min(df[column])
        vmax = max(df[column])
        fig = plt.figure(figsize=(6 * 1.2, 6))
        gs = fig.add_gridspec(ncols=2, width_ratios=[1, 0.2], wspace=0.25)
        ax = fig.add_subplot(gs[0], projection="ternary")
        plot_ternary(
            df,
            column,
            ax,
            corner_columns=corner_columns,
            corner_labels=corner_labels,
            cmap=cmap,
            vmax=vmax,
            vmin=vmin
        )
    else:
        if anchor_values is None:
            anchor_values = [0, 10, 20]
        ncols = len(anchor_values)
        width_ratios = [1 - a / 100 for a in anchor_values]
        fig = plt.figure(figsize=(6 * (sum(width_ratios) + 0.2), 6))
        gs = fig.add_gridspec(
            ncols=ncols + 1, width_ratios=width_ratios + [0.2], wspace=0.25
        )
        vmax = 0
        for val in anchor_values:
            sub_df = df[
                (df[anchor_column] < val + 1e-6) & (df[anchor_column] > val - 1e-6)
            ]
            _vmax = max(
                #sub_df[sub_df["total_filter_score"] >= 0]["salt_cation_loading"]
                sub_df[column]
            )
            _vmin = min(
                #sub_df[sub_df["total_filter_score"] >= 0]["salt_cation_loading"]
                sub_df[column]
            )

            if _vmax > vmax:
                vmax = _vmax
            if _vmin > vmin:
                vmin = _vmin

        for idx in range(ncols):
            print(f"subplot #{idx}, anchor_value: {anchor_values[idx]}")
            ax = fig.add_subplot(
                gs[idx], projection="ternary", ternary_scale=100 - anchor_values[idx]
            )
            # filter according to anchor
            sub_df = df[
                (df[anchor_column] < anchor_values[idx] + 1e-6)
                & (df[anchor_column] > anchor_values[idx] - 1e-6)
            ]
            print(sub_df)
            title = f"{latexify(anchor_column.split()[0])}: {anchor_values[idx]} at%"
            plot_ternary(
                sub_df,
                ax,
                corner_columns=corner_columns,
                corner_labels=corner_labels,
                cmap=cmap,
                vmax=vmax,
                vmin=vmin,
                title=title,
            )

    # color bar
    # cax = ax.inset_axes([1.025, 0.25, 0.025, 0.7], transform=ax.transAxes)
    cax = fig.add_subplot(gs[-1])
    cax.set_position([0.875, 0.2, 0.05, 0.65])
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap)
    mycolorbar = plt.colorbar(sm, cax=cax, norm=norm, aspect=100, fraction=0.046, pad=0.04)
    #mycolorbar = plt.colorbar(cax, cmap=cmap, norm=norm, aspect=50)
    print( column_label )
    mycolorbar.set_label(f"{column_label}", rotation=270, va="baseline")

    if save_to_disk:
        if filename is None:
            if anchor_column is None:
                filename = f"{Path(csv_path).stem}.png"
            else:
                tag = anchor_column.split()[0]
                filename = f"{Path(csv_path).stem}_ANCHOR_{tag}.png"
        if not Path(output_dir).is_dir():
            os.mkdir(output_dir)
        fig_path = Path(output_dir) / filename
        print(fig_path)
        fig.savefig(fig_path, bbox_inches="tight", pad_inches=0.2, dpi=300)
