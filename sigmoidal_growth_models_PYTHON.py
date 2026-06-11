import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

# ═══════════════════════════════════════════════════════════════════
#  MODEL FUNCTIONS
# ═══════════════════════════════════════════════════════════════════

def logistic(t, K, r, N0):
    """
    Logistic Growth Model
      N(t)  = K / (1 + ((K-N0)/N0) * exp(-r*t))
      dN/dt = r * N * (1 - N/K)
    """
    return K / (1 + ((K - N0) / N0) * np.exp(-r * t))

def logistic_dNdt(N, K, r):
    return r * N * (1 - N / K)


def gompertz(t, K, r, N0):
    """
    Gompertz Growth Model
      N(t)  = K * exp( -ln(K/N0) * exp(-r*t) )
      dN/dt = r * N * ln(K/N)
    """
    return K * np.exp(-np.log(K / N0) * np.exp(-r * t))

def gompertz_dNdt(N, K, r):
    return r * N * np.log(K / N)


def baranyi_A(t, mu, h0):
    """
    A(t) = t + (1/mu) * ln( exp(-mu*t) + exp(-h0) - exp(-mu*t - h0) )
    """
    return t + (1.0 / mu) * np.log(
        np.exp(-mu * t) + np.exp(-h0) - np.exp(-mu * t - h0)
    )

def baranyi(t, N_max, mu, N0, h0):
    """
    Baranyi-Roberts Model
      N(t)  = N_max * N0 * exp(mu*A(t)) / (N_max - N0 + N0 * exp(mu*A(t)))
      dN/dt = alpha(t) * mu * (1 - N/N_max) * N
      alpha(t) = exp(-mu*t) / (exp(-mu*t) + exp(-h0) - exp(-mu*t-h0))
    """
    At  = baranyi_A(t, mu, h0)
    eAt = np.exp(mu * At)
    return (N_max * N0 * eAt) / (N_max - N0 + N0 * eAt)

def baranyi_dNdt(t, N, N_max, mu, h0):
    num   = np.exp(-mu * t)
    den   = np.exp(-mu * t) + np.exp(-h0) - np.exp(-mu * t - h0)
    alpha = num / den
    return alpha * mu * (1 - N / N_max) * N


def richards(t, K, r, N0, theta):
    """
    Richard's (Generalised Logistic) Model
      a     = (K/N0)^theta - 1
      N(t)  = K * [1 + a * exp(-r*theta*t)]^(-1/theta)
      dN/dt = r * N * [1 - (N/K)^theta]
      Note: when theta=1 this reduces exactly to the Logistic model.
    """
    a = (K / N0) ** theta - 1
    return K * (1 + a * np.exp(-r * theta * t)) ** (-1.0 / theta)

def richards_dNdt(N, K, r, theta):
    return r * N * (1 - (N / K) ** theta)


def von_bertalanffy(t, N_inf, K, N0):
    """
    Von Bertalanffy Model
      N(t)  = N_inf - (N_inf - N0) * exp(-K*t)
      dN/dt = K * (N_inf - N)
    """
    return N_inf - (N_inf - N0) * np.exp(-K * t)

def von_bertalanffy_dNdt(N, N_inf, K):
    return K * (N_inf - N)


# ═══════════════════════════════════════════════════════════════════
#  PARAMETERS  (matching updated Excel workbook — all N0 = 5)
# ═══════════════════════════════════════════════════════════════════

L_K,  L_r,  L_N0                = 1000.0, 1.0, 5.0   # Logistic
G_K,  G_r,  G_N0                = 1000.0, 1.0, 5.0   # Gompertz
BR_Nmax, BR_mu, BR_N0, BR_h0    = 1000.0, 1.0, 5.0, 5.0  # Baranyi-Roberts
RI_K,  RI_r,  RI_N0, RI_theta   = 1000.0, 1.0, 5.0, 1.0  # Richard's
VB_Ninf, VB_K, VB_N0            = 1000.0, 1.0, 5.0   # Von Bertalanffy

t_cont     = np.linspace(0, 30, 2000)
t_discrete = np.arange(0, 31)


# ═══════════════════════════════════════════════════════════════════
#  COMPUTE CURVES
# ═══════════════════════════════════════════════════════════════════

N_log  = logistic(t_cont,        L_K,  L_r,  L_N0)
N_gomp = gompertz(t_cont,        G_K,  G_r,  G_N0)
N_br   = baranyi(t_cont,         BR_Nmax, BR_mu, BR_N0, BR_h0)
N_ri   = richards(t_cont,        RI_K,  RI_r,  RI_N0, RI_theta)
N_vb   = von_bertalanffy(t_cont, VB_Ninf, VB_K, VB_N0)

dN_log  = logistic_dNdt(N_log,  L_K, L_r)
dN_gomp = gompertz_dNdt(N_gomp, G_K, G_r)
dN_br   = baranyi_dNdt(t_cont,  N_br, BR_Nmax, BR_mu, BR_h0)
dN_ri   = richards_dNdt(N_ri,   RI_K, RI_r, RI_theta)
dN_vb   = von_bertalanffy_dNdt(N_vb, VB_Ninf, VB_K)

Nd_log  = logistic(t_discrete,        L_K,  L_r,  L_N0)
Nd_gomp = gompertz(t_discrete,        G_K,  G_r,  G_N0)
Nd_br   = baranyi(t_discrete,         BR_Nmax, BR_mu, BR_N0, BR_h0)
Nd_ri   = richards(t_discrete,        RI_K,  RI_r,  RI_N0, RI_theta)
Nd_vb   = von_bertalanffy(t_discrete, VB_Ninf, VB_K, VB_N0)


# ═══════════════════════════════════════════════════════════════════
#  PRINT SIMULATION TABLES
# ═══════════════════════════════════════════════════════════════════

def print_table(name, params, t_arr, Nt_fn, Nt1_fn, dNdt_fn):
    print("=" * 72)
    print(f"{name}  |  {params}")
    print("=" * 72)
    print(f"{'t':>4}  {'N(t)':>14}  {'N(t+1)':>14}  {'dN/dt':>14}")
    print("-" * 54)
    for ti in t_arr[:-1]:
        Ni   = Nt_fn(ti)
        Ni1  = Nt1_fn(ti)
        dNi  = dNdt_fn(ti, Ni)
        print(f"{ti:>4.0f}  {Ni:>14.6f}  {Ni1:>14.6f}  {dNi:>14.6f}")
    print()

print_table(
    "LOGISTIC", f"K={L_K:.0f}  r={L_r}  N(0)={L_N0}", t_discrete,
    lambda t: logistic(t, L_K, L_r, L_N0),
    lambda t: logistic(t+1, L_K, L_r, L_N0),
    lambda t, N: logistic_dNdt(N, L_K, L_r)
)
print_table(
    "GOMPERTZ", f"K={G_K:.0f}  r={G_r}  N(0)={G_N0}", t_discrete,
    lambda t: gompertz(t, G_K, G_r, G_N0),
    lambda t: gompertz(t+1, G_K, G_r, G_N0),
    lambda t, N: gompertz_dNdt(N, G_K, G_r)
)
print_table(
    "BARANYI-ROBERTS",
    f"N_max={BR_Nmax:.0f}  mu_max={BR_mu}  N(0)={BR_N0}  h0={BR_h0}", t_discrete,
    lambda t: baranyi(t, BR_Nmax, BR_mu, BR_N0, BR_h0),
    lambda t: baranyi(t+1, BR_Nmax, BR_mu, BR_N0, BR_h0),
    lambda t, N: baranyi_dNdt(t, N, BR_Nmax, BR_mu, BR_h0)
)
print_table(
    "RICHARD'S", f"K={RI_K:.0f}  r={RI_r}  N(0)={RI_N0}  theta={RI_theta}", t_discrete,
    lambda t: richards(t, RI_K, RI_r, RI_N0, RI_theta),
    lambda t: richards(t+1, RI_K, RI_r, RI_N0, RI_theta),
    lambda t, N: richards_dNdt(N, RI_K, RI_r, RI_theta)
)
print_table(
    "VON BERTALANFFY", f"N_inf={VB_Ninf:.0f}  K={VB_K}  N(0)={VB_N0}", t_discrete,
    lambda t: von_bertalanffy(t, VB_Ninf, VB_K, VB_N0),
    lambda t: von_bertalanffy(t+1, VB_Ninf, VB_K, VB_N0),
    lambda t, N: von_bertalanffy_dNdt(N, VB_Ninf, VB_K)
)


# ═══════════════════════════════════════════════════════════════════
#  PLOT
# ═══════════════════════════════════════════════════════════════════

BG     = "#0f1117"
PANEL  = "#1a1d27"
GRID_C = "#2a2d3a"
TEXT_C = "#e0e0e0"

COLORS = {
    "Logistic":         "#4c9be8",
    "Gompertz":         "#4ecb71",
    "Baranyi-Roberts":  "#f5a623",
    "Richard's":        "#e84c8b",
    "Von Bertalanffy":  "#b57bee",
}

def style_ax(ax, title, xlabel="Time (t)", ylabel=""):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=TEXT_C, labelsize=8.5)
    ax.set_title(title, color=TEXT_C, fontsize=10, fontweight="bold", pad=7)
    ax.set_xlabel(xlabel, color=TEXT_C, fontsize=8.5)
    ax.set_ylabel(ylabel, color=TEXT_C, fontsize=8.5)
    ax.spines[:].set_color(GRID_C)
    ax.grid(color=GRID_C, linewidth=0.55, linestyle="--", alpha=0.7)

fig = plt.figure(figsize=(16, 18), facecolor=BG)
gs  = gridspec.GridSpec(4, 2, figure=fig, hspace=0.52, wspace=0.32)

# ── Panel 1: All 5 models overlaid ──────────────────────────────
ax0 = fig.add_subplot(gs[0, :])
overlay = [
    ("Logistic",        N_log,  Nd_log),
    ("Gompertz",        N_gomp, Nd_gomp),
    ("Baranyi-Roberts", N_br,   Nd_br),
    ("Richard's",       N_ri,   Nd_ri),
    ("Von Bertalanffy", N_vb,   Nd_vb),
]
for label, Nc, Ndc in overlay:
    c = COLORS[label]
    ax0.plot(t_cont, Nc, color=c, lw=2, label=label)
    ax0.scatter(t_discrete, Ndc, color=c, s=18, zorder=5, alpha=0.8)
ax0.axhline(1000, color="white", lw=0.7, ls=":", alpha=0.35, label="Asymptote = 1000")
style_ax(ax0, "All Five Sigmoidal Growth Models — N(t)  [N(0)=5 for all models]",
         ylabel="N(t)")
ax0.legend(facecolor=PANEL, labelcolor=TEXT_C, fontsize=9,
           framealpha=0.85, ncol=3, loc="lower right")

# ── Panels 2–6: Individual N(t) + dN/dt on twin axis ────────────
individual = [
    ("Logistic",        N_log,  dN_log,  Nd_log,  gs[1, 0]),
    ("Gompertz",        N_gomp, dN_gomp, Nd_gomp, gs[1, 1]),
    ("Baranyi-Roberts", N_br,   dN_br,   Nd_br,   gs[2, 0]),
    ("Richard's",       N_ri,   dN_ri,   Nd_ri,   gs[2, 1]),
    ("Von Bertalanffy", N_vb,   dN_vb,   Nd_vb,   gs[3, 0]),
]
for label, Nc, dNc, Ndc, gspec in individual:
    c  = COLORS[label]
    ax = fig.add_subplot(gspec)
    ax.plot(t_cont, Nc, color=c, lw=2, label="N(t)")
    ax.fill_between(t_cont, Nc, alpha=0.10, color=c)
    ax.scatter(t_discrete, Ndc, color=c, s=18, zorder=5, alpha=0.85)
    ax.axhline(1000, color="white", lw=0.6, ls=":", alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(t_cont, dNc, color=c, lw=1.3, ls="--", alpha=0.7, label="dN/dt")
    ax2.set_ylabel("dN/dt", color=c, fontsize=7.5)
    ax2.tick_params(colors=c, labelsize=7.5)
    ax2.spines[:].set_color(GRID_C)
    style_ax(ax, f"{label} — N(t)  &  dN/dt", ylabel="N(t)")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, facecolor=PANEL, labelcolor=TEXT_C,
              fontsize=7.5, framealpha=0.8)

# ── Panel 7: dN/dt comparison ────────────────────────────────────
ax_dN = fig.add_subplot(gs[3, 1])
for label, Nc, dNc, _, __ in individual:
    ax_dN.plot(t_cont, dNc, color=COLORS[label], lw=1.8, label=label)
style_ax(ax_dN, "Growth Rate Comparison — dN/dt", ylabel="dN/dt")
ax_dN.legend(facecolor=PANEL, labelcolor=TEXT_C, fontsize=7.5, framealpha=0.85)

fig.suptitle(
    "Sigmoidal Growth Models  ·  Logistic · Gompertz · Baranyi-Roberts · Richard's · Von Bertalanffy\n"
    "Parameters: K=1000, r=1.0, N(0)=5 for all models",
    color="white", fontsize=12, fontweight="bold", y=0.998)

out_png = Path(__file__).resolve().parent / "sigmoidal_growth_models_plot.png"
plt.savefig(out_png, dpi=155, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close()
print(f"Plot saved: {out_png}")
