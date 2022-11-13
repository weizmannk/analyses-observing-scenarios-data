from matplotlib import pyplot as plt
import numpy as np

ALPHA_1 = -2.16
ALPHA_2 = -1.46
A = 0.97
M_GAP_LO = 2.72
M_GAP_HI = 6.13
ETA_GAP_LO = 50
ETA_GAP_HI = 50
ETA_MIN = 50
ETA_MAX = 4.91
BETA = 1.89
M_MIN = 1.16
M_MAX = 54.38


def lopass(m, m_0, eta):
   return 1 / (1 + (m / m_0)**eta)


def hipass(m, m_0, eta):
   return 1 - lopass(m, m_0, eta)


def bandpass(m, m_lo, m_hi, eta_lo, eta_hi, A):
   return 1 - A * hipass(m, m_lo, eta_lo) * lopass(m, m_hi, eta_hi)


def pairing_function(m1, m2):
   m1, m2 = np.maximum(m1, m2), np.minimum(m1, m2)
   return np.where((m1 <= 60) | (m2 >= 2.5), (m2 / m1) ** BETA, 0)


def mass_distribution_1d(m):
   return (
      bandpass(m, M_GAP_LO, M_GAP_HI, ETA_GAP_LO, ETA_GAP_HI, A) *
      hipass(m, M_MIN, ETA_MIN) *
      lopass(m, M_MAX, ETA_MAX) *
      (m / M_GAP_HI) ** np.where(m < M_GAP_HI, ALPHA_1, ALPHA_2)
   )


def mass_distribution_2d(m1, m2):
   return (
      mass_distribution_1d(m1) *
      mass_distribution_1d(m2) *
      pairing_function(m1, m2)
   )


# Plot 1D distribution of component mass.

m = np.geomspace(1, 100, 200)
fig, ax = plt.subplots()
ax.set_xscale('log')
ax.plot(m, m * mass_distribution_1d(m))
ax.set_xlim(1, 100)
ax.set_ylim(0, None)
ax.set_xlabel(r'Component mass, $m$ ($M_\odot$)')
ax.set_ylabel(r'$m \, p(m|\lambda)$')
ax.set_yticks([])
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xscale(ax.get_xscale())
ax2.set_xticks([M_MIN, M_GAP_LO, M_GAP_HI, M_MAX])
ax2.set_xticklabels([r'$M_\mathrm{min}$',
                     r'$M^\mathrm{gap}_\mathrm{low}$',
                     r'$M^\mathrm{gap}_\mathrm{high}$',
                     r'$M_\mathrm{max}$'])
ax2.set_xticks([], minor=True)
ax2.grid(axis='x')
fig.show()

# Plot joint 2D distribution of m1, m2.

m1, m2 = np.meshgrid(m, m)
fig, ax = plt.subplots(subplot_kw=dict(aspect=1))
ax.set_xlim(1, 100)
ax.set_ylim(1, 100)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'Primary mass, $m_1$ ($M_\odot$)')
ax.set_ylabel(r'Secondary mass, $m_2$ ($M_\odot$)')
img = ax.pcolormesh(m, m, m1 * m2 * mass_distribution_2d(m1, m2),
                    vmin=0, vmax=25, shading='gouraud', rasterized=True)
cbar = plt.colorbar(img, ax=ax)
cbar.set_label(r'$m_1 \, m_2 \, p(m_1, m_2 | \Lambda)$')
cbar.set_ticks([])

ax.fill_between([1, 100],
                [1, 100],
                [100, 100],
                color='white', linewidth=0, alpha=0.75, zorder=1.5)
ax.plot([1, 100], [1, 100], '--k')

ax.annotate('',
            xy=(0.975, 1.025), xycoords='axes fraction',
            xytext=(1.025, 0.975), textcoords='axes fraction',
            ha='center', va='center',
            arrowprops=dict(
                arrowstyle='->', shrinkA=0, shrinkB=0,
                connectionstyle='angle,angleA=90,angleB=180,rad=7'))
ax.text(0.975, 1.025, '$m_1 \geq m_2$ by definition  ',
        ha='right', va='center', transform=ax.transAxes, fontsize='small')

fig.show()