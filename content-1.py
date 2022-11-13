from matplotlib import pyplot as plt
from matplotlib import patheffects
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter
import seaborn

def get_center(bbox):
    return 0.6 * (bbox.x0 + bbox.x1), 0.6 * (bbox.y0 + bbox.y1)

min_mass = 1
ns_max_mass = 3
bh_min_mass = 3
max_mass = 6
ax = plt.axes(aspect='equal')
ax.set_xlim(min_mass, max_mass)
ax.set_ylim(min_mass, max_mass)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ticks = [min_mass, ns_max_mass, max_mass]
ax.set_xticks(ticks)
ax.set_yticks(ticks)

ticklabels = [r'{} M$_\odot$'.format(tick) for tick in ticks]

ticklabels[2] = r'{} M$_\odot$'.format(50)

#xticklabels[1] =''
ax.set_xticklabels(ticklabels, fontsize=14, fontname="Times New Roman",fontweight="bold")

ticklabels[0] = ''



ax.set_yticklabels(ticklabels, fontsize=14, fontname="Times New Roman",fontweight="bold")

ax.set_xlabel(r'$m_1$', fontsize=20, fontname="Times New Roman",fontweight="bold")
ax.set_ylabel(r'$m_2$', rotation=0, ha='right', fontsize=20, fontname="Times New Roman",fontweight="bold")
ax.xaxis.set_label_coords(1.2, -0.025)
ax.yaxis.set_label_coords(-0.025, 1.1)

bns_color, nsbh_color, bbh_color = seaborn.color_palette(
    'rocket', 3)

p = ax.add_patch(Rectangle((min_mass, min_mass),
                           ns_max_mass - min_mass, ns_max_mass - min_mass,
                           color=bns_color, linewidth=0))
ax.text(0.25 * min_mass + 0.75 * ns_max_mass, 0.5 * min_mass + 0.5 * ns_max_mass,
        'BNS', ha='center', va='center', fontsize=16, fontname="Times New Roman",fontweight="bold", 
        path_effects=[patheffects.Stroke(linewidth=2, foreground=bns_color),
                      patheffects.Normal()])

p = ax.add_patch(Rectangle((bh_min_mass, bh_min_mass),
                           max_mass - bh_min_mass, max_mass - bh_min_mass,
                           color=bbh_color, linewidth=0))
ax.text(0.5 * (bh_min_mass + max_mass), 0.75 * bh_min_mass + 0.25 * max_mass,
        'BBH', ha='center', va='center', fontsize=16, fontname="Times New Roman",fontweight="bold")

p = ax.add_patch(Rectangle((min_mass, bh_min_mass),
                           ns_max_mass - min_mass, max_mass - bh_min_mass,
                           color=nsbh_color, linewidth=0))

p = ax.add_patch(Rectangle((bh_min_mass, min_mass),
                           max_mass - bh_min_mass, ns_max_mass - min_mass,
                           color=nsbh_color, linewidth=0))
ax.text(*get_center(p.get_bbox()), 'NSBH', ha='center', va='center', fontsize=16, fontname="Times New Roman",fontweight="bold")

ax.fill_between([min_mass, max_mass],
                [min_mass, max_mass],
                [max_mass, max_mass],
                color='white', linewidth=0, alpha=0.75, zorder=1.5)
ax.plot([min_mass, max_mass], [min_mass, max_mass], '--k')

#ax.annotate('',
#            xy=(0.975, 1.025), xycoords='axes fraction',
#            xytext=(1.025, 0.975), textcoords='axes fraction',
#            ha='center', va='center',
#            arrowprops=dict(
#                arrowstyle='->', shrinkA=0, shrinkB=0,
#                connectionstyle='angle,angleA=90,angleB=180,rad=7'))

ax.text(0.6, 1.025, '$m_2 \leq m_1$  ',   ha='right', va='center', transform=ax.transAxes, fontsize='large', fontname="Times New Roman")

for args in [[1, 0, 0.20, 0], [0, 1, 0, 0.115]]:
    ax.arrow(*args,
             transform=ax.transAxes, clip_on=False,
             head_width=0.025, head_length=0.025, width=0,
             linewidth=ax.spines['bottom'].get_linewidth(),
             edgecolor=ax.spines['bottom'].get_edgecolor(),
             facecolor=ax.spines['bottom'].get_edgecolor())
    
plt.show()