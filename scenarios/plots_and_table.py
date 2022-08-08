import io
from pathlib import Path

# gwpy hijacks Matplotlib's default axes class
# (https://github.com/gwpy/gwpy/issues/1187).
# Take it back and restore default style.
import matplotlib
import gwpy
matplotlib.projections.register_projection(matplotlib.axes.Axes)
matplotlib.rcdefaults()

from astropy.table import join, Table, vstack
from astropy import units as u
from astropy.cosmology import default_cosmology
from astropy.utils.data import download_file
from matplotlib import pyplot as plt
import gwpy.table  # to read ligolw documents as Astropy tables
from IPython.display import Markdown
import numpy as np
from scipy import integrate
from scipy import optimize
from scipy import special
from scipy import stats
from scikits import bootstrap
import seaborn
from tqdm.auto import tqdm

plt.style.use('seaborn-paper')




# Common functions

def betabinom_k_n(k, n):
    return stats.betabinom(n, k + 1, n - k + 1)


@np.vectorize
def poisson_lognormal_rate_cdf(k, mu, sigma):
    lognorm_pdf = stats.lognorm(s=sigma, scale=np.exp(mu)).pdf

    def func(lam):
        prior = lognorm_pdf(lam)
        poisson_pdf = np.exp(special.xlogy(k, lam) - special.gammaln(k + 1) - lam)
        poisson_cdf = special.gammaincc(k + 1, lam)
        return poisson_cdf * prior

    # Marginalize over lambda.
    #
    # Note that we use scipy.integrate.odeint instead
    # of scipy.integrate.quad because it is important for the stability of
    # root_scalar below that we calculate the pdf and the cdf at the same time,
    # using the same exact quadrature rule.
    cdf, _ = integrate.quad(func, 0, np.inf, epsabs=0)
    return cdf



@np.vectorize
def poisson_lognormal_rate_quantiles(p, mu, sigma):
    """Find the quantiles of a Poisson distribution with
    a log-normal prior on its rate.

    Parameters
    ----------
    p : float
        The quantiles at which to find the number of counts.
    mu : float
        The mean of the log of the rate.
    sigma : float
        The standard deviation of the log of the rate.

    Returns
    -------
    k : float
        The number of events.

    Notes
    -----
    This algorithm treats the Poisson count k as a continuous
    real variable so that it can use the scipy.optimize.root_scalar
    root finding/polishing algorithms.
    """
    def func(k):
        return poisson_lognormal_rate_cdf(k, mu, sigma) - p

    if func(0) >= 0:
        return 0

    result = optimize.root_scalar(func, bracket=[0, 1e6])
    return result.root

def format_with_errorbars(mid, lo, hi):
    plus = hi - mid
    minus = mid - lo
    smallest = min(max(0, plus), max(0, minus))

    if smallest == 0:
        return str(mid)
    decimals = 1 - int(np.floor(np.log10(smallest)))

    if all(np.issubdtype(type(_), np.integer) for _ in (mid, lo, hi)):
        decimals = min(decimals, 0)

    plus, minus, mid = np.round([plus, minus, mid], decimals)
    if decimals > 0:
        fstring = '%%.0%df' % decimals
    else:
        fstring = '%d'
    return [fstring % _ for _ in [mid, minus, plus]]


# Settings

alpha = 0.9  # Confidence band for histograms
run_names = run_dirs = ['O4', 'O5']
pops = ['BNS', 'NSBH', 'BBH']  # Populations
classification_names = pops
classification_colors = seaborn.color_palette(n_colors=len(classification_names))
fieldnames = ['area(90)', 'vol(90)', 'distance']
fieldlabels = ['90% cred. area (deg²)',
               '90% cred. comoving volume (10⁶ Mpc³)',
               'Luminosity distance (Mpc)']

cosmo = default_cosmology.get_cosmology_from_string('Planck15')




## Fiducial rates in Gpc$^{-3}$ yr$^{-1}$


# Lower 5% and upper 95% quantiles of log normal distribution
rates_table = Table(
    [
        # BNS rate from GWTC-3, PDB (pair) model
        # https://arxiv.org/pdf/2111.03634.pdf
        {'population': 'BNS', 'lower': 50.0, 'mid': 170.0, 'upper': 440.0},
        # NSBH rate from GWTC-3, PDB (pair) model
        # https://arxiv.org/pdf/2111.03634.pdf
        {'population': 'NSBH', 'lower': 10.0, 'mid': 27.0, 'upper': 58.0},
        # BBH rate from GWTC-3, PDB (pair), model
        # https://arxiv.org/pdf/2111.03634.pdf
        {'population': 'BBH', 'lower': 18.0, 'mid': 25.0, 'upper': 35.0}
    ]
)

standard_90pct_interval, = np.diff(stats.norm.interval(0.9))
rates_table['mu'] = np.log(rates_table['mid'])
rates_table['sigma'] = (np.log(rates_table['upper']) - np.log(rates_table['lower'])) / standard_90pct_interval

print(rates_table)


fiducial_log_rates = np.asarray(rates_table['mu'])
fiducial_log_rate_errs = np.asarray(rates_table['sigma'])


## Load all data sets

tables = {}
with tqdm(total=len(run_dirs) * len(pops)) as progress:
    for run_name, run_dir in zip(run_names, run_dirs):
        for pop in pops:
            path = Path('/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/New-sim-lower-upper-limit-3/runs') / run_dir / (pop.lower() + '_astro')

            allsky = Table.read(str(path / 'allsky.dat'), format='ascii.fast_tab')
            injections = Table.read(str(path / 'injections.dat'), format='ascii.fast_tab')
            allsky.rename_column('coinc_event_id', 'event_id')
            injections.rename_column('simulation_id', 'event_id')
            table = join(allsky, injections)

            # Convert from Mpc^3 to 10^6 Mpc^3
            for colname in ['searched_vol', 'vol(20)', 'vol(50)', 'vol(90)']:
                table[colname] *= 1e-6

            # Get simulated rate from LIGO-LW process table
            process_table = Table.read(str(path / 'events.xml.gz'),
                                       format='ligolw', tablename='process')
            table.meta['rate'] = u.Quantity(process_table[0]['comment'])
            table.meta['network'] = process_table[1]['ifos'].replace('1', '').replace(',', '')

            # Get number of Monte Carlo samples from LIGO-LW process_params table
            process_params_table = Table.read(str(path / 'events.xml.gz'),
                                              format='ligolw', tablename='process_params')
            table.meta['nsamples'], = process_params_table[process_params_table['param'] == '--nsamples']['value'].astype(int)

            tables.setdefault(run_name, {})[pop] = table
            del allsky, injections, table, process_table, process_params_table
            progress.update()
            
## Load old Living Review data sets

old_tables = {}

for pop in tqdm(['BNS', 'NSBH', 'BBH']):
    url_root = f'https://git.ligo.org/emfollow/obs-scenarios-2019-fits-files/-/raw/master/O3_HLV/{pop.lower()}_astro/'
    allsky = Table.read(f'{url_root}/allsky.dat', format='ascii')
    injections = Table.read(f'{url_root}/injections.dat', format='ascii')
    coincs = Table.read(f'{url_root}/coincs.dat', format='ascii')
    table = join(allsky, coincs, 'coinc_event_id')
    table = join(table, injections, 'simulation_id')
    table['vol(90)'] *= 1e-6
    old_tables[pop] = table
    del allsky, injections, coincs, table
    

## Cumulative histograms

run_names_weiz = ['O4/HKLV', 'O5/HKLV']
axs = [plt.subplots()[1] for _ in range(len(fieldnames))]
colors = seaborn.color_palette('colorblind', len(tables))
linestyles = ['-', '--', ':']

for ax, fieldlabel in zip(axs, fieldlabels):
    ax.set_xlabel(fieldlabel)

for ax in axs:
    ax.set_xscale('log')
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.25, 0.50, 0.75, 1])
    ax.set_ylabel('Cumulative fraction of events')
    ax.legend(
        [plt.Line2D([], [], linestyle=linestyle, color='black') for linestyle in linestyles] +
        [plt.Rectangle((0, 0), 0, 0, facecolor=color) for color in colors],
        pops + run_names_weiz)

axs[0].set_xlim(1e-1, 86400)
axs[1].set_xlim(1e-5, 1e4)
axs[2].set_xlim(1e1, 1e4)

ax = axs[2]
zax = ax.twiny()
zax.set_xlim(*ax.get_xlim())
zax.set_xscale(ax.get_xscale())

zax.minorticks_off()
n = np.arange(2, 10)
z = np.concatenate([0.001 * n, 0.01 * n, 0.1 * n, n])
minor = cosmo.luminosity_distance(z).to_value(u.Mpc)
minor = minor[minor > ax.get_xlim()[0]]
minor = minor[minor < ax.get_xlim()[1]]
zax.set_xticks(minor, minor=True)
zax.set_xticklabels([], minor=True)

z = [0.01, 0.1, 1]
zax.set_xticks(cosmo.luminosity_distance(z).to_value(u.Mpc))
zax.set_xticklabels([str(_) for _ in z])
zax.set_xlabel('Redshift')

for irun, (run_name, tables1) in enumerate(tables.items()):
    for ipop, (pop, table) in enumerate(tables1.items()):
        for ifield, fieldname in enumerate(fieldnames):
            data = table[fieldname]
            data = data[np.isfinite(data)]
            ax = axs[ifield]
            t = np.geomspace(*ax.get_xlim(), 100)
            kde = stats.gaussian_kde(np.asarray(np.log(data)))
            (std,), = np.sqrt(kde.covariance)
            y = stats.norm(kde.dataset.ravel(), std).cdf(np.log(t)[:, np.newaxis]).mean(1)
            ax.plot(t, y, color=colors[irun], linestyle=linestyles[ipop])

for ax, fieldname in zip(axs, fieldnames):
    ax.figure.savefig(f'/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/Results/Figures/leo/{fieldname}_new.pdf')
    ax.figure.savefig(f'/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/Results/Figures/leo/{fieldname}_new.svg')
    

"""

## Comparisons with O3 public alerts

o3_data = Table.read('public-alerts.dat', format='ascii')
o3_data['vol(90)'] *= 1e-6

o3_data_by_classification = o3_data.group_by(o3_data['classification']).groups
o3_data_by_classification = dict(zip(o3_data_by_classification.keys, o3_data_by_classification))

print(o3_data)


run_name = 'O3'

fig, axs = plt.subplots(len(['BNS', 'NSBH', 'BBH']), len(fieldnames),
    sharex='col', sharey=True,
    gridspec_kw=dict(bottom=0.08, left=0.08, top=0.92, right=0.95),
    figsize=(7.3, 6))

for ax, fieldlabel in zip(axs[-1], fieldlabels):
    ax.set_xlabel(fieldlabel)
    ax.set_xscale('log')

ax = axs[1][0]
ax.set_ylim(0, 1)
ax.set_yticks([0, 0.25, 0.50, 0.75, 1])
ax.set_ylabel('Cumulative fraction of events')

axs[0, 0].set_xlim(1e0, 86400)
axs[0, 1].set_xlim(1e-3, 1e4)
axs[0, 2].set_xlim(1e1, 1e4)

for pop, color, ax in zip(['BNS', 'NSBH', 'BBH'], classification_colors, axs[:, 0]):
    ax.text(0.05, 0.95, pop, transform=ax.transAxes, color=color, va='top')

for ax in axs[::-1, fieldnames.index('distance')]:
    ax2 = ax.twiny()
    ax2.set_xlim(*ax.get_xlim())
    ax2.set_xscale(ax.get_xscale())

    ax2.minorticks_off()
    n = np.arange(2, 10)
    z = np.concatenate([0.001 * n, 0.01 * n, 0.1 * n, n])
    minor = cosmo.luminosity_distance(z).to_value(u.Mpc)
    minor = minor[minor > ax.get_xlim()[0]]
    minor = minor[minor < ax.get_xlim()[1]]
    ax2.set_xticks(minor, minor=True)
    ax2.set_xticklabels([], minor=True)

    z = [0.01, 0.1, 1]
    ax2.set_xticks(cosmo.luminosity_distance(z).to_value(u.Mpc))
ax2.set_xticklabels([f'$z$={_}' for _ in z])

for ax in axs[::-1, fieldnames.index('area(90)')]:
    ax2 = ax.twiny()
    ax2.set_xlim(*ax.get_xlim())
    ax2.set_xscale(ax.get_xscale())

    ax2.minorticks_off()
    ticks = [3, 9.6, 47]
    ticklabels = ['DECam', 'VRO', 'ZTF']
    ax2.set_xticks(ticks)
ax2.set_xticklabels(ticklabels)
label1, *_ = ax2.xaxis.get_ticklabels()
label1.set_ha('right')

for pop, color, axrow in zip(['BNS', 'NSBH', 'BBH'], classification_colors, axs):
    for fieldname, ax in zip(fieldnames, axrow):
        
        medians = []
        for data, label, linewidth in [
                    [old_tables[pop][fieldname], 'LRR', 0.5 * plt.rcParams['lines.linewidth']],
                    [tables[run_name][pop][fieldname], 'this work', plt.rcParams['lines.linewidth']]
                ]:
            data = data[np.isfinite(data)]
            medians.append(np.median(data))
            kde = stats.gaussian_kde(np.asarray(np.log(data)))
            (std,), = np.sqrt(kde.covariance)
            t = np.geomspace(*ax.get_xlim(), 100)
            y = stats.norm(kde.dataset.ravel(), std).cdf(np.log(t)[:, np.newaxis]).mean(1)
            ax.plot(t, y, color=color, linewidth=linewidth, label=label)

        if fieldname != 'distance':
            ax.annotate(
                '', (medians[1], 0.5), (medians[0], 0.5),
                arrowprops=dict(
                    arrowstyle='-|>', color=color,
                    shrinkA=4, shrinkB=4,
                    linewidth=plt.rcParams['lines.linewidth']))

        data = o3_data_by_classification[pop][fieldname]
        t = np.minimum(np.concatenate(((-np.inf,), np.sort(data))), 10 * ax.get_xlim()[-1])
        y = np.arange(len(data) + 1) / len(data)
        ax.plot(t, y, color='black', drawstyle='steps-post', label='O3 alerts')

axs[-1, -1].legend(frameon=False)
fig.align_labels()
fig.savefig('/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/Results/Figures/leo/o3-comparison.pdf')



## Projections for next several observing runs



linestyles = ['-', '--', ':']

fig, axs = plt.subplots(
    len(pops), len(fieldnames),
    sharex='col', sharey=True,
    gridspec_kw=dict(bottom=0.08, left=0.08, top=0.92, right=0.95),
    figsize=(7.3, 6))

for ax, fieldlabel in zip(axs[-1], fieldlabels):
    ax.set_xlabel(fieldlabel)
    ax.set_xscale('log')

ax = axs[1][0]
ax.set_ylim(0, 1)
ax.set_yticks([0, 0.25, 0.50, 0.75, 1])
ax.set_ylabel('Cumulative fraction of events')

axs[0, 0].set_xlim(1e0, 86400)
axs[0, 1].set_xlim(1e-3, 1e4)
axs[0, 2].set_xlim(1e1, 1e4)

for pop, color, ax in zip(pops, classification_colors, axs[:, 0]):
    ax.text(
        0.05, 0.95, pop,
        transform=ax.transAxes,
        color=color, va='top')

for ax in axs[::-1, fieldnames.index('distance')]:
    ax2 = ax.twiny()
    ax2.set_xlim(*ax.get_xlim())
    ax2.set_xscale(ax.get_xscale())

    ax2.minorticks_off()
    n = np.arange(2, 10)
    z = np.concatenate([0.001 * n, 0.01 * n, 0.1 * n, n])
    minor = cosmo.luminosity_distance(z).to_value(u.Mpc)
    minor = minor[minor > ax.get_xlim()[0]]
    minor = minor[minor < ax.get_xlim()[1]]
    ax2.set_xticks(minor, minor=True)
    ax2.set_xticklabels([], minor=True)

    z = [0.01, 0.1, 1]
    ax2.set_xticks(cosmo.luminosity_distance(z).to_value(u.Mpc))
ax2.set_xticklabels([f'$z$={_}' for _ in z])

for ax in axs[::-1, fieldnames.index('area(90)')]:
    ax2 = ax.twiny()
    ax2.set_xlim(*ax.get_xlim())
    ax2.set_xscale(ax.get_xscale())

    ax2.minorticks_off()
    ticks = [3, 9.6, 47]
    ticklabels = ['DECam', 'VRO', 'ZTF']
    ax2.set_xticks(ticks)
ax2.set_xticklabels(ticklabels)
label1, *_ = ax2.xaxis.get_ticklabels()
label1.set_ha('right')

for pop, color, axrow in zip(pops, classification_colors, axs):
    for fieldname, ax in zip(fieldnames, axrow):
        for run_name, linestyle in zip(run_names, linestyles):
            data = tables[run_name][pop][fieldname]

            data = data[np.isfinite(data)]
            medians.append(np.median(data))
            kde = stats.gaussian_kde(np.asarray(np.log(data)))
            (std,), = np.sqrt(kde.covariance)
            t = np.geomspace(*ax.get_xlim(), 100)
            y = stats.norm(kde.dataset.ravel(), std).cdf(np.log(t)[:, np.newaxis]).mean(1)
            ax.plot(t, y, color=color, linestyle=linestyle, linewidth=linewidth, label=run_name)

axs[-1, -1].legend()
fig.align_labels()
fig.savefig('/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/Results/Figures/leo/predictions.pdf')




linestyles = ['-', '--', ':']

fig, axs = plt.subplots(
    len(pops), len(fieldnames),
    sharex='col', sharey=True,
    gridspec_kw=dict(bottom=0.08, left=0.08, top=0.92, right=0.95),
    figsize=(7.3, 6))

for ax, fieldlabel in zip(axs[-1], fieldlabels):
    ax.set_xlabel(fieldlabel)
    ax.set_xscale('log')

ax = axs[1][0]
ax.set_ylim(1, 1000)
ax.set_yscale('log')
ax.set_ylabel('Cumulative detection rate (events / year)')

axs[0, 0].set_xlim(1e0, 86400)
axs[0, 1].set_xlim(1e-3, 1e4)
axs[0, 2].set_xlim(1e1, 1e4)

for pop, color, ax in zip(pops, classification_colors, axs[:, 0]):
    ax.text(
        0.05, 0.95, pop,
        transform=ax.transAxes,
        color=color, va='top')

for ax in axs[::-1, fieldnames.index('distance')]:
    ax2 = ax.twiny()
    ax2.set_xlim(*ax.get_xlim())
    ax2.set_xscale(ax.get_xscale())

    ax2.minorticks_off()
    n = np.arange(2, 10)
    z = np.concatenate([0.001 * n, 0.01 * n, 0.1 * n, n])
    minor = cosmo.luminosity_distance(z).to_value(u.Mpc)
    minor = minor[minor > ax.get_xlim()[0]]
    minor = minor[minor < ax.get_xlim()[1]]
    ax2.set_xticks(minor, minor=True)
    ax2.set_xticklabels([], minor=True)

    z = [0.01, 0.1, 1]
    ax2.set_xticks(cosmo.luminosity_distance(z).to_value(u.Mpc))
ax2.set_xticklabels([f'$z$={_}' for _ in z])

for ax, pop, fiducial_log_rate in zip(axs[::-1, fieldnames.index('distance')], reversed(pops), reversed(fiducial_log_rates)):
    ax3 = ax.twinx()
    ax3.set_ylim(*ax.get_ylim())
    ax3.set_yscale(ax.get_yscale())
    ax3.set_yticks([len(tables[run_name][pop]) * np.exp(fiducial_log_rate) / tables[run_name][pop].meta['rate'].to_value(u.Gpc**-3 * u.yr**-1) for run_name in run_names])
    ax3.set_yticklabels(run_names)
    ax3.tick_params(length=0)
    ax3.minorticks_off()

for ax in axs[::-1, fieldnames.index('area(90)')]:
    ax2 = ax.twiny()
    ax2.set_xlim(*ax.get_xlim())
    ax2.set_xscale(ax.get_xscale())

    ax2.minorticks_off()
    ticks = [3, 9.6, 47]
    ticklabels = ['DECam', 'VRO', 'ZTF']
    ax2.set_xticks(ticks)
ax2.set_xticklabels(ticklabels)
label1, *_ = ax2.xaxis.get_ticklabels()
label1.set_ha('right')

for pop, color, axrow in zip(pops, classification_colors, axs):
    for fieldname, ax in zip(fieldnames, axrow):
        for run_name, linestyle in zip(reversed(run_names), reversed(linestyles)):
            data = tables[run_name][pop][fieldname]
            rate_row, = rates_table[rates_table['population'] == pop]

            data = data[np.isfinite(data)]
            medians.append(np.median(data))
            kde = stats.gaussian_kde(np.asarray(np.log(data)))
            (std,), = np.sqrt(kde.covariance)
            t = np.geomspace(*ax.get_xlim(), 100)
            y = stats.norm(kde.dataset.ravel(), std).cdf(np.log(t)[:, np.newaxis]).mean(1)
            scale = len(data) / tables[run_name][pop].meta['rate'].to_value(u.Gpc**-3 * u.yr**-1)
            ymid = y * scale * rate_row['mid']
            ylo = y * scale * rate_row['lower']
            yhi = y * scale * rate_row['upper']
            ax.plot(t, ymid, color=color, linestyle=linestyle, linewidth=linewidth, label=run_name)
            ax.fill_between(t, ylo, yhi, color=color, alpha=0.25)

fig.align_labels()
fig.savefig('/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/Results/Figures/leo/annual-predictions.pdf')


## Tabulated statistics

np.random.seed(150914)
statfunc = np.nanmedian
f = io.StringIO()
with open('/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/Results/Figures/leo/summary.rst', 'w') as f_rst:
    print('+-----------+-----------+---------------+---------------+---------------+', file=f_rst)
    print('|           |           | Source class                                  |', file=f_rst)
    print('| Observing |           +---------------+---------------+---------------+', file=f_rst)
    print('| run       | Network   | BNS           | NSBH          | BBH           |', file=f_rst)
    print('+===========+===========+===============+===============+===============+', file=f_rst)
    for fieldname, fieldlabel in zip(fieldnames + ['volume', 'rate'], fieldlabels + ['Sensitive volume (Gpc$^3$)', 'Annual number of detections']):
        print('| {:69s} |'.format(fieldlabel), file=f_rst)
        print('<table>', file=f)
        print('<caption>', fieldlabel, '</caption>', file=f)
        print('<thead>', file=f)
        print('<tr>', file=f)
        print('<th>', 'Run', '</th>', file=f)
        for pop in pops:
            print('<th>', pop, '</th>', file=f)
        print('</tr>', file=f)
        print('</thead>', file=f)
        print('<tbody>', file=f)
        for irun, (run, tables1) in enumerate(tables.items()):
            print('<tr>', file=f)
            print('<th>', run, '</th>', file=f)

            results = {}
            for ipop, (pop, table) in enumerate(tables1.items()):
                rate = table.meta['rate'].to_value(u.Gpc**-3 * u.yr**-1)
                nsamples = table.meta['nsamples']
                fiducial_log_rate = fiducial_log_rates[ipop]
                fiducial_log_rate_err = fiducial_log_rate_errs[ipop]
                mu = fiducial_log_rate + np.log(len(table) / rate)
                sigma = fiducial_log_rate_err

                quantiles = [0.05, 0.5, 0.95]
                if fieldname == 'volume':
                    lo, mid, hi = betabinom_k_n(len(table), nsamples).ppf(quantiles) / rate
                elif fieldname == 'rate':
                    lo, mid, hi = poisson_lognormal_rate_quantiles(quantiles, mu, sigma)
                    lo = int(np.floor(lo))
                    mid = int(np.round(mid))
                    hi = int(np.ceil(hi))
                else:
                    data = table[fieldname]
                    lo, mid, hi = bootstrap.ci(data, statfunc, quantiles)

                mid, lo, hi = format_with_errorbars(mid, lo, hi)
                mathtext = '{}^{{+{}}}_{{-{}}}'.format(mid, hi, lo)
                print('<td>${}$</td>'.format(mathtext), file=f)

                results.setdefault('lo', {})[pop] = lo
                results.setdefault('mid', {})[pop] = mid
                results.setdefault('hi', {})[pop] = hi

            print('</tr>', file=f)
            print('+-----------+-----------+---------------+---------------+---------------+', file=f_rst)
            print('| {:9s} | {:9s} '.format(run, table.meta['network']) + ('| :math:`{:7s}' * 3).format(*results['mid'].values()) + '|', file=f_rst)
            print('|           |           ' + ('| ^{:13s}' * 3).format(*('{{+{}}}'.format(_) for _ in results['hi'].values())) + '|', file=f_rst)
            print('|           |           ' + ('| _{:13s}' * 3).format(*('{{-{}}}`'.format(_) for _ in results['lo'].values())) + '|', file=f_rst)
        print('</tbody>', file=f)
        print('</table>', file=f)
        print('+-----------+-----------+---------------+---------------+---------------+', file=f_rst)
Markdown(f.getvalue())


np.random.seed(150914)
statfunc = np.nanmedian
f = io.StringIO()
with open('/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/Results/Figures/leo/extremes.rst', 'w') as f_rst:
    print('+-----------+-----------+---------------+---------------+---------------+', file=f_rst)
    print('|           |           | Source class                                  |', file=f_rst)
    print('| Observing |           +---------------+---------------+---------------+', file=f_rst)
    print('| run       | Network   | BNS           | NSBH          | BBH           |', file=f_rst)
    print('+===========+===========+===============+===============+===============+', file=f_rst)
    for fieldlabel, statfunc in [
                [
                    'Percentage of events with area(90) <= 5 deg2',
                    lambda _: stats.percentileofscore(_['area(90)'], 5)
                ],
                [
                    'Percentage of events with area(90) <= 20 deg2',
                    lambda _: stats.percentileofscore(_['area(90)'], 20)
                ],
                [
                    'Percentage of events with vol(90) <= 1e3 Mpc3',
                    lambda _: stats.percentileofscore(_['vol(90)'], 1)
                ],
                [
                    'Percentage of events with vol(90) <= 1e4 Mpc3',
                    lambda _: stats.percentileofscore(_['vol(90)'], 10)
                ]
            ]:
        print('| {:69s} |'.format(fieldlabel), file=f_rst)
        print('<table>', file=f)
        print('<caption>', fieldlabel, '</caption>', file=f)
        print('<thead>', file=f)
        print('<tr>', file=f)
        print('<th>', 'Run', '</th>', file=f)
        for pop in pops:
            print('<th>', pop, '</th>', file=f)
        print('</tr>', file=f)
        print('</thead>', file=f)
        print('<tbody>', file=f)
        for irun, (run, tables1) in enumerate(tables.items()):
            print('<tr>', file=f)
            print('<th>', run, '</th>', file=f)
            for ipop, (pop, table) in enumerate(tables1.items()):
                quantiles = [0.05, 0.5, 0.95]
                lo, mid, hi = bootstrap.ci(table, statfunc, quantiles)
                mid, lo, hi = format_with_errorbars(mid, lo, hi)
                mathtext = '{}^{{+{}}}_{{-{}}}'.format(mid, hi, lo)
                print('<td>${}$</td>'.format(mathtext), file=f)

                results.setdefault('lo', {})[pop] = lo
                results.setdefault('mid', {})[pop] = mid
                results.setdefault('hi', {})[pop] = hi

            print('</tr>', file=f)
            print('+-----------+-----------+---------------+---------------+---------------+', file=f_rst)
            print('| {:9s} | {:9s} '.format(run, table.meta['network']) + ('| :math:`{:7s}' * 3).format(*results['mid'].values()) + '|', file=f_rst)
            print('|           |           ' + ('| ^{:13s}' * 3).format(*('{{+{}}}'.format(_) for _ in results['hi'].values())) + '|', file=f_rst)
            print('|           |           ' + ('| _{:13s}' * 3).format(*('{{-{}}}`'.format(_) for _ in results['lo'].values())) + '|', file=f_rst)
        print('</tbody>', file=f)
        print('</table>', file=f)
        print('+-----------+-----------+---------------+---------------+---------------+', file=f_rst)
Markdown(f.getvalue())





"""
