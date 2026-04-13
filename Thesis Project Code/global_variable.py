import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd

# PLOT FUNCTIONS:
def create_folder(func_folder_path):
    if not os.path.exists(func_folder_path):
        os.makedirs(func_folder_path, exist_ok=True)

def single_histogram(x, xlabel, savepath, normalize=False, ylim=(), grid=True, legendsize=20,
                        title='', label='', xlogscale=False, bins=50, ylogscale=False,
                     x_edges=(), ylabel='Count', markmedian=False, figsize=(6, 6),
                     x_ticks_numbers=0, y_ticks_numbers=0, save_format='jpg', xline=(), line_label='', fontsize=34):
    plt.figure(figsize=figsize)
    if not x_edges:
        x_edges = (min(x), max(x))
    else:
        x = x[(x >= x_edges[0]) & (x <= x_edges[1])]
    if xlogscale:
        bins = np.logspace(float(np.log10(float(x_edges[0]))), float(np.log10(float(x_edges[1]))), num=bins)
        plt.xscale('log')
    if not label:
        label = f'N = {sum(~np.isnan(x))}'
    if normalize:
        counts, bins, patches = plt.hist(x, label=label, bins=bins, density=True)
        plt.ylabel(ylabel=normalized_histogram_label, fontsize=fontsize)
    else:
        counts, bins, patches = plt.hist(x, label=label, bins=bins)
        plt.ylabel(ylabel=ylabel, fontsize=fontsize)
    if xline:
        plt.plot(xline, [0, max(counts)], label=line_label, linestyle='--', color='black')
    plt.xlabel(xlabel, fontsize=fontsize)
    if grid:
        plt.grid()
    plt.xlim(x_edges)
    if ylim:
        plt.ylim(ylim)
    plt.tick_params(axis='both', which='both', labelsize=fontsize)
    if title:
        plt.title(title, fontsize=fontsize)
    if ylogscale:
        plt.yscale('log')
    if markmedian:
        ymin, ymax = plt.gca().get_ylim()
        median =  np.median(x[~np.isnan(x)])
        rounded_median = str(round(median, 2))
        plt.plot([median, median], [max(counts), ymax], color='black', linestyle='--',
                 linewidth=3, label=f'Median = {rounded_median}')
    if type(x_ticks_numbers) is not int:
        x_ticks_string = [str(tick) for tick in x_ticks_numbers]
        plt.xticks(x_ticks_numbers, x_ticks_string, fontsize=fontsize)
    if type(y_ticks_numbers) is not int:
        y_ticks_string = [str(tick) for tick in y_ticks_numbers]
        plt.yticks(y_ticks_numbers, y_ticks_string, fontsize=fontsize)
    plt.legend(fontsize=legendsize)
    plt.savefig(savepath, format=save_format, bbox_inches='tight')
    plt.close()


def single_scatter_plot(x, y, xlabel, ylabel, savepath, z=None, zlabel='', cmap='gist_rainbow',
                        title='', xlogscale=False, ylogscale=False, pointsize=15,
                        flip_y_axis=False, alpha=1, x_line=None, y_line=None, line_label='', x_line2=None, y_line2=None,
                        xlim=[], ylim=[], x_ticks_numbers=[], y_ticks_numbers=[], flip_x_axis=False,
                        linewidth=2, grid=True, save_format='png', figsize=(8, 8), line_label2='',
                        current_fontsize=16, x_point=(), y_point=(), pointlabel='', current_ticksize=14):
    if xlim:
        if ylim:
            label = f'N = {sum(~np.isnan(x[(x >= xlim[0]) & (x <= xlim[1]) & (y >= ylim[0]) & (y <= ylim[1])]))}'
        else:
            label = f'N = {sum(~np.isnan(x[(x >= xlim[0]) & (x <= xlim[1])]))}'
    else:
        if ylim:
            label = f'N = {sum(~np.isnan(y[(y >= ylim[0]) & (y <= ylim[1])]))}'
        else:
            label = f'N = {sum(~np.isnan(x))}'
    plt.figure(figsize=figsize)
    if z is None:
        plt.scatter(x, y, label=label, s=pointsize, alpha=alpha)
    else:
        scatter = plt.scatter(x, y, c=z, cmap=cmap, s=pointsize, alpha=alpha,
                              label=label)
        colorbar = plt.colorbar(scatter)
        colorbar.set_label(zlabel, fontsize=current_fontsize)
        colorbar.ax.tick_params(labelsize=current_fontsize)
    if flip_x_axis:
        plt.gca().invert_xaxis()
    plt.xlabel(xlabel, fontsize=current_fontsize)
    plt.ylabel(ylabel, fontsize=current_fontsize)
    plt.tick_params(axis='both', which='both', labelsize=current_fontsize)
    if x_point:
        plt.scatter(x_point, y_point, s=40, alpha=alpha, color='red', label=pointlabel)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    if xlogscale:
        plt.xscale('log')
    if ylogscale:
        plt.yscale('log')
    if title:
        plt.title(title, fontsize=26)
    if x_line2 is not None:
        plt.plot(x_line2, y_line2, linestyle='--', color='red', label=line_label2, linewidth=linewidth)
    if x_line is not None:
        plt.plot(x_line, y_line, linestyle='--', color='black', label=line_label, linewidth=linewidth)
    if x_ticks_numbers:
        x_ticks_string = [str(tick) for tick in x_ticks_numbers]
        plt.xticks(x_ticks_numbers, x_ticks_string, fontsize=current_ticksize)
    if y_ticks_numbers:
        y_ticks_string = [str(tick) for tick in y_ticks_numbers]
        plt.yticks(y_ticks_numbers, y_ticks_string, fontsize=current_ticksize)
    if flip_y_axis:
        plt.gca().invert_yaxis()
    if grid:
        plt.grid()
    plt.legend(fontsize=current_fontsize)
    plt.savefig(savepath, format=save_format, bbox_inches='tight')
    plt.close()

# LABELS:
depth_ratio_label = 'Depth Ratio'
distance_label = 'Distance [pc]'
vbroad_significance_label = 'Vbroad / $\u03C3_{vbroad}$'
normalized_histogram_label = 'Frequency'
abs_mag_label = '$M_{G,0}$ [mag]'
bp_rp_label = r'$(G_\mathrm{BP}-G_\mathrm{RP})_0$ [mag]'
regular_histogram_label = 'sources count'
vsini_label = 'v*sin(i) [km/s]'
vbroad_label = 'Vbroad [km/s]'
vbroad_error_label = 'Vbroad Error [km/s]'
temperature_label = r'Temperature$^{\circ}\mathrm{[k]}$'
cdf_label = 'Cumulative\nDistribution Function'
period_label = 'Period [day]'
period_error_label = 'Period Error [day]'
mass_label = 'Stellar Mass [m$_{\u2299}$]'
radius_label = 'Stellar Radius [R$_{\u2299}$]'
expected_velocity_label = 'Fitted vbroad [km/s]'
global_ranking_label = 'Global Ranking'
logpost_MSC_label = 'logposterior_MSC'
vbroad_over_velocity_label = 'Vbroad / Fitted vbroad'
ruwe_label = 'RUWE'
g_excess_label = 'G Excess [mag]'
apparent_g_label = 'Apparent G Magnitude [mag]'
rvgof_label = 'rv_renormalised_gof'
rvchi2_label = 'rv_chisq_pvalue'
rvamp_label = 'rv_amplitude_robust [km/s]'
e_cosw_label = r'e * cos($\omega$)'
period_ratio_label = 'Rotational Period / Orbital Period'
metalicity_label = '[Fe/H]'
vbroad_ratio_label = 'Gaia / GALAH Vbroad'
rotational_velocity_label = 'Estimated Rotational Velocity [km/s]'
vbroad_over_rotational_velocity_label = 'Vbroad / Rotational Velocity'
wavelength_angstrom_label = 'Wavelength $[\AA]$'
normalised_flux_label = 'Normalised Flux'