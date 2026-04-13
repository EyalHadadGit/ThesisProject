from global_variable import *
from Gaussian_Function_Fitting import create_fitting

# variables range: Temp in [5600, 8000]; rad in [1.25, 3]; q < 0.7 and ecosw < 0.1

# general
depth_Ratio_threshold = 0.7
current_savepth = r'C:\Users\ASUS\Desktop\Thesis Project Code\Figures'
folder_name = fr'depth ratio below {depth_Ratio_threshold} temperature in 5600 8000 rad in 1.25 3'
nummy_datatable = pd.read_csv(f'EBs with vbroad on MS.csv')
current_fontsize = 24
vbroad_ticks = np.linspace(30, 300, 100)

# variables from datatable (used for filtering conditions)
depth_ratio = nummy_datatable['derived_secondary_ecl_depth'] / nummy_datatable['derived_primary_ecl_depth']
teff, rad = nummy_datatable['teff_gspphot'], nummy_datatable['radius_gspphot']
period, vbroad = nummy_datatable['frequency'] ** -1, nummy_datatable['vbroad']
ecosw = abs(np.pi * 0.5 * (abs(nummy_datatable['derived_secondary_ecl_phase'] -
                nummy_datatable['derived_primary_ecl_phase']) - 0.5))

# the two filtering conditions:
cond = ((teff < 8000) & (teff > 5600) & (rad < 3) & (rad > 1.25) & (depth_ratio < depth_Ratio_threshold) & (ecosw < 0.1) &
        (period < 3) & (period > 0.3) & (vbroad < 300) & (vbroad > 30))
cond_with_ecosw_outliers = ((teff < 8000) & (teff > 5600) & (rad < 3) & (rad > 1.25)
                            & (depth_ratio < depth_Ratio_threshold) &
                            (period < 3) & (period > 0.3) & (vbroad < 300) & (vbroad > 30))

# variables with outlier cond:
vbroad_with_outlier, period_with_outliers, ecosw_with_outliers, source_id_with_outliers =\
    (nummy_datatable[cond_with_ecosw_outliers]['vbroad'],
     nummy_datatable[cond_with_ecosw_outliers]['frequency'] ** -1,
     ecosw[cond_with_ecosw_outliers],
     nummy_datatable[cond_with_ecosw_outliers]['source_id_1'])

# variables of outliers:
vbroad_of_outliers, period_of_outliers = (nummy_datatable[cond_with_ecosw_outliers & ~cond]['vbroad'],
                                          nummy_datatable[cond_with_ecosw_outliers & ~cond]['frequency'] ** -1)

# reducing the datatable to cond without outliers:
nummy_datatable_with_cond = nummy_datatable[cond]

# variables without outliers:
vbroad, period, g_ranking, teff, rad = (nummy_datatable_with_cond['vbroad'], 1 / nummy_datatable_with_cond['frequency'],
                                        nummy_datatable_with_cond['global_ranking'], nummy_datatable_with_cond['teff_gspphot'],
                                        nummy_datatable_with_cond['radius_gspphot'])
ecosw = np.pi * 0.5 * abs(abs(nummy_datatable_with_cond['derived_secondary_ecl_phase'] -
                              nummy_datatable_with_cond['derived_primary_ecl_phase']) - 0.5)
vbroad_significance = vbroad / nummy_datatable_with_cond['e_Vbroad']
radius = nummy_datatable_with_cond['radius_gspphot']
radius_MAD_uncertainty = 1.46 * np.median(np.abs(np.median(radius) - radius))


# call 'create fitting'. generates the algorithm with its solution:
slope, intercept, scatter_with_solution_fig = create_fitting(vbroad, period, g_ranking, teff, '', f'{current_savepth}/{folder_name}',
                                                             ecosw, vbroad_with_outlier, period_with_outliers, ecosw_with_outliers)

# ecosw histogram:
single_histogram(ecosw[cond_with_ecosw_outliers], r'$|ecos\omega|$',
                 xlogscale=True, savepath=f'{current_savepth}/ecosw histogram with grid.jpg',
                 bins=25, figsize=(12, 12), xline=(0.1, 0.1), line_label=r'$|ecos\omega|$ = 0.7')

# prints the data with solution, and with a marked outlier:
plt.subplots(figsize=(12, 12))
plt.scatter(period[cond], vbroad[cond], alpha=0.3, label=f'N = {sum(cond)}')
plt.scatter(period_of_outliers, vbroad_of_outliers, label=f'N = {len(vbroad_of_outliers)}',
            marker='s', color='black', s=60)
plt.scatter(1, 2 * np.pi * np.median(radius) * (696000 / (24 * 60 * 60)), color='green', s=140)
plt.errorbar(1, 2 * np.pi * np.median(radius) * (696000 / (24 * 60 * 60)), color='green',
             label=r'$(\frac{2\pi \mathrm{R_{med}}}{\mathrm{1 [day]}}); \mathrm{R_{med}} = 2.018$',
             yerr=2 * np. pi * radius_MAD_uncertainty * (696000 / (24 * 60 * 60)), fmt='o', capsize=14, elinewidth=4)

plt.plot(10 ** ((np.log10(vbroad_ticks) - intercept) / slope), vbroad_ticks, color='red', linestyle='--', label=f'logV=AlogP+B')
plt.ylabel(vbroad_label, fontsize=current_fontsize)
plt.xlabel(period_label, fontsize=current_fontsize)
plt.legend(fontsize=current_fontsize)
plt.xscale('log')
plt.yscale('log')
plt.xticks([0.3, 0.4, 0.6, 1, 2, 3], ['0.3', '0.4', '0.6', '1', '2', '3'], rotation=45, fontsize=current_fontsize)
plt.yticks([30, 40, 60, 100, 200, 300], ['30', '40', '60', '100', '200', '300'], fontsize=current_fontsize)
plt.grid()
plt.savefig(fr'{current_savepth}\vbroad vs period with solution.jpg', bbox_inches='tight')

# Now I want the same figure, with a marked brown dot and a slightly different legend:

plt.gca().get_legend().remove()
plt.scatter(2.745, 128.678, marker='s', color='brown', s=60, alpha=1, label=f'N = 1', zorder=1)
plt.legend(["N = 977", "N = 72", r'$(\frac{2\pi \mathrm{R_{med}}}{\mathrm{1 [day]}}); \mathrm{R_{med}} = 2.018$',
            'logV=AlogP+B', 'N = 1'], fontsize=current_fontsize)

plt.savefig(fr'{current_savepth}\vbroad vs period with outlier.jpg', bbox_inches='tight')
