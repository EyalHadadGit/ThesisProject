from emcee.autocorr import AutocorrError
from global_variable import *
import os
import h5py
import corner
import emcee

# FUNCTIONS:

# CURRENTLY UNUSED: this function computes the number of points within ellipse.
# def count_points_in_elipse(x_data, y_data, x0, slope, intercept, y_length, x_length):
#     major = max(y_length, x_length)
#     minor = min(y_length, x_length)
#     c = np.sqrt(major ** 2 - minor ** 2) / 2
#     logx0 = np.log10(x0)
#     y0 = slope * logx0 + intercept
#     if y_length <= x_length:
#         alpha_angle = np.arctan(slope)
#     else:
#         alpha_angle = np.arctan(-slope ** -1)
#     left_foci_x = logx0 - c * np.cos(alpha_angle)
#     right_foci_x = logx0 + c * np.cos(alpha_angle)
#     left_foci_y = y0 - c * np.sin(alpha_angle)
#     right_foci_y = y0 + c * np.sin(alpha_angle)
#     logx_data = np.log10(x_data)
#     logy_data = np.log10(y_data)
#     sum_dist_from_foci = np.sqrt((logx_data - right_foci_x) ** 2 + (logy_data - right_foci_y) ** 2) +\
#                          np.sqrt((logx_data - left_foci_x) ** 2 + (logy_data - left_foci_y) ** 2)
#     print(f'{sum(sum_dist_from_foci <= major)} / {len(sum_dist_from_foci)}')

# this function calculates the probability score of a certain configuration
def prob_calc(x, y, x0, y0, theta, sigmax, sigmay, p):
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    sin_2_theta = np.sin(2 * theta)
    a = cos_theta ** 2 / (2 * sigmax ** 2) + sin_theta ** 2 / (2 * sigmay ** 2)
    b = -sin_2_theta * (1 / (4 * sigmax ** 2) - 1 / (4 * sigmay ** 2))
    c = sin_theta ** 2 / (2 * sigmax ** 2) + cos_theta ** 2 / (2 * sigmay ** 2)
    gauss_normalization = 2 * np.pi * sigmay * sigmax
    q = (1 - gauss_normalization * p) / graph_area
    prob_result = (p * np.exp(- (a * (x - x0) ** 2 + 2 * b * (x - x0) * (y - y0) + c * (y - y0) ** 2))
                   + q)
    return np.log(prob_result)

# this function takes a given set of parameters {parameters_set} and the data (period, vbroad) and deals with edge cases (sigma < 0; non-finite probability scores)
def log_prob(parameters_set, log_prob_vbroad, log_prob_period):
    x0, y0, theta, sigmax, sigmay, p = parameters_set
    if sigmax < 0 or sigmay < 0 or p <= 0 or 2 * p * np.pi * sigmay * sigmax >= 1:
        return -np.inf
    try:
        probabilities_vector = prob_calc(log_prob_period,
                                         log_prob_vbroad, x0, y0, theta, sigmax, sigmay, p)
    except RuntimeWarning:
        return -np.inf
    if sum(~np.isfinite(probabilities_vector)): # this checks if any of the probabilities vector is non-finite
        return -np.inf
    log_prob_result = sum(probabilities_vector)
    if ~np.isfinite(log_prob_result): # again checks that log_prob doesn't return -inf (the non of the probabilities is negative)
        return -np.inf
    return log_prob_result

# this function makes a scatter plot with the solution of the MCMC algorithm:
def scatter_plot_with_solution(x_func, y_func, x0, y0, slopes,
                               l, w, savepath='', title='', cosw='', markoutliers=False):
    create_folder(savepath)
    fig = plt.figure(figsize=(12, 12))
    if markoutliers:
        plt.scatter(x_func[cosw < 0.1], y_func[cosw < 0.1], alpha=0.3,
                              color='blue', s=30, label=f'N = {len(x_func[cosw < 0.1])}')
        plt.scatter(x_func[cosw > 0.1], y_func[cosw > 0.1],
                    marker='s', color='black', s=60, alpha=1, label=f'N = {len(x_func[cosw > 0.1])}', zorder=1)
        plt.scatter(x_func[cosw > 0.1], y_func[cosw > 0.1] - 0.05,
                    marker='_', color='black', s=60, alpha=1, zorder=2)
        plt.scatter(x_func[cosw > 0.1], y_func[cosw > 0.1] + 0.05,
                    marker='_', color='black', s=60, alpha=1, zorder=2)
    else:
        scatter = plt.scatter(x_func, y_func, c=cosw, cmap='viridis_r', s=30, alpha=1, label=f'N = {len(x_func)}')
        colorbar = plt.colorbar(scatter)
        colorbar.set_label(r'$ecos\omega$', fontsize=30)
        colorbar.ax.tick_params(labelsize=30)

    plt.xlabel(f'{period_label}', fontsize=24)
    plt.ylabel(vbroad_label, fontsize=24)
    if title:
        plt.title(title, fontsize=24)
    plt.xscale('log')
    plt.yscale('log')
    x_ticks_numbers = [0.3, 0.4, 0.6, 1, 1.5, 2, 2.5, 3]
    x_ticks_string = [str(tick) for tick in x_ticks_numbers]
    y_ticks_numbers = [40, 60, 100, 150, 200, 250, 300]
    y_ticks_string = [str(tick) for tick in y_ticks_numbers]
    plt.xticks(x_ticks_numbers, x_ticks_string, fontsize=20, rotation=45)
    plt.yticks(y_ticks_numbers, y_ticks_string, fontsize=20)
    for ind in range(len(x0)): # draws lines and ellipses
        current_x0 = np.log10(x0[ind])
        current_y0 = np.log10(y0[ind])
        current_slope = slopes[ind]
        current_intercept = current_y0 - current_slope * current_x0
        x0_vec = np.linspace(period_lower_cutoff, period_upper_cutoff, 100)
        y0_vec = 10 ** (current_slope * np.log10(x0_vec) + current_intercept)
        plt.plot(x0_vec, y0_vec, linewidth=2, linestyle='--', color='red', alpha=0.1)
        current_sigmay = w[ind]
        current_sigmax = l[ind]
        num_sigma = 3
        x_edges, y_edges, x_ellipse_points, y_ellipse_points = get_ellipse_edges(current_slope, current_intercept,
                                                                                 num_sigma * current_sigmay,
                                                                                 num_sigma * current_sigmax, 10 ** current_x0,
                                                                                 x_vec_num=100)
        plt.plot(x_ellipse_points, y_ellipse_points, '', fillstyle='none', color='red', alpha=0.05)
    plt.xlim([period_lower_cutoff, period_upper_cutoff])
    plt.ylim([vbroad_lower_cutoff, vbroad_upper_cutoff])
    plt.legend(loc='upper right', fontsize=18)
    plt.savefig(f'{savepath}/scatter with solution.pdf')
    plt.close()
    return fig

def get_ellipse_edges(slope, intercept, y_length, x_length, x0, x_vec_num=1000):
    logx0 = np.log10(x0)
    y_length_half = y_length / 2
    x_length_half = x_length / 2
    y0 = slope * logx0 + intercept
    alpha_angle = np.arctan(slope)
    c0 = (np.cos(alpha_angle) ** 2) / x_length_half ** 2 + (np.sin(alpha_angle) ** 2) / y_length_half ** 2
    c1 = np.sin(alpha_angle) ** 2 / x_length_half ** 2 + np.cos(alpha_angle) ** 2 / y_length_half ** 2
    c2 = np.sin(2 * alpha_angle) / x_length_half ** 2 - np.sin(2 * alpha_angle) / y_length_half ** 2
    c3 = -2 * logx0 * np.cos(alpha_angle) ** 2 / x_length_half ** 2 - y0 * np.sin(
        2 * alpha_angle) / x_length_half ** 2 - 2 * logx0 * np.sin(alpha_angle) ** 2 / y_length_half ** 2 + y0 * np.sin(
        2 * alpha_angle) / y_length_half ** 2
    c4 = - logx0 * np.sin(2 * alpha_angle) / x_length_half ** 2 - 2 * y0 * np.sin(
        alpha_angle) ** 2 / x_length_half ** 2 + logx0 * np.sin(2 * alpha_angle) / y_length_half ** 2 - 2 * y0 * np.cos(
        alpha_angle) ** 2 / y_length_half ** 2
    c5 = logx0 ** 2 * np.cos(alpha_angle) ** 2 / x_length_half ** 2 + logx0 * y0 * np.sin(
        2 * alpha_angle) / x_length_half ** 2 + y0 ** 2 * np.sin(alpha_angle) ** 2 / x_length_half ** 2 \
         + logx0 ** 2 * np.sin(alpha_angle) ** 2 / y_length_half ** 2 - logx0 * y0 * np.sin(
        2 * alpha_angle) / y_length_half ** 2 + y0 ** 2 * np.cos(alpha_angle) ** 2 / y_length_half ** 2 - 1

    # this treats ellipse equation as a pollinom of y and finds the x values that give discriminant=0;
    # a_x, b_x, c_x are the coefficients of the pollinom in x that is inside the discriminant:
    a_x, b_x, c_x = c2 ** 2 - 4 * c1 * c0, 2 * (c2 * c4 - 2 * c1 * c3), c4 ** 2 - 4 * c1 * c5
    x_discriminant = np.sqrt(b_x ** 2 - 4 * a_x * c_x)
    x_low, x_high = min((x_discriminant - b_x) / (2 * a_x), (-x_discriminant - b_x) / (2 * a_x)), max(
        (x_discriminant - b_x) / (2 * a_x), (-x_discriminant - b_x) / (2 * a_x))

    # this treats ellipse equation as a pollinom of x and finds the y values that give discriminant=0
    # a_y, b_y, c_y are the coefficients of the pollinom in y that is inside the discriminant:
    a_y, b_y, c_y = (c2 ** 2 - 4 * c0 * c1), (2 * c2 * c3 - 4 * c0 * c4), (c3 ** 2 - 4 * c0 * c5)
    y_discriminant = np.sqrt(b_y ** 2 - 4 * a_y * c_y)
    y_low, y_high = min((y_discriminant - b_y) / (2 * a_y), (-y_discriminant - b_y) / (2 * a_y)), \
        max((y_discriminant - b_y) / (2 * a_y), (-y_discriminant - b_y) / (2 * a_y))

    # this gives the frame:
    x = np.linspace(x_low, x_high, x_vec_num)
    a_y, b_y, c_y = c1, c2 * x + c4, c0 * x ** 2 + c3 * x + c5
    y_discriminant = np.sqrt(b_y ** 2 - 4 * a_y * c_y)
    y_1, y_2 = (-b_y - y_discriminant) / (2 * a_y), (-b_y + y_discriminant) / (2 * a_y)
    x = np.hstack((x, np.flip(x)))
    y = np.hstack((y_1, np.flip(y_2)))
    return (10 ** x_low, 10 ** x_high), (10 ** y_low, 10 ** y_high), np.array([10 ** x_point for x_point in x]),\
        np.array([10 ** y_point for y_point in y])

# this function gets called by 'Activate_Gaussian_Function'.
# It runs the algorithm and prints (most of) the relevant figures.
# other figures are generated in 'Activate_Gaussian_Function'
def create_fitting(vbroad, period, global_ranking, teff, scatter_title, current_folder_path, ecosw,
                   vbroad_with_outliers='', period_with_outliers='', ecosw_with_outliers=''):
    log_vbroad = np.log10(vbroad)
    log_period = np.log10(period)
    create_folder(current_folder_path)

    backend_filename = f'{current_folder_path}/sampler data.h5'

    if not os.path.exists(backend_filename):
        backend = emcee.backends.HDFBackend(backend_filename)
        np.random.seed(42)

        p0 = np.zeros((nwalkers, ndim))
        for walker_ind in range(nwalkers):
            while True:
                current_walker = [np.random.uniform(low, high) for (low, high) in initial_guess]
                p0[walker_ind, :] = current_walker
                if np.isfinite(log_prob(p0[walker_ind], np.log10(vbroad), np.log10(period))):
                    break

        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob,
                                        args=[np.log10(vbroad), np.log10(period)], backend=backend)
        p0, prob, state = sampler.run_mcmc(p0, n_burnin)
        sampler.reset()
        sampler.run_mcmc(p0, n_samples_MCMC)

    with h5py.File(backend_filename, "r") as backend_file:
        sampler_chain = np.array(backend_file['mcmc']['chain'])

        (median_logx0, median_logy0, median_angle,
         median_std_x, median_std_y, median_p) = \
            (np.median(sampler_chain[:, :, 0]), np.median(sampler_chain[:, :, 1]),
             np.median(sampler_chain[:, :, 2]), np.median(sampler_chain[:, :, 3]),
             np.median(sampler_chain[:, :, 4]), np.median(sampler_chain[:, :, 5]))

        sampler_chain[:, :, 0] = 10 ** sampler_chain[:, :, 0] # you turn log(x0) into x0
        sampler_chain[:, :, 1] = 10 ** sampler_chain[:, :, 1] # same for y0
        sampler_chain[:, :, 2] = -np.tan(sampler_chain[:, :, 2]) # because you want the slope instead of the angle
        q_vec = (1 - 2 * np.pi * sampler_chain[:, :, 3:4] *
                 sampler_chain[:, :, 4:5] * sampler_chain[:, :, 5:6]) / graph_area
        sampler_chain = np.concatenate((sampler_chain, q_vec), axis=2)
        sampler_ln_probability = np.array(backend_file['mcmc']['log_prob'])

    try:
        auto_correlation_length = emcee.autocorr.integrated_time(sampler_chain)
    except AutocorrError:
        auto_correlation_length = 'Not Converging'

    flattened_samples_chain = sampler_chain.reshape(sampler_chain.shape[0]
                                                    * sampler_chain.shape[1], sampler_chain.shape[2])
    sampler_chain_after_convergence = sampler_chain[points_before_convergence:, :, :].copy()
    flattened_samples_chain_after_convergence = sampler_chain_after_convergence.reshape(
        sampler_chain_after_convergence.shape[0] * sampler_chain_after_convergence.shape[1],
        sampler_chain_after_convergence.shape[2]).copy()
    random_rows_indices = np.random.choice(flattened_samples_chain_after_convergence.shape[0], frames_num,
                                           replace=False)

    randomized_1000_list = flattened_samples_chain_after_convergence[random_rows_indices, :]

    slopes = flattened_samples_chain[:, 2]
    intercepts = np.log10(flattened_samples_chain[:, 1]) - slopes * np.log10(flattened_samples_chain[:, 0])
    mean_s, median_s, std_s = round(np.mean(slopes), 3), round(np.median(slopes), 3), round(np.std(slopes), 3)
    mean_i, median_i, std_i = round(np.mean(intercepts), 3), round(np.median(intercepts), 3), round(
        np.std(intercepts), 3)

    logprob_of_median_solution = prob_calc(log_period, log_vbroad, median_logx0, median_logy0, median_angle,
                                           median_std_x, median_std_y, median_p)
    # SI hist:
    fig, axes = plt.subplots(1, 2, figsize=(16, 14))
    axes[0].hist(slopes, bins=np.linspace(median_s - 5 * std_s, median_s + 5 * std_s, 50), density=True,
                 label=f'Median = {median_s}\nSTD = {std_s}')
    axes[0].set_xlabel('Slope', fontsize=20)
    axes[0].tick_params(labelsize=16)
    axes[1].hist(intercepts,
                 bins=np.linspace(median_i - 5 * std_i, median_i + 5 * std_i, 50), density=True,
                 label=f'Median = {median_i}\nSTD = {std_i}')
    axes[1].set_xlabel('Intercept', fontsize=20)
    axes[1].tick_params(labelsize=16)

    for ax in axes.flatten():
        ax.grid(True)
        ax.set_ylabel('Normalized Frequency', fontsize=20)
        ax.legend(loc='upper right', fontsize=18)
    plt.suptitle(scatter_title, fontsize=20, y=0.92)
    plt.savefig(fr'{current_folder_path}\SI hist.jpg')
    plt.close()

    # logprob vs global ranking
    plt.figure(figsize=(14, 12))
    plt.scatter(global_ranking, logprob_of_median_solution, s=10, label=f'N = {len(global_ranking)}')
    plt.ylabel('Log Prob of Median of Solution', fontsize=18)
    plt.xlabel('Global Ranking', fontsize=18)
    plt.tick_params(labelsize=17)
    plt.legend(fontsize=16)
    plt.grid()
    plt.savefig(f'{current_folder_path}/logprob vs g ranking')
    plt.close()

# logprob chain:
    plt.figure(figsize=(16, 14))
    plt.plot(sampler_ln_probability)
    plt.grid()
    plt.ylabel('ln Probability', fontsize=20)
    plt.xlabel('Step', fontsize=20)
    plt.title(f'N Samples = {n_samples_MCMC}\n{auto_correlation_length}', fontsize=20)
    plt.tick_params(labelsize=17)
    plt.savefig(fr'{current_folder_path}\log.jpg')
    plt.close()

# parameters chain:
    fig, axes = plt.subplots(2, 3, figsize=(16, 12))
    plt.suptitle('Parameters vs. Step', fontsize=20)
    for ind in range(nwalkers):
        axes[0, 0].plot(sampler_chain[:, ind, 0])
        axes[0, 1].plot(sampler_chain[:, ind, 1])
        axes[0, 2].plot(-np.tan(sampler_chain[:, :, 2]))
        axes[1, 0].plot(sampler_chain[:, ind, 3])
        axes[1, 1].plot(sampler_chain[:, ind, 4])
        axes[1, 2].plot(sampler_chain[:, ind, 5])

    axes[0, 0].set_title(f'{label_names[0]}')
    axes[0, 1].set_title(f'{label_names[1]}')
    axes[0, 2].set_title(f'{label_names[2]}')
    axes[1, 0].set_title(f'{label_names[3]}')
    axes[1, 1].set_title(f'{label_names[4]}')
    axes[1, 2].set_title(f'{label_names[5]}')

    axes[0, 0].tick_params(labelsize=17)
    axes[0, 1].tick_params(labelsize=17)
    axes[0, 2].tick_params(labelsize=17)
    axes[1, 0].tick_params(labelsize=17)
    axes[1, 1].tick_params(labelsize=17)
    axes[1, 2].tick_params(labelsize=17)

    for col_index in range(int(ndim / 2)):
        axes[0, col_index].tick_params(axis='x', which='both', length=0, labelsize=0)
    for col_index in range(int(ndim / 2)):
        axes[1, col_index].tick_params(axis='x', rotation=45)
    for ax in axes.flatten():
        ax.grid(True)
    plt.subplots_adjust(wspace=0.3)
    plt.savefig(fr'{current_folder_path}\par.jpg')
    plt.close()

# corner plot:
    fig = corner.corner(flattened_samples_chain_after_convergence[:, 0:6],
                        labels=label_names, label_kwargs={"fontsize": 20})
    for ax in fig.get_axes():
        for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
            tick_label.set_rotation(30)
        ax.tick_params(labelsize=14)

    plt.savefig(f'{current_folder_path}/corner plot.jpg')
    plt.close()

    # Histogram of parameters:
    fig, axes = plt.subplots(3, 3, figsize=(16, 12))
    plt.suptitle('Parameters Histogram')
    hist_bins = 50
    axes[0, 0].hist(flattened_samples_chain[:, 0], bins=hist_bins, color="blue", density=True,
                    label=f'{np.median(flattened_samples_chain[:, 0])} +/- {np.std(flattened_samples_chain[:, 0])}')
    axes[0, 1].hist(flattened_samples_chain[:, 1], bins=hist_bins, color="red", density=True,
                    label=f'{np.median(flattened_samples_chain[:, 1])} +/- {np.std(flattened_samples_chain[:, 1])}')
    axes[0, 2].hist(flattened_samples_chain[:, 2], bins=hist_bins, color="green", density=True,
                    label=f'{np.median(flattened_samples_chain[:, 2])} +/- {np.std(flattened_samples_chain[:, 2])}')
    axes[1, 0].hist(flattened_samples_chain[:, 3], bins=hist_bins, color="orange", density=True,
                    label=f'{np.median(flattened_samples_chain[:, 3])} +/- {np.std(flattened_samples_chain[:, 3])}')
    axes[1, 1].hist(flattened_samples_chain[:, 4], bins=hist_bins, color="black", density=True,
                    label=f'{np.median(flattened_samples_chain[:, 4])} +/- {np.std(flattened_samples_chain[:, 4])}')
    axes[1, 2].hist(flattened_samples_chain[:, 5], bins=hist_bins, color="brown", density=True,
                    label=f'{np.median(flattened_samples_chain[:, 5])} +/- {np.std(flattened_samples_chain[:, 5])}')
    axes[2, 1].hist(flattened_samples_chain[:, 6], bins=hist_bins, color="pink", density=True,
                    label=f'{np.median(flattened_samples_chain[:, 6])} +/- {np.std(flattened_samples_chain[:, 6])}')

    axes[0, 0].set_title(f'{label_names[0]}', fontsize=20)
    axes[0, 1].set_title(f'{label_names[1]}', fontsize=20)
    axes[0, 2].set_title(f'{label_names[2]}', fontsize=20)
    axes[1, 0].set_title(f'{label_names[3]}', fontsize=20)
    axes[1, 1].set_title(f'{label_names[4]}', fontsize=20)
    axes[1, 2].set_title(f'{label_names[5]}', fontsize=20)
    axes[2, 1].set_title(f'{label_names[6]}', fontsize=20)

    axes[0, 0].tick_params(labelsize=17)
    axes[0, 1].tick_params(labelsize=17)
    axes[0, 2].tick_params(labelsize=17)
    axes[1, 0].tick_params(labelsize=17)
    axes[1, 1].tick_params(labelsize=17)
    axes[1, 2].tick_params(labelsize=17)
    axes[2, 1].tick_params(labelsize=17)

    for ax in axes.flatten():
        ax.grid(True)
    for (row_ind, col_ind) in [(2, 0), (2, 2)]:
        axes[row_ind, col_ind].remove()
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.savefig(f'{current_folder_path}/par hist')
    plt.close()

    # color coded g ranking:
    single_scatter_plot(log_period, log_vbroad, period_label, vbroad_label,
                        f'{current_folder_path}/color coded for g ranking.jpg', save_format='jpg',
                        z=global_ranking,
                        zlabel=global_ranking_label)

    # color coded temperature:
    single_scatter_plot(log_period, log_vbroad, period_label, vbroad_label,
                        f'{current_folder_path}/color coded for temperature.jpg', save_format='jpg',
                        z=teff,
                        zlabel=temperature_label, title=scatter_title)

    # scatter without solution:
    single_scatter_plot(period, vbroad, period_label, vbroad_label,
                        xlogscale=True, ylogscale=True, savepath=f'{current_folder_path}/scatter without solution.jpg')

    # scatter with solution:
    scatter_plot_with_solution(period, vbroad, randomized_1000_list[:, 0], randomized_1000_list[:, 1],
                                                           randomized_1000_list[:, 2], randomized_1000_list[:, 3], randomized_1000_list[:, 4],
                                                           savepath=current_folder_path,
                                                           title=scatter_title, cosw=ecosw)

    # scatter with solution (removed outliers)
    scatter_with_solution_fig = scatter_plot_with_solution(period_with_outliers, vbroad_with_outliers,
                                                           randomized_1000_list[:, 0],
                                                           randomized_1000_list[:, 1],
                                                           randomized_1000_list[:, 2],
                                                           randomized_1000_list[:, 3],
                                                           randomized_1000_list[:, 4],
                                                           savepath=f'{current_folder_path}/period vbroad with outliers',
                                                           title=scatter_title, markoutliers=True,
                                                           cosw=ecosw_with_outliers)

    return median_s, median_i, scatter_with_solution_fig

# parameters:
initial_guess = [(-0.1, -0.07), (2, 2.03), (0.6, 0.7)
    , (0.25, 0.28), (0.12, 0.14), (4, 5)]
MCMC_parameters_ranges_string = f'{initial_guess}'
period_lower_cutoff, period_upper_cutoff = 0.3, 3
vbroad_lower_cutoff, vbroad_upper_cutoff = 30, 300
nwalkers = 32
n_samples_MCMC = 10 ** 4
n_burnin = 10 ** 3
ndim = 6
label_names = ['X0 [day]', 'y0 [km/s]', 'Slope', 'STD (x)', 'STD (y)', 'p', 'q']
points_before_convergence = int(2.5 * 10 ** 3)
frames_num = 100
graph_area = (np.log10(vbroad_upper_cutoff) - np.log10(vbroad_lower_cutoff)) * (np.log10(period_upper_cutoff) - np.log10(period_lower_cutoff))


