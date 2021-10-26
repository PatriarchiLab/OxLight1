import matplotlib.pyplot as plt
import numpy as np

from typing import NamedTuple, Tuple, List

from tail_pick_data_utils import load_metadata, load_fp_data, get_zscore, load_behaviour_data, extract_data_event_period, \
    get_interpolated_data, get_data_distribution, get_zscore_from_distribution, get_dff_poly, show_p_value, \
    smooth_data, Measurement

__author__ = "Anna Zych"
__maintainer__ = "Anna Zych"
__email__ = "annazych@neuro.mpg.de"

# setting necessary variables for the data (groups and animals used) and figures (fonts and colours)

font = {'family': 'arial',
        'size': 22}
plt.rc('font', **font)

DATA_PATH = r"C:\Users\annazych\Documents\oxLightpy\data"

PALE_GREEN = '#6ABD45'
GREEN = '#0F8140'

ExperimentConfig = NamedTuple('ExperimentConfig', [('sheet_name', str),
                                                   ('behaviour_sheet_name', str),
                                                   ('animals', List[int]),
                                                   ('name', str),
                                                   ('plot_colors', Tuple[str, str])
                                                   ])

EXPERIMENT_CONFIG = ExperimentConfig(sheet_name="all_OxL_pIC",
                                     behaviour_sheet_name="OxLight_pIC_all",
                                     animals=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                                     name="experimental",
                                     plot_colors=('black', PALE_GREEN))
CONTROL_CONFIG = ExperimentConfig(sheet_name="orexin_stimuli_controls_pIC",
                                  behaviour_sheet_name="OxLight_controls_pIC",
                                  animals=[1, 3, 6, 8],
                                  name="control",
                                  plot_colors=('black', GREEN))

def load_data(experiment_config):
    # loading meta data
    meta_file_path = r"{}\meta_data.xlsx".format(DATA_PATH)
    behaviour_file_path = r"{}\Scoring_all.xlsx".format(DATA_PATH)
    animal, data_folder, fp_file, auto, gcamp, time_stamp = load_metadata(meta_file_path,
                                                                          experiment_config.sheet_name)
    data = [{'ani': k} for k in np.arange(animal.size)]

    behaviour = load_behaviour_data(behaviour_file_path,
                                    experiment_config.behaviour_sheet_name)

    # loading raw data from 1 rec at a time

    for animal_idx in np.arange(animal.size, dtype=int):
        fp_times, auto_data, gcamp_data = load_fp_data(data_folder,
                                                       fp_file,
                                                       auto,
                                                       gcamp,
                                                       time_stamp,
                                                       animal_idx,
                                                       time_offset=30)

        # processing OxLight signal using functions from the data_utils

        auto_data = auto_data
        gcamp_data = gcamp_data
        control_fit, dff = get_dff_poly(gcamp_data, fp_times, 2)
        zscore = get_zscore(dff)

        data[animal_idx]['ts'] = fp_times
        data[animal_idx]['auto'] = auto_data
        data[animal_idx]['gcamp'] = gcamp_data
        data[animal_idx]['dff'] = dff
        data[animal_idx]['zscore'] = zscore
        data[animal_idx]['behaviour'] = behaviour[animal_idx]

    return data


# get the data for plotting the Oxlight signal coupled to the tailpick
def plot_episodes(data, animal_slice, behaviour_name, group_name, colors):
    # setting the timeframe for the tail pick (tst) analysis: 5 seconds before the event + 10 second of event duration
    time_period = (-5, 10)
    behaviour_episodes_dff = []
    animal_batch_name = ''
    for animal_idx in animal_slice:
        animal_batch_name += str(animal_idx)

        for episode in data[animal_idx]['behaviour'][behaviour_name]:
            period = (episode[0] + time_period[0],
                      episode[0] + time_period[1])
            period_dff, period_time = extract_data_event_period(data[animal_idx]['dff'],
                                                                data[animal_idx]['ts'],
                                                                period)

            period_time = period_time - period_time[0] + time_period[0]  # all should start at 0
            period_measurement = Measurement(data_series=period_dff, time_series=period_time)
            behaviour_episodes_dff.append(period_measurement)

    # get the statistics for each group and plot graphs for OxLight signal and statistics for both groups

    all_data = []
    pre_data = []
    post_data = []
    for episode in behaviour_episodes_dff:
        time_series = episode.time_series
        data_series = episode.data_series

        data_series_for_distribution = data_series[time_series < 0]
        distribution = get_data_distribution(data_series_for_distribution)
        data_series = get_zscore_from_distribution(data_series, distribution)
        measurement = Measurement(time_series=time_series, data_series=data_series)
        all_data.append(measurement)
        animal_min_pre = np.mean(measurement.data_series[measurement.time_series < 0])
        animal_min_post = np.mean(measurement.data_series[measurement.time_series > 0])
        pre_data.append(animal_min_pre)
        post_data.append(animal_min_post)
    show_p_value(pre_data, post_data, "Final_episode_{}_{}".format(group_name, behaviour_name), colors[1])

    plt.figure(behaviour_name)
    if all_data != []:
        mean_data, std_data, sem, all_time_series = get_interpolated_data(all_data)
        plt.plot(all_time_series, smooth_data(mean_data, window_length=201), color=colors[0])
        plt.fill_between(all_time_series,
                         smooth_data(mean_data - sem, window_length=201),
                         smooth_data(mean_data + sem, window_length=201),
                         alpha=0.3,
                         color=colors[1])
    for episode in all_data:
        plt.plot(episode.time_series,
                 smooth_data(episode.data_series, window_length=201),
                 color='gray',
                 alpha=0.3,
                 linewidth=0.9)

    plt.ylabel("ΔF/F (%)")
    plt.xlabel("Time (s)")
    plt.xlim(time_period)
    plt.axvline(x=0, color='black', linewidth=0.5)
    plt.ylim(-5, 16)
    filepath = r"{}\{}{}animal{}.svg".format(DATA_PATH, group_name, behaviour_name, animal_batch_name)
    plt.savefig(filepath)
    plt.show()


if __name__ == '__main__':
    for config in [CONTROL_CONFIG, EXPERIMENT_CONFIG]:
        data = load_data(config)
        plot_episodes(data, config.animals, 'tst', config.name, config.plot_colors)
