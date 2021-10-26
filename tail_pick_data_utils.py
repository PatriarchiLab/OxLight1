import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pingouin as pt
from scipy.signal import savgol_filter
from collections import defaultdict
from typing import Tuple, NamedTuple

from scipy.stats import stats

Distribution = NamedTuple('Distribution', [('mean', float),
                                           ('std', float)])

Measurement = NamedTuple('Measurement', [('data_series', np.ndarray),
                                         ('time_series', np.ndarray)])


def get_data_distribution(data):
    # type: (np.ndarray) -> Distribution
    mean = np.mean(data)
    std = np.std(data)
    return Distribution(mean=mean, std=std)


def get_zscore_from_distribution(data, distribution):
    # type: (np.ndarray, Distribution) -> np.ndarray
    zscore = (data - distribution.mean) / distribution.std
    return zscore


def load_metadata(meta_file_path, sheet_name):
    # type: (str, str) -> ...
    data_list = pd.read_excel(meta_file_path, sheet_name=sheet_name)

    animal = np.array(data_list.values[:, 0])
    data_folder = np.array(data_list.values[:, 1])
    fp_file = np.array(data_list.values[:, 2])
    auto = np.array(data_list.values[:, 3])
    gcamp = np.array(data_list.values[:, 4])
    time_stamp = np.array(data_list.values[:, 5])

    return animal, data_folder, fp_file, auto, gcamp, time_stamp


def load_fp_data(data_folder, fp_file, auto, gcamp, time_stamp, animal_idx, time_offset=None):
    # load the corresponding fp file (ignore the first raw with text)
    file_name = os.path.join(data_folder[animal_idx], fp_file[animal_idx])
    csv = pd.read_csv(file_name, skiprows=(1))  # , skiprows = (0), skipfooter =(10))
    fp_times = csv.values[:, time_stamp[animal_idx]]
    auto = csv.values[:, auto[animal_idx]]
    gcamp = csv.values[:, gcamp[animal_idx]]

    if time_offset is not None:
        auto = auto[fp_times > time_offset]
        gcamp = gcamp[fp_times > time_offset]
        fp_times = fp_times[fp_times > time_offset]

    gcamp = validate_data(gcamp)
    auto = validate_data(auto)

    return fp_times, auto, gcamp


def my_default_dict():
    return defaultdict(list)


def load_behaviour_data(file_path, sheet_name):
    data_list = pd.read_excel(file_path, sheet_name=sheet_name)
    behaviour_data = defaultdict(my_default_dict)
    counter = {"eating": 0,
               "grooming": 0,
               "tst": 0,
               }
    for row_id in range(data_list.shape[0]):
        animal_id = data_list.values[row_id, 0]
        eating = tuple(data_list.values[row_id, 2:4])
        grooming = tuple(data_list.values[row_id, 4:6])
        tst = tuple(data_list.values[row_id, 6:8])

        if not np.isnan(animal_id):
            counter["eating"] = add_data_to_dict(animal_id, 'eating', behaviour_data, eating, counter["eating"])
            counter["grooming"] = add_data_to_dict(animal_id, 'grooming', behaviour_data, grooming, counter["grooming"])
            counter["tst"] = add_data_to_dict(animal_id, 'tst', behaviour_data, tst, counter["tst"])


    return behaviour_data


def add_data_to_dict(animal_id, key, behaviour_dict_data, data, counter):
    if not np.isnan(data[0]) and not np.isnan(data[1]):
        behaviour_dict_data[animal_id][key].append(data)
        return counter + 1
    return counter


def get_zscore(dff):
    # zscore whole data set with overall mean
    distribution = get_data_distribution(dff)
    dffzscore = get_zscore_from_distribution(dff, distribution)
    return dffzscore


def normalize(data):
    data_min = np.nanmin(data)
    data = data - data_min
    data_max = np.nanmax(data)
    data = data / data_max
    return data


def get_dff(auto, gcamp):
    # fitting like in LERNER paper
    # https://github.com/talialerner
    reg = np.polyfit(auto, gcamp, 1)
    a = reg[0]
    b = reg[1]
    controlFit = a * auto + b
    dff = (gcamp - controlFit) / controlFit
    dff = dff * 100
    return auto, gcamp, controlFit, dff

def get_dff_poly(gcamp, ts, order=2):
    # compute the dff like in Cai et al. (2020), Witten lab
    #
    # dF/F = (F-F0)/ F0, where F0 is the second order polynomial fit to the gcamp signal
    poly_gcamp = np.polyfit(ts, gcamp, order)
    poly_gcamp = np.polyval(poly_gcamp, ts)
    poly_dff = (gcamp - poly_gcamp) / poly_gcamp
    poly_dff = poly_dff * 100
    return poly_gcamp, poly_dff

def get_dff_mean(gcamp):
    # compute the dff like Alex suggested with mean
    mean_gcamp = np.mean(gcamp)
    stand_dff = (gcamp - mean_gcamp) / mean_gcamp
    stand_dff = stand_dff * 100
    control_fit = stand_dff * 0 + mean_gcamp
    return control_fit, stand_dff

def smooth_data(data, window_length=21, polyorder=2):
    # smoothing the data by applying filter
    data = savgol_filter(data, window_length, polyorder)
    return data


def validate_data(data):
    # convert array of object to array of float
    valid_data = data.astype(np.float)

    # replace NaN's with closest non-NaN
    mask = np.isnan(valid_data)
    time_with_nans = np.flatnonzero(mask)
    time_with_data = np.flatnonzero(~mask)
    valid_data[mask] = np.interp(time_with_nans, time_with_data, valid_data[~mask])
    return valid_data


def is_data_valid(data):
    """returns an error when data contains nans"""
    nan_mask = np.isnan(data)
    assert not nan_mask.any()


def extract_data_event_period(data_series, time_series, period):
    # type: (np.ndarray, np.ndarray, Tuple[float, float]) -> Tuple[np.ndarray, np.ndarray]
    condition_1 = time_series >= period[0]
    condition_2 = time_series <= period[1]
    mask = np.logical_and(condition_1, condition_2)
    period_data_series = data_series[mask]
    period_time_series = time_series[mask]

    return period_data_series, period_time_series


def get_interpolated_data(all_data):
    # type: (List[Measurement]) -> List[Tuple[np.ndarray, np.ndarray]]
    all_time_series = all_data[0].time_series
    all_data_series = np.array([all_data[0].data_series])
    for series in all_data[1:]:
        interpolated_data = np.interp(all_time_series, series.time_series, series.data_series)
        all_data_series = np.vstack((all_data_series, np.array([interpolated_data])))
    mean_data = np.mean(all_data_series, axis=0)
    std_data = np.std(all_data_series, axis=0)
    sem = stats.sem(all_data_series, axis=0)

    return mean_data, std_data, sem, all_time_series


def get_data_for_stats(animal_behaviours, animal_slice):
    x = []
    y = []
    each_x = []
    each_y = []
    for animal_idx in animal_slice:
        if animal_behaviours[animal_idx] != []:
            mean_data, std_data, _, all_time_series = get_interpolated_data(animal_behaviours[animal_idx])
            animal_min_pre = np.mean(mean_data[all_time_series < 0])
            animal_min_post = np.mean(mean_data[all_time_series > 0])
            x.append(animal_min_pre)
            y.append(animal_min_post)
            for episode in animal_behaviours[animal_idx]:
                measurement = episode[0]
                episode_pre = np.mean(measurement.data_series[measurement.time_series < 0])
                episode_post = np.mean(measurement.data_series[measurement.time_series > 0])
                each_x.append(episode_pre)
                each_y.append(episode_post)
    return each_x, each_y, x, y

# sample use
# each_x, each_y, x, y = get_data_for_stats(animal_behaviours, animal_slice)
# show_p_value(x, y, "Each animal", behaviour_name, group_name)
# show_p_value(each_x, each_y, "Each_episode_{}".format(group_name), behaviour_name, group_name)


def show_p_value(x, y, name, color):

    ttest_result = pt.ttest(x, y, paired=True)
    print(ttest_result)
    p_val = ttest_result['p-val']["T-test"]
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    x_sem = stats.sem(x)
    y_sem = stats.sem(y)

    print(x_mean, y_mean)

    fig = plt.figure("mean-p_val", figsize=(3, 4.8))
    plt.subplots_adjust(wspace=0.6, hspace=0.6, left=0.2, bottom=0.2, right=0.9, top=0.9)
    # plt.xlabel("Time [s]")

    pre_label = 'a'
    post_label = 'c'
    for x_point, y_point in zip(x, y):
        plt.plot([pre_label, post_label], [x_point, y_point], marker="o", fillstyle='none', color='gray', alpha=0.3)

    plt.plot([pre_label, post_label], [x_mean, y_mean], color=color, alpha=0.6)
    plt.scatter([pre_label, post_label], [x_mean, y_mean], marker="o", color=color)
    plt.plot([pre_label, pre_label], [x_mean - x_sem, x_mean + x_sem], marker="_", color=color)
    plt.plot([post_label, post_label], [y_mean - y_sem, y_mean + y_sem], marker="_", color=color)
    plt.ylim(-2, 12)
    fig.suptitle("p = {:.3f}".format(p_val))
    plt.ylabel("ΔF/F (%)")

    filepath = r"C:\Users\annazych\Documents\oxLightpy\data\{}p-3-value-animal.svg".format(name)
    plt.savefig(filepath)
    plt.show()
