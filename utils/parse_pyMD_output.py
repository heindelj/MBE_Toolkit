import numpy as np
import matplotlib.pyplot as plt
import sys

try:
    ifile = sys.argv[1]
except:
    print("Did not get a pyMD output file. Please provide this as an argument.")
    sys.exit(1)

def parse_md_file_to_dict(ifile: str):
    """
    Takes a pyMD output md file and returns a dictionary where each key is labeled
    by the column labels printed at the top of the file.
    """
    md_data = {}
    num_header_lines = 0
    with open(ifile, 'r') as f:
        reading_header = True
        while reading_header:
            line = f.readline().split()
            if line and line[0] == "Column":
                md_data[line[-1]] = []
                num_header_lines += 1
            else:
                reading_header = False
    data = np.loadtxt(ifile, skiprows=num_header_lines)
    for i, key in enumerate(md_data.keys()):
        md_data[key] = data[:,i]
    return md_data

def get_running_average(data, window_width):
    cumsum_vec = np.cumsum(np.insert(data, 0, 0))
    return (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width

def plot_md_data(md_data: dict, y_data_key: str, x_axis="step", running_average=True, x_label=None, y_label=None):
    """
    Plots the md data collected from parse_md_file_to_dict using md_data[y_data_key]
    as the y-axis data and md_data[x_axis] as the x-axis data.
    By default x_axis="step". The other logical choice would be "time". If
    x_axis is not found in the dictionary, then a simple count will be used as a fallback.
    """
    if y_data_key in md_data:
        y_data = md_data[y_data_key]
        x_data = np.arange(len(md_data[y_data_key]))
        if running_average:
            y_data = get_running_average(y_data, 5000)
        if x_axis in md_data:
            x_data = md_data[x_axis]


        plt.plot(x_data[(len(x_data)-len(y_data)):], y_data)
        plt.hlines(np.mean(md_data[y_data_key]), min(md_data[x_axis]), max(md_data[x_axis]), colors='black')
        if x_label:
            plt.xlabel(x_label)
        else:
            plt.xlabel(f"{x_axis}")
        if y_label:
            plt.ylabel(y_label)
        else:
            plt.ylabel(f"{y_data_key}")
        plt.show()
    else:
        print(f"{y_data_key} was not found in the dictionary. Not plotting anything.")

def get_averages(md_data: dict, lag=0):
    """
    Compues the averages of all quantities (except "step" and "time") from the
    md_data dictionary. Starts the average from index lag.
    """
    keys_to_exclude = ["step", "time"]
    print("Observable Averages:")
    for key in md_data.keys():
        if key not in keys_to_exclude:
            print(f"{key}: {np.mean(md_data[key][lag:])}")

if __name__ == '__main__':
    md_data = parse_md_file_to_dict(ifile)
    get_averages(md_data)
    plot_md_data(md_data, "temperature", 
                x_label="Step Number", 
                y_label="Kinetic Energy (a.u.)")