import uproot

def print_histogram_bins(file_name, hist_name):
    with uproot.open(file_name) as root_file:
        hist = root_file[hist_name]

        bin_contents = hist.values()

        x_axis = hist.member("fXaxis")
        if hasattr(x_axis, "labels") and callable(x_axis.labels):
            bin_labels = x_axis.labels()
        else:
            bin_labels = [f"Bin{i}" for i in range(len(bin_contents))]

        print(f"Histogram: {hist_name}")
        for i, (content, label) in enumerate(zip(bin_contents, bin_labels)):
            print(f"Bin: {i+1} : {label} : {content}")

file_name = "skimmed_nano.root"
hist_name = "cutFlow"
print_histogram_bins(file_name, hist_name)
