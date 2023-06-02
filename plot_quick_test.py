#!/usr/bin/env python3

import csv
import matplotlib.pyplot as plt
import sys

# Function to load data from a file
def load_data(filename):
    data = []
    with open(filename, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            data.append(row)
    return data

# Function to plot the data
def plot_comparison(data_dict):
    # Plot search times
    fig, ax1 = plt.subplots(figsize=(10, 5))
    for plot_name, file_name in data_dict.items():
        data = load_data(file_name)
        file_names = [row['File'].split('/')[-1].split('-')[0] + '-' + row['File'].split('-')[-2] + '-' + row['File'].split('-')[-1].split('.')[0] for row in data]
        search_times = [float(row['Search Time']) for row in data]
        ax1.plot(file_names, search_times, label=plot_name)

    ax1.set_ylabel('Search Time (seconds)')
    ax1.legend(loc='upper left')

    # Rotate x-axis labels for better visibility
    plt.xticks(rotation=45, ha='right')

    # Add a subtext
    plt.figtext(0.5, 0.02, 'Average taken from 5 runs', ha='center', fontsize=12)

    # Display the chart
    plt.suptitle('Search Time Comparison', fontsize=16)
    plt.show()
    if save:
        plt.savefig('search_time_comparison.png')

    # Plot peak memory
    fig, ax2 = plt.subplots(figsize=(10, 5))
    for plot_name, file_name in data_dict.items():
        data = load_data(file_name)
        file_names = [row['File'].split('/')[-1].split('-')[0] + '-' + row['File'].split('-')[-2] + '-' + row['File'].split('-')[-1].split('.')[0] for row in data]
        peak_memory = [float(row['Peak Memory']) for row in data]
        ax2.plot(file_names, peak_memory, label=plot_name)

    ax2.set_xlabel('File')
    ax2.set_ylabel('Peak Memory (KB)')
    ax2.legend(loc='upper right')

    # Rotate x-axis labels for better visibility
    plt.xticks(rotation=45, ha='right')

    # Add a subtext
    plt.figtext(0.5, 0.02, 'Average taken from 5 runs', ha='center', fontsize=12)

    # Display the chart
    plt.suptitle('Peak Memory Comparison', fontsize=16)
    plt.show()
    if save:
        plt.savefig('peak_memory_comparison.png')

# Define the files to load with their plot names
data_dict = {
    'Classical': 'quick_test_og.csv',
    'Symbolic': 'quick_test_symb.csv',
    'Symbolic (set ordering)': 'quick_test_symb_set_ordering.csv',
}

save = len(sys.argv) > 1 and sys.argv[1] == 'save'
# check if file is called with save option

# Plot the comparison
plot_comparison(data_dict)
