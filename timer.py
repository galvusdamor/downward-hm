#! /usr/bin/env python3

import subprocess
import glob
import csv
import time

# Global variables
timeout = 30  # Timeout in seconds
file_path = "classical-domains/classical-sas/blocks-probBLOCKS*.sas"
# file_path = "classical-domains/classical-sas/blocks-probBLOCKS-6-0.sas"
search_heuristic = "hm(2)"
# search_heuristic = "symb_h2()"
csv_file_name = "quick_test_class.csv"
num_runs = 5  # Number of test runs

def run_command_for_files():
    # build
    process = subprocess.Popen("./build.py", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    # wait for the process to complete
    process.wait()

    files = glob.glob(file_path)
    # sort files first on the name of the problem then on the number of blocks
    files = sorted(files, key=lambda x: ''.join(x.split('-')[:-2]) + x.split('-')[-2].zfill(3) + x.split('-')[-1].split('.')[0].zfill(3))
    total_files = len(files)
    estimated_runtime = (total_files * num_runs * timeout) / 3600  # Estimation in hours
    
    print(f"Estimated runtime: {estimated_runtime} hours")
    print()
    
    processed_files = []
    
    # Check if CSV file exists and load processed files
    if csv_file_exists():
        processed_files = load_processed_files()
        print(f"Resuming from file {len(processed_files) + 1}/{total_files}")
    

    problem_to_skip = ""
    for idx, file in enumerate(files, 1):
        if file in processed_files:
            continue
        if ''.join(file.split('-')[:-2]) == problem_to_skip:
            continue
        
        print(f"Processing file {idx}/{total_files}: {file}")
        
        search_times = []
        total_times = []
        peak_memories = []
        elapsed_times = []
        
        for _ in range(num_runs):
            # Prepare the command
            command = f"./fast-downward.py {file} --search \"astar({search_heuristic})\""
            
            try:
                # Run the command with a timeout
                process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
                
                try:
                    # Wait for the process to complete or timeout
                    start_time = time.time()
                    output, _ = process.communicate(timeout=timeout)
                    end_time = time.time()
                    elapsed_time = end_time - start_time
                    
                    # Extract relevant information from the output
                    search_time = extract_value_from_output(output, "Search time:")
                    total_time = extract_value_from_output(output, "Total time:")
                    peak_memory = extract_value_from_output(output, "Peak memory:")
                    
                    # Add the results to the lists
                    search_times.append(float(search_time[:-1]))
                    total_times.append(float(total_time[:-1]))
                    peak_memories.append(float(peak_memory[:-3]))
                    elapsed_times.append(elapsed_time)
                    
                except subprocess.TimeoutExpired:
                    process.kill()
                    print(f"Search for file {file} timed out and was skipped.")
                    problem_to_skip = ''.join(file.split('-')[:-2])
                    break
                    
                except Exception as e:
                    print(f"An error occurred while running the command for file {file}: {e}")
                
            except subprocess.CalledProcessError as e:
                print(f"An error occurred while running the command: {e}")
        
        if len(search_times) == 0:
            continue

        # Calculate averages
        avg_search_time = sum(search_times) / len(search_times)
        avg_total_time = sum(total_times) / len(total_times)
        avg_peak_memory = sum(peak_memories) / len(peak_memories)
        avg_elapsed_time = sum(elapsed_times) / len(elapsed_times)
        
        # Print the average results
        print(f"Avg. Search time: {avg_search_time}")
        print(f"Avg. Total time: {avg_total_time}")
        print(f"Avg. Peak memory: {avg_peak_memory}")
        print(f"Avg. Elapsed time: {avg_elapsed_time} seconds")
        print()
        
        # Save the average results to a CSV file
        save_to_csv(file, avg_search_time, avg_total_time, avg_peak_memory, avg_elapsed_time)
        
        # Add the processed file to the list
        processed_files.append(file)
    
    print("Processing completed.")

def extract_value_from_output(output, keyword):
    value = None
    lines = output.split('\n')
    
    for line in lines:
        if keyword in line:
            value = line.split(keyword)[1].strip()
    
    return value

def save_to_csv(file, search_time, total_time, peak_memory, elapsed_time):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    data = [timestamp, file, search_time, total_time, peak_memory, elapsed_time]
    columns = ['Timestamp', 'File', 'Search Time', 'Total Time', 'Peak Memory', 'Elapsed Time']
    
    with open(csv_file_name, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write column names if the file is empty
        if csvfile.tell() == 0:
            writer.writerow(columns)
        
        writer.writerow(data)

def csv_file_exists():
    try:
        with open(csv_file_name, 'r') as csvfile:
            return True
    except FileNotFoundError:
        return False

def load_processed_files():
    processed_files = []
    
    with open(csv_file_name, 'r') as csvfile:
        reader = csv.reader(csvfile)
        
        # Skip the header row
        next(reader)
        
        for row in reader:
            file = row[1]
            processed_files.append(file)
    
    return processed_files



# Run the command for each file
run_command_for_files()
