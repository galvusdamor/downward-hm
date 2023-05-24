#! /usr/bin/env python3

import subprocess
import glob
import csv
import time

# Global variables
timeout = 10  # Timeout in seconds
file_path = "classical-domains/classical-sas/blocks-probBLOCKS*.sas"
csv_file_name = "search_results.csv"

def run_command_for_files():
    files = glob.glob(file_path)
    files.sort()  # Sort files alphabetically
    total_files = len(files)
    estimated_runtime = (total_files * timeout) / 3600  # Estimation in hours
    
    print(f"Estimated runtime: {estimated_runtime} hours")
    print()
    
    processed_files = []
    
    # Check if CSV file exists and load processed files
    if csv_file_exists():
        processed_files = load_processed_files()
        print(f"Resuming from file {len(processed_files) + 1}/{total_files}")
    
    for idx, file in enumerate(files, 1):
        if file in processed_files:
            continue
        
        print(f"Processing file {idx}/{total_files}: {file}")
        
        # Prepare the command
        command = f"./build.py && ./fast-downward.py {file} --search \"astar(symb_h2())\""
        
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
                
                # Print the extracted information
                print(f"Search time: {search_time}")
                print(f"Total time: {total_time}")
                print(f"Peak memory: {peak_memory}")
                print(f"Elapsed time: {elapsed_time} seconds")
                print()
                
                # Save the information to a CSV file
                save_to_csv(file, search_time, total_time, peak_memory, elapsed_time)
                
                # Add the processed file to the list
                processed_files.append(file)
                
            except subprocess.TimeoutExpired:
                process.kill()
                print(f"Search for file {file} timed out and was skipped.")
                
            except Exception as e:
                print(f"An error occurred while running the command for file {file}: {e}")
            
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running the command: {e}")
    
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
