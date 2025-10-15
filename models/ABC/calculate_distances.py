# Note: This script is called as a subprocess from calculate_distances.R

import os
import csv
import math

def euclidean_distance(x, y):
    return math.sqrt(sum((xi - yi) ** 2 for xi, yi in zip(x, y)))

def sum_absolute_difference(x, y):
    return sum(abs(xi - yi) for xi, yi in zip(x, y))

def rmse(x, y):
    return math.sqrt(sum((xi - yi) ** 2 for xi, yi in zip(x, y)) / len(x))

def process_files(directory, file_prefix, data):
    target_prefix = f"{file_prefix}_"
    output_rows = []

    times = []
    prevalences = []
    data_csv_path = data
    with open(data_csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(int(row['time'].strip()))
            prevalences.append(float(row['actual'].strip()))
    
    all_filenames = sorted([f for f in os.listdir(directory) if f.startswith(target_prefix) and f.endswith('.csv')])
        
    
    for file_index, filename in enumerate(all_filenames, 1):
        file_path = os.path.join(directory, filename)
        
        with open(file_path, 'r') as file:
            lines = file.readlines()
            second_line = lines[1].strip()
            values = second_line.split(',')

            a, b, c, d, e, f, g, h, i, j, k, l = values[:12]
            m = ','.join(values[12:])
            
            list_of_trials = [[float(val.strip().strip('"').strip("'")) for val in segment.split(',') if val.strip()] for segment in m.split(';') if segment.strip()]

            list_of_trials = [sublist for sublist in list_of_trials if sublist and sublist[-1] != 0]

            filtered_trials = [[sublist[i] for i in times if i < len(sublist)] for sublist in list_of_trials]
            
            euclidean_vals = []
            abs_diff_vals = []
            rmse_vals = []

            for sublist in filtered_trials:
                if not sublist:
                    continue
            
                euclidean_vals.append(euclidean_distance(sublist, prevalences))
                abs_diff_vals.append(sum_absolute_difference(sublist, prevalences))
                rmse_vals.append(rmse(sublist, prevalences))
                
            if not euclidean_vals:
                continue

            mean_distance = sum(euclidean_vals) / len(euclidean_vals)
            mean_abs_diff = sum(abs_diff_vals) / len(abs_diff_vals)
            mean_rmse = sum(rmse_vals) / len(rmse_vals)

            output_rows.append([a, b, c, d, e, f, g, h, i, j, k, l,
                                mean_distance, mean_abs_diff, mean_rmse])

        print(f"File {file_index} of {len(all_filenames)} processed")

    output_file = os.path.join(directory, f"dist_{file_prefix}.csv")
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['R0_wuhan', 'R0_alpha', 'R0_delta', 'R0_omicron', 
                  'recovery_wuhan', 'recovery_alpha', 'recovery_delta', 'recovery_omicron',
                  'waning_wuhan', 'waning_alpha', 'waning_delta', 'waning_omicron',
                  'mean_distance', 'mean_abs_diff', 'mean_rmse']
        writer.writerow(header)
        writer.writerows(output_rows)

    print(f"\nOutput written to: {output_file}")
