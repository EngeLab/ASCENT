#!/usr/bin/env python3
"""
Find FASTQ files for a list of cells using parallel processing.
Usage: find_normal_fastqs.py <cell_list> <fastq_dir> <output_file> [threads]

Output format:
cell    fq1    fq2
Multiple rows will be created for cells with multiple matching FASTQ pairs.
"""

import sys
import os
import subprocess
from pathlib import Path
import multiprocessing as mp
from threading import Lock
import time

# Global variables for progress tracking
processed_cells = 0
found_pairs = 0
progress_lock = Lock()

def find_fastqs_for_cell(args):
    """Find R1 and R2 FASTQ files for a cell using find command"""
    cell_id, fastq_dir = args
    patterns = [
        f"{cell_id}_S*_R1_001.fastq.gz",  # Standard Illumina naming
        f"{cell_id}_*_1.fq.gz"            # Alternative naming
    ]
    
    found_pairs = []
    
    for pattern in patterns:
        cmd = f"find {fastq_dir} -type f -name '{pattern}'"
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            r1_files = [f for f in result.stdout.strip().split('\n') if f]
            
            if not r1_files:
                continue
            
            # Get corresponding R2 files and check they exist
            for r1_file in r1_files:
                r2_file = r1_file.replace("_R1_", "_R2_").replace("_1.fq", "_2.fq")
                if os.path.exists(r2_file):
                    found_pairs.append((r1_file, r2_file))
                
        except subprocess.SubprocessError:
            continue
            
    return cell_id, found_pairs

def update_progress(result, total_cells):
    """Update progress counters"""
    global processed_cells, found_pairs
    cell_id, pairs = result
    
    with progress_lock:
        processed_cells += 1
        if pairs:
            found_pairs += len(pairs)
        
        # Print progress every 10 cells or at completion
        if processed_cells % 10 == 0 or processed_cells == total_cells:
            print(f"Progress: {processed_cells}/{total_cells} cells processed ({found_pairs} FASTQ pairs found)", 
                  end='\r', flush=True)

def process_results(result, output_file):
    """Write results to output file - one row per FASTQ pair"""
    cell_id, pairs = result
    if pairs:
        with open(output_file, 'a') as out:
            for r1_file, r2_file in pairs:
                out.write(f"{cell_id}\t{r1_file}\t{r2_file}\n")

def main():
    if len(sys.argv) not in [4, 5]:
        print(__doc__)
        sys.exit(1)
        
    cell_list_file = sys.argv[1]
    fastq_dir = sys.argv[2]
    output_file = sys.argv[3]
    num_threads = int(sys.argv[4]) if len(sys.argv) == 5 else mp.cpu_count()
    
    # Read cell IDs
    with open(cell_list_file) as f:
        cells = [line.strip() for line in f if line.strip()]
    
    total_cells = len(cells)
    print(f"Processing {total_cells} cells using {num_threads} threads...")
    
    # Write header to output file
    with open(output_file, 'w') as out:
        out.write("cell\tfq1\tfq2\n")
    
    # Prepare arguments for parallel processing
    args = [(cell, fastq_dir) for cell in cells]
    
    # Process cells in parallel
    start_time = time.time()
    with mp.Pool(num_threads) as pool:
        for result in pool.imap_unordered(find_fastqs_for_cell, args):
            update_progress(result, total_cells)
            process_results(result, output_file)

    # Final summary
    elapsed_time = time.time() - start_time
    print(f"\nCompleted: {total_cells}/{total_cells} cells processed")
    print(f"Found {found_pairs} FASTQ pairs across {total_cells} cells")
    print(f"Time elapsed: {elapsed_time:.1f} seconds")
    print(f"Average time per cell: {(elapsed_time/total_cells):.2f} seconds")

if __name__ == '__main__':
    main()