# AutoGapfiller.py - Automated gap filling pipeline for DEGAP v2.0
# Accepted parameters:
# 1.1, --reads: HiFi reads file path
# 1.2, --genome | -g: Genome file path to process
# 1.3, --out | -o: Output directory path
# 1.4, --flag: Gap filling direction (left, right, or all), default: left

import argparse
import os
import sys
import subprocess
import time
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import json

# Add exit function to ensure complete program termination
def _exit_program(code=0):
    """Force program exit, ensuring all threads are terminated"""
    print("\nForcing program exit...")
    os._exit(code)  # Use os._exit instead of sys.exit to ensure all threads terminate

# Create argument parser
parser = argparse.ArgumentParser(description='Gap filling script')

# Add parameter definitions
parser.add_argument('--reads', type=str, help='HiFi reads file path')
parser.add_argument('-g', '--genome', required=True, help='Genome file path')
parser.add_argument('-o', '--out', required=True,default='output', help='Output directory path')
parser.add_argument('--flag', type=str, default='all', help='Choose seqleft or seqright as seed to fill the gap: left, right, or all (both)')
parser.add_argument('--kmer_size', '-ks', type=int, default=41, help='k-mer size for filtering reads')
parser.add_argument('--kmer_num', '-kn', type=int, default=10, help='number of k-mers to use for filtering reads')
parser.add_argument('--kmer_len', '-kl', type=float, default=0.1,
                   help='proportion of mean read length for k-mer extraction (0.1 = 10%% of mean length)')
parser.add_argument('--mode', type=str, help='Operation mode: gapfiller or ctglinker')
parser.add_argument('-j', type=int, default=100, help='Number of parallel jobs for processing reads')
parser.add_argument('-t', '--thread', type=str, help='Number of threads')
parser.add_argument('--remove', type=int, default=2,
	                   help='1: only keep final result; 2: keep every round basic result; 3: keep all files')
parser.add_argument('--edge', type=int, default=500, help='Edge Controller set max Edge length')
parser.add_argument('--filterDepth', type=float, help='Filter HiFi reads by mapped depth')
parser.add_argument('--MaximumExtensionLength', type=int, default=1000000, help='Stop Extension when reach this length (default: 1000000 bp)')
parser.add_argument('--batch', type=int, default=None, help='Number of batches for task execution, default None (no batching, all jobs submitted in parallel)')
parser.add_argument('--post_process', action='store_true', help='Execute only post-processing steps for dependent jobs')
parser.add_argument('-q', '--queue', type=str, default='normal', help='Specify queue for gap job submission, default is normal')




# Parse arguments
args = parser.parse_args()


def usage():
	print ("--reads HiFi_reads.fasta")
	print ("-o | --out ./path/")
	print ("--remove 1 | 2 | 3 default:2 1:only keep final result; 2: keep every round basic result ; 3 : keep all files")
	print ("--edge Edge Controller set max Edge length(missequening)")
	print ("--filterDepth num default:None. You can filtered HiFi reads by mapped depth on contig set. if num==0.3 means: mapped Hifi reads on depth>=0.3*avgdepth and depth<=(2-0.2)*avgdepth will be filtered and will not be used in whole project")
	print ("--MaximumExtensionLength num default:1000000. Stop Extension when reach the num (in bp)")
	print ("--kmer_size | -ks num default:41. k-mer size for filtering reads")
	print ("--kmer_num | -kn num default:10. number of k-mers to use for filtering reads")
	print ("--kmer_length | -kl num default:0.1. proportion of mean read length for k-mer extraction (0.1 = 10% of mean length)")
	print ("-j num default:70. number of parallel jobs for processing reads")
	print ("--mode gapfiller | ctglinker")
	print ("--resume num Resume from specified round (e.g.: --resume 118 continues from round118)")
	print ("--resume_auto Automatically resume from last interrupted round")
	print ("\n\ngapfiller\n")
	print ("\t--seqleft sequence before GAP")
	print ("\t--seqright sequence after GAP")
	print ("\t--flag left | right | all Choose filling direction: left uses only left sequence, right uses only right sequence, all uses both directions and prioritizes left direction results (default: all)")


# Reuse core logic from getIndex.py
def process_hifi(reads_path, output_dir):
    """
    Create reads index and generate statistics (with intelligent skip mechanism)
    """
    import math
    from Bio import SeqIO
    
    idx_path = os.path.join(output_dir, "reads.idx")
    stats_path = os.path.join(output_dir, "HiFi.reads.stat")

    # Existence check logic
    if os.path.exists(idx_path) and \
       os.path.exists(stats_path) and \
       os.path.getsize(stats_path) > 0:

        print("Detected existing index and statistics files, skipping processing")

        # Read existing statistics
        stats = {}
        with open(stats_path) as f:
            for line in f:
                key, value = line.strip().split('\t')
                try:
                    stats[key] = float(value)
                except ValueError:
                    stats[key] = value

        # Ensure index file is valid
        try:
            SeqIO.index_db(idx_path)
            return stats
        except Exception as e:
            print(f"Index file corrupted, regenerating: {e}")

    # Continue processing if not exists
    print("\n=== Starting HiFi file processing ===")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create index
    idx_path = os.path.join(output_dir, "reads.idx")
    reads_path_abs = os.path.abspath(reads_path)  # New: get absolute path
    print("Creating reads index file:", idx_path)
    reads_dict = SeqIO.index_db(idx_path, reads_path_abs, "fasta")




    # Generate statistics
    stats_path = os.path.join(output_dir, "HiFi.reads.stat")
    print("Generating statistics file:", stats_path)

    total_len = 0
    max_len = 0
    read_count = len(reads_dict)

    # Use created index for statistics
    print(f"Starting precise statistics calculation for all {read_count} reads...")
    count = 0

    with open(stats_path, 'w') as f:
        # Calculate detailed statistics
        for read_id in reads_dict:
            seq_len = len(reads_dict[read_id].seq)
            total_len += seq_len
            if seq_len > max_len:
                max_len = seq_len
            count += 1
            if count % 100000 == 0:
                print(f"Processed {count}/{read_count} reads ({count/read_count*100:.1f}%)...")

        mean_len = total_len / read_count if read_count > 0 else 0
        print(f"Statistics complete: average read length {mean_len:.2f} bp, maximum read length {max_len} bp")

        # Calculate seed length - use same method as demo.py
        a = 10**(int(math.log(max_len, 10)))
        b = max_len / a + 1
        seed_len = a * b

        # Round to integers
        mean_len_int = int(round(mean_len))
        seed_len_int = int(round(seed_len))

        # Write statistics file - keep original typo TolalLenth for compatibility, but use rounded values
        f.write(f"Number\t{read_count}\n"
                f"TolalLenth\t{total_len}\n"
                f"MaxLength\t{max_len}\n"
                f"MeanLength\t{mean_len_int}\n"
                f"SeedLength\t{seed_len_int}\n")

    return {
        'Number': read_count,
        'TotalLength': total_len,
        'MaxLength': max_len,
        'MeanLength': mean_len_int,
        'SeedLength': seed_len_int
    }

# Execute processing
stats = process_hifi(args.reads, args.out)
print("\nStatistics results:")
for k, v in stats.items():
    print(f"{k}: {v}")

# Part 3: Gap information extraction


def extract_gap_info(genome_path, output_dir, SeedLength):
    # Create output directory and subdirectories
    os.makedirs(output_dir, exist_ok=True)
    gap_data_dir = os.path.join(output_dir, 'gapDataBase')
    os.makedirs(gap_data_dir, exist_ok=True)

    # Open log file for writing
    log_path = os.path.join(output_dir, 'gapbase.log')
    with open(log_path, 'w') as log_file:
        log_file.write("Chromosome\tId\tStart\tEnd\tGapLen\tLeftLen\tRightLen\n")

        # Iterate through each chromosome sequence
        for seq_record in SeqIO.parse(genome_path, "fasta"):
            chrom = seq_record.id
            seq = seq_record.seq
            gaps = []
            in_gap = False
            gap_start = 0

            # Detect all N regions (gaps)
            for i, c in enumerate(seq):
                if c.upper() == 'N':
                    if not in_gap:
                        in_gap = True
                        gap_start = i
                else:
                    if in_gap:
                        in_gap = False
                        gaps.append((gap_start, i - 1))
            # Handle possible gap at the end
            if in_gap:
                gaps.append((gap_start, len(seq) - 1))

            # Process each gap
            for gap_id, (start, end) in enumerate(gaps, 1):
                # Convert to chromosome coordinates (1-based)
                chr_start = start + 1
                chr_end = end + 1
                gap_len = end - start + 1

                # Get left flanking sequence
                left_chars = []
                current_pos = start - 1
                while current_pos >= 0 and len(left_chars) < SeedLength:
                    char = seq[current_pos].upper()
                    if char == 'N':
                        break
                    left_chars.append(char)
                    current_pos -= 1
                left_seq = ''.join(reversed(left_chars))
                left_len = len(left_seq)

                # Get right flanking sequence
                right_chars = []
                current_pos = end + 1
                while current_pos < len(seq) and len(right_chars) < SeedLength:
                    char = seq[current_pos].upper()
                    if char == 'N':
                        break
                    right_chars.append(char)
                    current_pos += 1
                right_seq = ''.join(right_chars)
                right_len = len(right_seq)

                # Write to log
                log_file.write(
                    f"{chrom}\t{gap_id}\t{chr_start}\t{chr_end}\t{gap_len}\t"
                    f"{left_len}\t{right_len}\n"
                )

                # Save sequence files
                if left_len > 0:
                    left_path = os.path.join(
                        gap_data_dir, f"{chrom}.{gap_id}.left.fasta")
                    with open(left_path, 'w') as f:
                        f.write(f">{chrom}.{gap_id}.left\n{left_seq}\n")
                if right_len > 0:
                    right_path = os.path.join(
                        gap_data_dir, f"{chrom}.{gap_id}.right.fasta")
                    with open(right_path, 'w') as f:
                        f.write(f">{chrom}.{gap_id}.right\n{right_seq}\n")


# Part 4: Gap visualization

def get_chromosome_lengths(genome_file):
    """Get chromosome lengths from genome file"""
    chromosome_lengths = {}
    try:
        for record in SeqIO.parse(genome_file, "fasta"):
            chromosome_lengths[record.id] = len(record.seq)
    except Exception as e:
        print(f"Error reading genome file: {e}")
        sys.exit(1)
    return chromosome_lengths

def parse_gap_log(log_file):
    """Parse gap log file and return gap regions for each chromosome"""
    chromosome_gaps = {}
    try:
        with open(log_file, 'r') as f:
            next(f)  # Skip header line
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 5:  # Ensure at least 5 fields (chromosome, ID, start, end, length)
                    continue
                # Parse coordinates correctly from start and end fields (columns 3 and 4)
                chrom, gap_id, start, end = parts[0], parts[1], int(parts[2]), int(parts[3])
                if chrom not in chromosome_gaps:
                    chromosome_gaps[chrom] = []
                chromosome_gaps[chrom].append((start, end))
    except Exception as e:
        print(f"Error reading log file: {e}")
        sys.exit(1)
    return chromosome_gaps

def generate_window_data(chrom_length, gaps, window_size=10000):
    """Generate binary data for each window"""
    num_windows = (chrom_length + window_size - 1) // window_size
    window_data = [0] * num_windows

    for gap_start, gap_end in gaps:
        # Convert to window indices, ensuring correct calculation using window size
        start_window = (gap_start - 1) // window_size  # Convert to 0-based and calculate window index
        end_window = min(num_windows-1, (gap_end - 1) // window_size)

        # Mark all windows containing gaps
        for i in range(start_window, end_window+1):
            window_data[i] = 1

    return window_data

def write_plot_data(output_dir, chromosome_lengths, chromosome_gaps, is_after_filling=False):
    """Write plotting data"""
    # Ensure output directory exists
    plot_dir = os.path.join(output_dir, 'gap_plots')
    os.makedirs(plot_dir, exist_ok=True)

    # Set output file path, distinguishing before/after filling
    file_suffix = "after" if is_after_filling else "before"
    output_file = os.path.join(plot_dir, f'plotData_{file_suffix}.txt')

    try:
        # Check if chromosome and gap data match
        if not chromosome_lengths:
            print("Warning: No chromosome length data found")
            return None

        if not chromosome_gaps:
            print("Warning: No gap data found")

        print(f"Starting to generate {file_suffix} plotting data for {len(chromosome_lengths)} chromosomes, {sum(len(gaps) for gaps in chromosome_gaps.values())} gap regions")

        with open(output_file, 'w') as f:
            for chrom, length in chromosome_lengths.items():
                gaps = chromosome_gaps.get(chrom, [])
                window_data = generate_window_data(length, gaps)

                # Ensure data validity
                if not window_data:
                    print(f"Warning: Window data for chromosome {chrom} is empty, skipped")
                    continue

                data_str = ''.join(map(str, window_data))
                f.write(f"{chrom}\t{length}\t{data_str}\n")
                print(f"Processed chromosome {chrom}, length {length} bp, contains {len(gaps)} gap regions, generated {len(window_data)} windows")

        print(f"Plotting data written to: {output_file}")
    except Exception as e:
        print(f"Error writing plotting data: {e}")
        import traceback
        traceback.print_exc()
        return None

    return output_file

def generate_plot(input_file, output_dir, is_after_filling=False):
    """Generate gap distribution plots for each chromosome based on plotData.txt"""
    # Ensure output directory exists
    plot_dir = os.path.join(output_dir, 'gap_plots')
    os.makedirs(plot_dir, exist_ok=True)

    # Record number of processed chromosomes
    processed_count = 0

    # Distinguish file suffix for before/after filling
    file_suffix = "after" if is_after_filling else "before"

    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()

        if not lines:
            print(f"Warning: Input file {input_file} is empty")
            return False

        print(f"Starting to plot {file_suffix} gap distribution for {len(lines)} chromosomes")

        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                print(f"Warning: Skipping invalid line: {line.strip()}")
                continue

            chrom, length_str, data_str = parts

            try:
                length = int(length_str)
                window_data = [int(x) for x in data_str]
            except ValueError as e:
                print(f"Warning: Data parsing error for chromosome {chrom}: {e}")
                continue

            if not window_data:
                print(f"Warning: Window data for chromosome {chrom} is empty")
                continue
                
            num_windows = len(window_data)
            
            # Set fixed segment length to 10000
            seg_size = 10000
            num_segments = (num_windows + seg_size - 1) // seg_size  # Round up to ensure all data is included

            segments = []
            starts = []

            # Segment by fixed length of 10000
            for i in range(num_segments):
                start = i * seg_size
                end = min((i + 1) * seg_size, num_windows)  # Ensure not exceeding data range
                segments.append(window_data[start:end])
                # Store 0-based indices, but convert to 1-based when displaying
                starts.append(start)

            # Determine maximum segment length for alignment
            max_seg_len = seg_size  # Most segments should be seg_size length, except the last one may be shorter

            plt.figure(figsize=(15, num_segments * 0.6))  # Dynamically adjust figure height based on number of segments

            for seg_idx in range(num_segments):
                seg_data = segments[seg_idx]
                seg_length = len(seg_data)

                # Modify Y coordinate calculation, place first segment at the top
                # Display from top to bottom, first segment (seg_idx=0) corresponds to maximum Y value
                y_coords = [num_segments - seg_idx] * seg_length

                # Calculate x coordinates (all segments start from x=0, normalized by maximum segment length)
                x_coords = [i / max_seg_len for i in range(seg_length)]

                # Plot data points
                for value, color in [(0, '#90EE90'), (1, '#FF3333')]:  # Light green and red
                    points = [i for i, v in enumerate(seg_data) if v == value]
                    if points:
                        plt.scatter(
                            [x_coords[i] for i in points],
                            [y_coords[i] for i in points],
                            s=5,
                            c=color,
                            marker='s',
                            edgecolors='none'
                        )

                # Draw connecting lines
                plt.plot(x_coords, y_coords,
                        color='#607D8B',
                        lw=0.8,
                        alpha=0.3,
                        solid_capstyle='round')

                # Mark start position of each segment - convert 0-based to 1-based, add bp unit
                start_pos = starts[seg_idx] * seg_size + 1
                plt.text(0, y_coords[0] + 0.05, f"{start_pos:,} bp",
                         fontsize=6, color='#607D8B', ha='right')

                # Mark end position of each segment - keep 1-based, add bp unit
                end_pos = min((starts[seg_idx] + seg_length) * seg_size, length)
                plt.text(seg_length / max_seg_len, y_coords[0] + 0.05,
                         f"{end_pos:,} bp", fontsize=6, color='#607D8B')
            
            # Visualization decoration
            title_suffix = "After Filling" if is_after_filling else "Before Filling"
            plt.title(f"{chrom} Gap Distribution ({title_suffix})", pad=15)
            plt.xlim(-0.05, 1.05)
            plt.ylim(0.5, num_segments + 0.5)  # Adjust Y-axis range based on number of segments

            # Hide default axes, keep bottom border for reference
            plt.axis('off')
            plt.axhline(y=0.5, color='#B0BEC5', lw=0.5, alpha=0.5)

            # Add segment marking lines
            for x in [i / max_seg_len for i in range(0, max_seg_len + 1, max_seg_len // 10)]:
                plt.axvline(x, color='#B0BEC5', lw=0.3, alpha=0.2, linestyle='--')

            # Legend settings
            legend_elements = [
                plt.Line2D([0], [0], marker='s', color='w',
                          markerfacecolor='#90EE90', markersize=8,
                          label='Continuous Sequence'),
                plt.Line2D([0], [0], marker='s', color='w',
                          markerfacecolor='#FF3333', markersize=8,
                          label='Gap Region')
            ]
            plt.legend(
                handles=legend_elements,
                loc='upper right',
                frameon=False,
                bbox_to_anchor=(1, 1.15),
                handletextpad=0.3
            )

            # Save image to specified directory, distinguishing before/after filling
            output_path = os.path.join(plot_dir, f"{chrom}_gap_{file_suffix}.png")
            plt.savefig(
                output_path,
                dpi=300,
                bbox_inches='tight',
                facecolor='white'
            )
            plt.close()
            processed_count += 1
            print(f"Generated {file_suffix} gap distribution plot for chromosome {chrom}, {num_segments} segments total")

        print(f"Plotting complete, processed {processed_count} chromosomes")
        return True
    except Exception as e:
        print(f"Error during plotting: {e}")
        import traceback
        traceback.print_exc()
        return False

def plot(genome_file, log_file, output_dir, is_after_filling=False):
    # Determine whether processing before or after filling data
    stage = "after filling" if is_after_filling else "before filling"
    print(f"Starting to process {stage} files: {genome_file} and {log_file}")

    # Ensure using absolute paths
    genome_file = os.path.abspath(genome_file)
    log_file = os.path.abspath(log_file)
    output_dir = os.path.abspath(output_dir)

    # Data preprocessing
    chr_lens = get_chromosome_lengths(genome_file)
    if not chr_lens:
        print("Error: Unable to get chromosome lengths from genome file")
        return False

    chr_gaps = parse_gap_log(log_file)
    if not chr_gaps:
        print(f"Warning: No valid gap data found in log file {log_file}")

    # Generate intermediate data, distinguishing before/after filling
    plot_data_file = write_plot_data(output_dir, chr_lens, chr_gaps, is_after_filling)
    if not plot_data_file:
        print("Error: Failed to generate plotting data")
        return False

    # Visualization, distinguishing before/after filling
    try:
        generate_plot(plot_data_file, output_dir, is_after_filling)
        return True
    except Exception as e:
        print(f"Error during plotting: {e}")
        import traceback
        traceback.print_exc()
        return False






# Part 5: Filter original gap information from Part 3
# Filtering logic:
# 1. Merge gaps with distance < SeedLength, use LeftLen from first gap and RightLen from last gap
# 2. Discard gaps at chromosome ends where either LeftLen or RightLen < SeedLength
# Filtered gap information saved to output_dir/gap.log, filtered gap data saved to output_dir/gapdata/
import os
from Bio import SeqIO

def filter_gaps(output_dir, SeedLength):
    # Create output directory structure
    filtered_dir = os.path.join(output_dir, 'gapData')
    os.makedirs(filtered_dir, exist_ok=True)

    # Read original gap information
    original_log = os.path.join(output_dir, 'gapbase.log')
    gaps_dict = {}
    with open(original_log) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            gap_len = int(parts[4])

            # Modified filtering condition: keep gaps with >= 20 Ns
            if gap_len < 20:
                continue
                
            gap_info = {
                'id': int(parts[1]),
                'start': int(parts[2]),
                'end': int(parts[3]),
                'left_len': int(parts[5]),
                'right_len': int(parts[6]),
                'left_seq': '',
                'right_seq': ''
            }
            
            # Read sequence files
            data_dir = os.path.join(output_dir, 'gapDataBase')
            prefix = f"{chrom}.{gap_info['id']}"
            for side in ['left', 'right']:
                file = os.path.join(data_dir, f"{prefix}.{side}.fasta")
                if os.path.exists(file):
                    with open(file) as fr:
                        fr.readline()  # Skip header line
                        seq = fr.readline().strip()
                        gap_info[f"{side}_seq"] = seq

            if chrom not in gaps_dict:
                gaps_dict[chrom] = []
            gaps_dict[chrom].append(gap_info)

    # Process each chromosome
    filtered_data = []
    print("\n=== Gap filtering details ===")
    for chrom, gaps in gaps_dict.items():
        print(f"Chromosome {chrom}: original gap count = {len(gaps)}")
        # Stage 1: Merge adjacent gaps
        merged_gaps = []
        for gap in sorted(gaps, key=lambda x: x['start']):
            if not merged_gaps:
                merged_gaps.append(gap)
            else:
                last = merged_gaps[-1]
                # Calculate interval distance (note coordinate conversion)
                distance = gap['start'] - last['end'] - 1
                if distance < SeedLength:
                    # Ensure both last and gap dictionaries have 'id' key
                    last_id = last.get('id', last.get('original_id', 'unknown'))
                    gap_id = gap.get('id', gap.get('original_id', 'unknown'))
                    print(f"  Merging gaps: {last_id}(right_len={last['right_len']}) and {gap_id}(left_len={gap['left_len']}, right_len={gap['right_len']}), distance={distance}bp")

                    # Use left_seq from first gap and right_seq from last gap
                    merged_gap = {
                        'start': last['start'],
                        'end': gap['end'],
                        'left_len': last['left_len'],
                        'right_len': gap['right_len'],
                        'left_seq': last['left_seq'],  # Keep left_seq from first gap
                        'right_seq': gap['right_seq']  # Use right_seq from last gap
                    }

                    # Copy other possible keys
                    for key in ['id', 'original_id']:
                        if key in last:
                            merged_gap[key] = last[key]

                    merged_gaps[-1] = merged_gap
                else:
                    merged_gaps.append(gap)

        print(f"  Gap count after merging = {len(merged_gaps)}")
        
        # Stage 2: Filtering logic
        total_gaps = len(merged_gaps)
        chrom_filtered_gaps = []  # Store filtered gaps for current chromosome
        
        for idx, gap in enumerate(merged_gaps, 1):
            # Determine if this is first or last gap after merging
            is_first_gap = (idx == 1)
            is_last_gap = (idx == total_gaps)

            # Check left length for first gap
            if is_first_gap and gap['left_len'] < SeedLength:
                gap_id = gap.get('id', gap.get('original_id', idx))
                print(f"  Skipping first gap (ID={gap_id}): left_len={gap['left_len']} < SeedLength={SeedLength}")
                continue

            # Check right length for last gap
            if is_last_gap and gap['right_len'] < SeedLength:
                gap_id = gap.get('id', gap.get('original_id', idx))
                print(f"  Skipping last gap (ID={gap_id}): right_len={gap['right_len']} < SeedLength={SeedLength}")
                continue

            # Save valid gap to current chromosome list
            chrom_filtered_gaps.append({
                'chrom': chrom,
                'original_id': idx,  # Keep original ID for debugging
                'start': gap['start'],
                'end': gap['end'],
                'gap_len': gap['end'] - gap['start'] + 1,
                'left_len': gap['left_len'],
                'right_len': gap['right_len'],
                'left_seq': gap['left_seq'],
                'right_seq': gap['right_seq']
            })
        
        # Reassign IDs for current chromosome gaps, starting from 1
        for new_id, gap in enumerate(chrom_filtered_gaps, 1):
            gap['id'] = new_id
            filtered_data.append(gap)

        print(f"  Gap count after filtering = {len(chrom_filtered_gaps)}")
        if len(chrom_filtered_gaps) > 0:
            print(f"  New ID range: {chrom}.1 - {chrom}.{len(chrom_filtered_gaps)}")
        print()

    print(f"Total: {len(filtered_data)} gaps retained after filtering\n")
    
    # Write final results
    log_path = os.path.join(output_dir, 'gap.log')
    with open(log_path, 'w') as fw:
        fw.write("Chromosome\tId\tStart\tEnd\tGapLen\tLeftLen\tRightLen\n")
        for item in filtered_data:
            # Write log entry
            fw.write(f"{item['chrom']}\t{item['id']}\t{item['start']}\t"
                     f"{item['end']}\t{item['gap_len']}\t"
                     f"{item['left_len']}\t{item['right_len']}\n")

            # Save sequence files
            for side in ['left', 'right']:
                seq = item[f"{side}_seq"]
                if seq:
                    filename = f"{item['chrom']}.{item['id']}.{side}.fasta"
                    with open(os.path.join(filtered_dir, filename), 'w') as fs:
                        fs.write(f">{filename[:-6]}\n{seq}\n")

def generate_job_script(reads_path, output_dir, optional_params=None, flag='all', mode='gapfiller'):
    """
    Generate job scripts for each gap

    Parameters:
    reads_path: HiFi file path
    output_dir: Output directory path
    optional_params: Optional parameters dictionary, including -ks, -kn, -kl, -j, --remove, --edge, --filterDepth, --MaximumExtensionLength etc.
    flag: Fill direction, default is 'all' (run both left and right)
    mode: Operation mode, default is 'gapfiller'
    """
    # Create DEGAP2.0_Output directory
    degap_output = os.path.join(output_dir, 'DEGAP2.0_Output')
    os.makedirs(degap_output, exist_ok=True)

    # Create jobScripts directory
    job_scripts_dir = os.path.join(output_dir, 'jobScripts')
    os.makedirs(job_scripts_dir, exist_ok=True)

    # Read filtered gap information
    gap_log = os.path.join(output_dir, 'gap.log')
    if not os.path.exists(gap_log):
        print(f"Error: Cannot find gap log file {gap_log}")
        return

    # Create main job script
    main_script = os.path.join(job_scripts_dir, 'all_gaps.sh')
    with open(main_script, 'w') as main_f:
        main_f.write("#!/bin/bash\n\n")

        # Read gap information and generate individual job scripts
        with open(gap_log) as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                chrom = parts[0]
                gap_id = parts[1]
                gap_name = f"{chrom}.{gap_id}"
                
                # Determine directions to run
                directions = []
                if flag == 'all':
                    directions = ['left', 'right']
                else:
                    directions = [flag]
                
                for direction in directions:
                    # Create gap output directory - differentiate by direction
                    gap_output = os.path.join(degap_output, f"{gap_name}.{direction}")
                    os.makedirs(gap_output, exist_ok=True)
                    
                    # Create symbolic links
                    idx_link = os.path.join(gap_output, 'reads.idx')
                    stat_link = os.path.join(gap_output, 'HiFi.reads.stat')
                    reads_part_link = os.path.join(gap_output, 'reads_part')
                    
                    # Create soft links
                    if not os.path.exists(idx_link):
                        os.symlink(os.path.abspath(os.path.join(output_dir, 'reads.idx')), idx_link)
                    if not os.path.exists(stat_link):
                        os.symlink(os.path.abspath(os.path.join(output_dir, 'HiFi.reads.stat')), stat_link)
                    if not os.path.exists(reads_part_link):
                        reads_part_dir = os.path.abspath(os.path.join(output_dir, 'reads_part'))
                        os.symlink(reads_part_dir, reads_part_link)
                    
                    # Build DEGAP.py command
                    left_seq = os.path.abspath(os.path.join(output_dir, 'gapData', f"{gap_name}.left.fasta"))
                    right_seq = os.path.abspath(os.path.join(output_dir, 'gapData', f"{gap_name}.right.fasta"))
                    
                    # Create individual gap job script - name by direction
                    gap_script = os.path.join(job_scripts_dir, f"run_{gap_name}.{direction}.sh")
                    with open(gap_script, 'w') as gap_f:
                        gap_f.write("#!/bin/bash\n\n")
                        
                        # Build basic command - ensure using absolute paths
                        # Get current script directory, assume DEGAP.py is in same directory as GetJobScript.py
                        script_dir = os.path.dirname(os.path.abspath(__file__))
                        degap_path = os.path.join(script_dir, 'DEGAP.py')
                        cmd = f"python {degap_path} --mode {mode} --reads {os.path.abspath(reads_path)} "
                        cmd += f"--seqleft {left_seq} --seqright {right_seq} --flag {direction} -o {os.path.abspath(gap_output)}"
                        
                        # Add optional parameters
                        if optional_params:
                            for param, value in optional_params.items():
                                if value is not None:
                                    cmd += f" {param} {value}"
                        
                        gap_f.write(cmd + "\n")
                    
                    # Set execution permissions
                    os.chmod(gap_script, 0o755)
                    
                    # Add individual gap script to main script
                    main_f.write(f"bash {os.path.abspath(gap_script)}\n")
    
    # Set main script execution permissions
    os.chmod(main_script, 0o755)
    
    print(f"Job scripts generated and saved in: {os.path.abspath(job_scripts_dir)}")
    print("Run the following command to execute all gap jobs:")
    print(f"bash {os.path.abspath(main_script)}")

def generate_batch_lsf_commands(output_dir, batch_count=10, cores=3, queue="normal", flag='all'):
    """
    Generate batched LSF batch submission scripts, tasks within each batch execute serially, batches execute in parallel

    Parameters:
    output_dir: Output directory path
    batch_count: Number of batches to divide into, default is 10
    cores: Number of cores per job
    queue: LSF queue name
    flag: Fill direction, default is all
    """
    # Input and output file paths
    job_scripts_dir = os.path.join(output_dir, 'jobScripts')
    run_all_gaps_path = os.path.join(job_scripts_dir, 'all_gaps.sh')
    batch_script_path = os.path.join(job_scripts_dir, f'lsf_all_jobs_batch_{batch_count}.sh')
    
    # Check if input file exists
    if not os.path.exists(run_all_gaps_path):
        print(f"Error: Cannot find job script file {run_all_gaps_path}")
        return [], None
    
    # Read run_all_gaps.sh file and get all tasks
    tasks = []
    gap_tasks = {}  # Tasks grouped by gap name
    
    with open(run_all_gaps_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and comment lines
            if not line or line.startswith('#'):
                continue
            
            # Extract script path and gap name
            if line.startswith('bash '):
                script_path = line.replace('bash ', '').strip()
                # Extract gap name from path (e.g. run_chr1A_TA299.1.sh -> chr1A_TA299.1)
                match = re.search(r'run_([^/]+)\.([^.]+)\.sh', script_path)
                if match:
                    gap_name = match.group(1)
                    direction = match.group(2)
                    
                    # If all mode, group tasks by gap name
                    if flag == 'all':
                        if gap_name not in gap_tasks:
                            gap_tasks[gap_name] = []
                        gap_tasks[gap_name].append((direction, f"{gap_name}.{direction}", script_path))
                    else:
                        tasks.append((f"{gap_name}.{direction}", script_path))
    
    # If all mode, prioritize left direction tasks
    if flag == 'all':
        for gap_name, gap_directions in gap_tasks.items():
            # Sort by direction, ensure left comes first
            sorted_directions = sorted(gap_directions, key=lambda x: 0 if x[0] == 'left' else 1)
            for _, task_name, script_path in sorted_directions:
                tasks.append((task_name, script_path))
    
    # If no tasks, exit
    if not tasks:
        print("Warning: No tasks found")
        return [], None
    
    # Calculate number of tasks per batch
    total_tasks = len(tasks)
    tasks_per_batch = (total_tasks + batch_count - 1) // batch_count  # Round up
    print(f"Total {total_tasks} tasks, divided into {batch_count} batches, approximately {tasks_per_batch} tasks per batch")
    
    # Create batch processing script
    all_job_names = []  # Store all job names for subsequent dependencies
    
    with open(batch_script_path, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Auto-generated batched LSF batch submission script\n")
        f.write(f"# Total {total_tasks} tasks, divided into {batch_count} batches, approximately {tasks_per_batch} tasks per batch\n\n")
        
        # If all mode, add check function
        if flag == 'all':
            f.write("# Function to check if gap is closed\n")
            f.write("check_gap_closed() {\n")
            f.write("    local gap_name=$1\n")
            f.write("    local direction=$2\n")
            f.write("    local result_dir=\"${DEGAP_OUTPUT}/${gap_name}.${direction}\"\n")
            f.write("    local process_log=\"${result_dir}/process.log\"\n")
            f.write("    \n")
            f.write("    if [ -f \"${process_log}\" ]; then\n")
            f.write("        if grep -q \"GAP can be closed\" \"${process_log}\"; then\n")
            f.write("            return 0  # Closed\n")
            f.write("        fi\n")
            f.write("    fi\n")
            f.write("    return 1  # Not closed\n")
            f.write("}\n\n")
            
            # Set output directory variable
            f.write("# Set DEGAP output directory\n")
            f.write(f"DEGAP_OUTPUT=\"{os.path.abspath(os.path.join(output_dir, 'DEGAP2.0_Output'))}\"\n\n")
        
        # Batch processing tasks
        for batch_idx in range(batch_count):
            batch_start = batch_idx * tasks_per_batch
            batch_end = min((batch_idx + 1) * tasks_per_batch, total_tasks)
            batch_tasks = tasks[batch_start:batch_end]
            
            if not batch_tasks:
                continue
                
            f.write(f"# Batch {batch_idx + 1}: Task {batch_start + 1}-{batch_end} (Total {len(batch_tasks)} tasks)\n")
            
            # First task in batch
            first_task = batch_tasks[0]
            job_scripts_dir_abs = os.path.abspath(job_scripts_dir)
            output_file = os.path.join(job_scripts_dir_abs, f"{first_task[0]}.out")
            error_file = os.path.join(job_scripts_dir_abs, f"{first_task[0]}.err")
            
            # Generate job name, ensure LSF naming convention
            job_name = f"batch{batch_idx+1}_1"
            all_job_names.append(job_name)
            
            # First task directly submitted
            f.write(f"bsub -J \"{job_name}\" -n {cores} -o {output_file} -e {error_file} -q {queue} \"bash {first_task[1]}\"\n")
            
            # Subsequent tasks in batch, each task depends on the previous task
            for i, (gap_name, script_path) in enumerate(batch_tasks[1:], 2):
                output_file = os.path.join(job_scripts_dir_abs, f"{gap_name}.out")
                error_file = os.path.join(job_scripts_dir_abs, f"{gap_name}.err")
                
                job_name = f"batch{batch_idx+1}_{i}"
                all_job_names.append(job_name)
                
                # If all mode, check if previous task is left direction of same gap
                if flag == 'all':
                    # Extract current task's gap name and direction
                    match = re.search(r'([^/]+)\.([^.]+)$', gap_name)
                    if match:
                        current_gap = match.group(1)
                        current_direction = match.group(2)
                        
                        # Extract previous task's gap name and direction
                        prev_gap_name = batch_tasks[i-2][0]  
                        prev_match = re.search(r'([^/]+)\.([^.]+)$', prev_gap_name)
                        
                        if prev_match:
                            prev_gap = prev_match.group(1)
                            prev_direction = prev_match.group(2)
                            
                            # If previous task is left direction of same gap and current is right direction, add check logic
                            if current_gap == prev_gap and prev_direction == 'left' and current_direction == 'right':
                                # Add conditional check: only execute right direction when left direction failed to close gap
                                f.write(f"# Check if {current_gap} left direction has closed gap\n")
                                f.write(f"if check_gap_closed \"{current_gap}\" \"left\"; then\n")
                                f.write(f"    echo \"Gap {current_gap} successfully closed by left direction, skipping right direction\"\n")
                                f.write(f"else\n")
                                f.write(f"    bsub -J \"{job_name}\" -w \"done(batch{batch_idx+1}_{i-1})\" -n {cores} -o {output_file} -e {error_file} -q {queue} \"bash {script_path}\"\n")
                                f.write(f"fi\n")
                                continue
                
                # Use -w option to specify dependency, execute only after previous task completes
                f.write(f"bsub -J \"{job_name}\" -w \"done(batch{batch_idx+1}_{i-1})\" -n {cores} -o {output_file} -e {error_file} -q {queue} \"bash {script_path}\"\n")
            
            f.write("\n")  # Add blank line between batches
    
    # Set execution permissions
    os.chmod(batch_script_path, 0o755)

    print(f"Batched LSF batch submission script generated: {batch_script_path}")
    print(f"Run the following command to submit all jobs:\nbash {batch_script_path}")
    print(f"Note: This will run {batch_count} batches simultaneously, tasks within each batch will execute serially")
    
    return all_job_names, batch_script_path

def split_reads(reads=args.reads,out=args.out):
    # Split reads by record count (100,000 records per file)
    import glob
    import subprocess
    
    # Use absolute paths
    reads_abs = os.path.abspath(reads)
    out_abs = os.path.abspath(out)
    reads_part_dir = os.path.join(out_abs, "reads_part")
    
    if not os.path.exists(reads_part_dir):
        os.makedirs(reads_part_dir)

    # Check if directory is empty, if not empty it means already split
    split_files = glob.glob(f"{reads_part_dir}/*.fa*")

    if split_files:
        print(f"Found {len(split_files)} existing split files in {reads_part_dir}, skipping split step")
    else:
        # Use Python's subprocess module to execute seqkit command
        print("No split files found, performing reads splitting...")
        cmd = ["seqkit", "split", reads, "-O", reads_part_dir, "--force", "--by-size", "100000", "--two-pass"]
        try:
            subprocess.run(cmd, check=True)
            print(f"Split completed, files saved in: {reads_part_dir}")
        except subprocess.CalledProcessError as e:
            print(f"Split failed: {e}")
            sys.exit(1)
    
    return reads_part_dir

def generate_lsf_commands(output_dir, cores=3, queue="normal", flag='all'):
    """
    Generate job submission script for LSF batch processing system

    Parameters:
    output_dir: Output directory path
    cores: Number of cores per job
    queue: LSF queue name
    flag: Fill direction, default is all
    """
    # Input and output file paths
    job_scripts_dir = os.path.join(output_dir, 'jobScripts')
    run_all_gaps_path = os.path.join(job_scripts_dir, 'all_gaps.sh')
    bsub_script_path = os.path.join(job_scripts_dir, 'lsf_all_jobs.sh')
    
    # Check if input file exists
    if not os.path.exists(run_all_gaps_path):
        print(f"Error: Cannot find job script file {run_all_gaps_path}")
        return
    
    # Read run_all_gaps.sh file
    with open(run_all_gaps_path, 'r') as f:
        lines = f.readlines()
    
    # Create bsub command script
    with open(bsub_script_path, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Auto-generated LSF batch submission script\n\n")
        
        # If all mode, add monitoring and control functionality
        if flag == 'all':
            # Add function to check if gap is closed
            f.write("# Function to check if gap is closed\n")
            f.write("check_gap_closed() {\n")
            f.write("    local gap_name=$1\n")
            f.write("    local direction=$2\n")
            f.write("    local result_dir=\"${DEGAP_OUTPUT}/${gap_name}.${direction}\"\n")
            f.write("    local process_log=\"${result_dir}/process.log\"\n")
            f.write("    \n")
            f.write("    if [ -f \"${process_log}\" ]; then\n")
            f.write("        if grep -q \"GAP can be closed\" \"${process_log}\"; then\n")
            f.write("            return 0  # Closed\n")
            f.write("        fi\n")
            f.write("    fi\n")
            f.write("    return 1  # Not closed\n")
            f.write("}\n\n")
            
            # Add monitoring and job stopping function
            f.write("# Monitor gap filling status and stop other direction job when necessary\n")
            f.write("monitor_and_stop() {\n")
            f.write("    local gap_name=$1\n")
            f.write("    local job_id_left=$2\n")
            f.write("    local job_id_right=$3\n")
            f.write("    \n")
            f.write("    # Wait for any job to complete\n")
            f.write("    while true; do\n")
            f.write("        # Check if left direction has completed and successfully closed gap\n")
            f.write("        if ! bjobs -a $job_id_left > /dev/null 2>&1 || [ \"$(bjobs -a $job_id_left | grep -c 'DONE')\" -gt 0 ]; then\n")
            f.write("            if check_gap_closed \"$gap_name\" \"left\"; then\n")
            f.write("                echo \"Gap $gap_name successfully closed by left direction, stopping right direction job $job_id_right\"\n")
            f.write("                bkill $job_id_right > /dev/null 2>&1\n")
            f.write("                # Create marker file indicating right direction was actively terminated\n")
            f.write("                touch \"${DEGAP_OUTPUT}/${gap_name}.right/TERMINATED_BY_LEFT_SUCCESS\"\n")
            f.write("                return 0\n")
            f.write("            fi\n")
            f.write("        fi\n")
            f.write("        \n")
            f.write("        # Check if right direction has completed and successfully closed gap\n")
            f.write("        if ! bjobs -a $job_id_right > /dev/null 2>&1 || [ \"$(bjobs -a $job_id_right | grep -c 'DONE')\" -gt 0 ]; then\n")
            f.write("            if check_gap_closed \"$gap_name\" \"right\"; then\n")
            f.write("                echo \"Gap $gap_name successfully closed by right direction, stopping left direction job $job_id_left\"\n")
            f.write("                bkill $job_id_left > /dev/null 2>&1\n")
            f.write("                # Create marker file indicating left direction was actively terminated\n")
            f.write("                touch \"${DEGAP_OUTPUT}/${gap_name}.left/TERMINATED_BY_RIGHT_SUCCESS\"\n")
            f.write("                return 0\n")
            f.write("            fi\n")
            f.write("        fi\n")
            f.write("        \n")
            f.write("        # Check if both jobs have ended\n")
            f.write("        if (! bjobs -a $job_id_left > /dev/null 2>&1 || [ \"$(bjobs -a $job_id_left | grep -c 'DONE\\|EXIT')\" -gt 0 ]) && \\\n")
            f.write("           (! bjobs -a $job_id_right > /dev/null 2>&1 || [ \"$(bjobs -a $job_id_right | grep -c 'DONE\\|EXIT')\" -gt 0 ]); then\n")
            f.write("            echo \"Both direction jobs for gap $gap_name have ended\"\n")
            f.write("            return 0\n")
            f.write("        fi\n")
            f.write("        \n")
            f.write("        # Sleep for a while before checking again\n")
            f.write("        sleep 60\n")
            f.write("    done\n")
            f.write("}\n\n")
            
            # Set output directory variable
            f.write("# Set DEGAP output directory\n")
            f.write(f"DEGAP_OUTPUT=\"{os.path.abspath(os.path.join(output_dir, 'DEGAP2.0_Output'))}\"\n\n")
        
        # If all mode, organize tasks by gap name
        if flag == 'all':
            # Parse all tasks and group by gap name
            gap_tasks = {}
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                if line.startswith('bash '):
                    script_path = line.replace('bash ', '').strip()
                    match = re.search(r'run_([^/]+)\.([^.]+)\.sh', script_path)
                    if match:
                        gap_name = match.group(1)
                        direction = match.group(2)
                        
                        if gap_name not in gap_tasks:
                            gap_tasks[gap_name] = {}
                        gap_tasks[gap_name][direction] = script_path
            
            # Submit left and right direction jobs for each gap and set up monitoring
            for gap_name, directions in gap_tasks.items():
                if 'left' in directions and 'right' in directions:
                    left_script = directions['left']
                    right_script = directions['right']
                    
                    # Build left direction bsub command
                    job_scripts_dir_abs = os.path.abspath(job_scripts_dir)
                    left_output = os.path.join(job_scripts_dir_abs, f"{gap_name}.left.out")
                    left_error = os.path.join(job_scripts_dir_abs, f"{gap_name}.left.err")
                    
                    # Submit left direction job and get job ID
                    f.write(f"# Submit left direction job for {gap_name}\n")
                    f.write(f"left_job_id=$(bsub -J {gap_name}.left -n {cores} -o {left_output} -e {left_error} -q {queue} \"bash {left_script}\" | grep -o '[0-9]\\+')\n")
                    
                    # Build right direction bsub command
                    right_output = os.path.join(job_scripts_dir_abs, f"{gap_name}.right.out")
                    right_error = os.path.join(job_scripts_dir_abs, f"{gap_name}.right.err")
                    
                    # Submit right direction job and get job ID
                    f.write(f"# Submit right direction job for {gap_name}\n")
                    f.write(f"right_job_id=$(bsub -J {gap_name}.right -n {cores} -o {right_output} -e {right_error} -q {queue} \"bash {right_script}\" | grep -o '[0-9]\\+')\n")
                    
                    # Start monitoring process
                    f.write(f"# Start monitoring process\n")
                    f.write(f"monitor_and_stop \"{gap_name}\" \"$left_job_id\" \"$right_job_id\" &\n\n")
                elif 'left' in directions:
                    # Only left direction
                    script_path = directions['left']
                    full_name = f"{gap_name}.left"
                    
                    job_scripts_dir_abs = os.path.abspath(job_scripts_dir)
                    output_file = os.path.join(job_scripts_dir_abs, f"{full_name}.out")
                    error_file = os.path.join(job_scripts_dir_abs, f"{full_name}.err")
                    
                    bsub_cmd = f"bsub -J {full_name} -n {cores} -o {output_file} -e {error_file} -q {queue} \"bash {script_path}\"\n"
                    f.write(bsub_cmd)
                elif 'right' in directions:
                    # Only right direction
                    script_path = directions['right']
                    full_name = f"{gap_name}.right"
                    
                    job_scripts_dir_abs = os.path.abspath(job_scripts_dir)
                    output_file = os.path.join(job_scripts_dir_abs, f"{full_name}.out")
                    error_file = os.path.join(job_scripts_dir_abs, f"{full_name}.err")
                    
                    bsub_cmd = f"bsub -J {full_name} -n {cores} -o {output_file} -e {error_file} -q {queue} \"bash {script_path}\"\n"
                    f.write(bsub_cmd)
        else:
            # Non-all mode, directly submit all jobs
            for line in lines:
                line = line.strip()
                # Skip empty lines and comment lines
                if not line or line.startswith('#'):
                    continue
                
                # Extract script path and gap name
                if line.startswith('bash '):
                    script_path = line.replace('bash ', '').strip()
                    # Extract gap name from path (e.g. run_chr1A_TA299.1.left.sh -> chr1A_TA299.1.left)
                    match = re.search(r'run_([^/]+)\.([^.]+)\.sh', script_path)
                    if match:
                        gap_name = match.group(1)
                        direction = match.group(2)
                        full_name = f"{gap_name}.{direction}"
                        
                        # Build bsub command
                        # Use absolute paths to specify output and error files
                        job_scripts_dir_abs = os.path.abspath(job_scripts_dir)
                        output_file = os.path.join(job_scripts_dir_abs, f"{full_name}.out")
                        error_file = os.path.join(job_scripts_dir_abs, f"{full_name}.err")
                        
                        bsub_cmd = f"bsub -J {full_name} -n {cores} -o {output_file} -e {error_file} -q {queue} \"bash {script_path}\"\n"
                        f.write(bsub_cmd)
    
    # Set execution permissions
    os.chmod(bsub_script_path, 0o755)

    print(f"LSF batch submission script generated: {bsub_script_path}")
    print(f"Run the following command to submit all jobs:\nbash {bsub_script_path}")

def extract_filled_gap_info(output_dir, genome_file, flag='all'):
    """
    Extract filled gap information and generate gapDegap.log file

    Parameters:
    output_dir: Output directory path
    genome_file: Original genome file path
    flag: Fill direction, default is all (check both left and right directions)

    Returns:
    success: Whether successful
    new_genome_sequences: New genome sequences
    filled_gaps_info: List of filled gap information
    unfilled_gaps_info: List of unfilled gap information
    """
    # Convert to absolute paths
    output_dir = os.path.abspath(output_dir)
    genome_file = os.path.abspath(genome_file)

    print(f"Using output directory absolute path: {output_dir}")
    print(f"Using genome file absolute path: {genome_file}")
    
    # Read configuration parameters
    config_file = os.path.join(output_dir, 'degap_config.json')
    seed_len = None  # Initialize as None, will try to get from HiFi.reads.stat file later

    # First try to read SeedLength from HiFi.reads.stat file
    stats_path = os.path.join(output_dir, "HiFi.reads.stat")
    if os.path.exists(stats_path):
        try:
            with open(stats_path) as f:
                print(f"Reading HiFi.reads.stat file: {stats_path}")
                for line in f:
                    print(f"  Reading line: {line.strip()}")
                    if line.startswith("SeedLength"):
                        seed_len = int(line.strip().split('\t')[1])
                        print(f"Read seed length from HiFi.reads.stat file: {seed_len}")
                        break
        except Exception as e:
            print(f"Error reading HiFi.reads.stat file: {e}")
    else:
        print(f"HiFi.reads.stat file does not exist: {stats_path}")

    # If reading from HiFi.reads.stat file failed, try reading from config file
    if seed_len is None and os.path.exists(config_file):
        try:
            with open(config_file, 'r') as f:
                config = json.load(f)
                seed_len = config.get('SeedLength')
                if seed_len:
                    print(f"Read seed length from config file: {seed_len}")
                else:
                    print(f"No SeedLength field in config file")
        except Exception as e:
            print(f"Error reading config file: {e}")
    elif seed_len is None:
        print(f"Config file does not exist: {config_file}")

    # If both methods failed, use default value
    if seed_len is None:
        seed_len = 1000  # Default seed length
        print(f"Unable to read seed length from files, using default value: {seed_len}")

    print(f"Final seed length used: {seed_len}")
    
    # Read gap.log, only process filtered and merged gaps
    gap_log_path = os.path.join(output_dir, 'gap.log')
    if not os.path.exists(gap_log_path):
        print(f"Error: Cannot find gap log file {gap_log_path}")
        return False, {}, [], []
    
    # Read chromosome information (without loading sequences)
    print("\n=== Reading genome chromosome information ===")
    chromosome_ids = []
    chromosome_lengths = {}
    
    try:
        # First only read chromosome ID and length information
        with open(genome_file, 'r') as f:
            current_chrom = None
            for line in f:
                if line.startswith('>'):
                    current_chrom = line.strip()[1:].split()[0]
                    chromosome_ids.append(current_chrom)
                    chromosome_lengths[current_chrom] = 0
                elif current_chrom is not None:
                    chromosome_lengths[current_chrom] += len(line.strip())
        
        print(f"Found {len(chromosome_ids)} chromosomes")
        for chrom, length in chromosome_lengths.items():
            print(f"  - {chrom}: {length} bp")
    except Exception as e:
        print(f"Error reading genome file: {e}")
        return False, {}, [], []
    
    # Read gapbase.log to get original gap count (for statistics only)
    original_gap_count = 0
    gapbase_log_path = os.path.join(output_dir, 'gapbase.log')
    if os.path.exists(gapbase_log_path):
        with open(gapbase_log_path, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5 and int(parts[4]) >= 20:  # Filter out gaps that are too small
                    original_gap_count += 1

    # Read filtered and merged gap information
    gaps_info = []
    with open(gap_log_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            
            chrom, gap_id, start, end = parts[0], parts[1], int(parts[2]), int(parts[3])
            gap_len, left_len, right_len = int(parts[4]), int(parts[5]), int(parts[6])
            
            gaps_info.append({
                'chrom': chrom,
                'id': gap_id,
                'start': start,
                'end': end,
                'gap_len': gap_len,
                'left_len': left_len,
                'right_len': right_len,
                'filled': False,
                'filled_seq': '',
                'seed_len': seed_len  # Add seed length information
            })
    
    # Save original gap count for subsequent statistics
    original_gaps_count = original_gap_count

    # Check if results directory exists
    degap_output_dir = os.path.join(output_dir, 'DEGAP2.0_Output')
    if not os.path.exists(degap_output_dir):
        print(f"Warning: Cannot find DEGAP2.0_Output directory: {degap_output_dir}")
        print("Trying to find other possible result directories...")
        # Try to find other possible result directories
        possible_dirs = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d)) and 'output' in d.lower()]
        if possible_dirs:
            print(f"Found possible result directories: {possible_dirs}")
            degap_output_dir = os.path.join(output_dir, possible_dirs[0])
            print(f"Using directory: {degap_output_dir}")
        else:
            print("No possible result directories found, will continue trying to use original path")
    
    # Process filling results for each gap
    filled_gaps = []
    unfilled_gaps = []

    # Information for detailed logging
    filled_gaps_info = []
    unfilled_gaps_info = []

    # Determine directions to check
    directions = []
    if flag == 'all':
        directions = ['left', 'right']
    else:
        directions = [flag]
    
    # Output debug information
    print(f"\nSearching for result files, base path: {degap_output_dir}")
    
    for gap in gaps_info:
        chrom = gap['chrom']
        gap_id = gap['id']
        gap_name = f"{chrom}.{gap_id}"
        
        # Check results for each direction
        best_fill = None
        fail_reasons = []

        for direction in directions:
            # Only use standard paths
            result_dir = os.path.join(degap_output_dir, f"{gap_name}.{direction}")

            print(f"Checking directory: {result_dir}")

            if not os.path.exists(result_dir):
                reason = f"Cannot find {direction} direction result directory"
                fail_reasons.append(reason)
                print(f"Warning: {reason}: {result_dir}")
                continue
            
            # Check if job was actively terminated
            terminated_by_other_success = False
            if direction == 'left':
                terminated_file = os.path.join(result_dir, "TERMINATED_BY_RIGHT_SUCCESS")
                if os.path.exists(terminated_file):
                    print(f"Gap {gap_name} left direction job was actively terminated after right direction successful filling")
                    terminated_by_other_success = True
                    # Skip checking this direction because the other direction has succeeded
                    continue
            elif direction == 'right':
                terminated_file = os.path.join(result_dir, "TERMINATED_BY_LEFT_SUCCESS")
                if os.path.exists(terminated_file):
                    print(f"Gap {gap_name} right direction job was actively terminated after left direction successful filling")
                    terminated_by_other_success = True
                    # Skip checking this direction because the other direction has succeeded
                    continue
            
            # Result file and log paths
            final_fa = os.path.join(result_dir, f"{gap_name}.final.fa")
            process_log = os.path.join(result_dir, "process.log")
            
            print(f"Checking file: {final_fa}")
            print(f"Checking file: {process_log}")
            
            # List files in directory for debugging
            if os.path.exists(result_dir):
                print(f"Files in directory {result_dir}:")
                files_in_dir = os.listdir(result_dir)
                for file in files_in_dir:
                    print(f"  - {file}")
                
                # If standard filename doesn't exist, try using correct filename format
                if not os.path.exists(final_fa):
                    # Try using filename format with direction suffix
                    alternative_final_fa = os.path.join(result_dir, f"{gap_name}.{direction}.final.fa")
                    if os.path.exists(alternative_final_fa):
                        print(f"Found alternative file: {alternative_final_fa}")
                        final_fa = alternative_final_fa
                    else:
                        # Look for files containing final and fa
                        final_files = [f for f in files_in_dir if 'final' in f.lower() and ('.fa' in f.lower() or '.fasta' in f.lower())]
                        if final_files:
                            final_fa = os.path.join(result_dir, final_files[0])
                            print(f"Using found final file: {final_fa}")
            
            if not os.path.exists(final_fa) or not os.path.exists(process_log):
                reason = f"Cannot find {direction} direction result file or log"
                fail_reasons.append(reason)
                print(f"Warning: {reason}")
                continue

            # Extract gap length information from log
            gap_length_from_log = None
            is_closed = False
            
            with open(process_log, 'r') as f:
                content = f.read()
                
                # Check if successfully filled
                if "GAP can be closed!" in content:
                    is_closed = True
                    
                    # Extract filled gap length
                    for line in content.split('\n'):
                        if "GAP Length:" in line:
                            try:
                                gap_length_from_log = int(line.strip().split(':')[1].strip())
                                break
                            except (ValueError, IndexError):
                                pass
            
            if not is_closed or gap_length_from_log is None:
                reason = f"{direction} direction not successfully filled"
                fail_reasons.append(reason)
                print(f"Gap {gap_name} {reason}")
                continue
            
            # Read filled sequence
            filled_seq = None
            try:
                seq_records = list(SeqIO.parse(final_fa, "fasta"))
                if not seq_records:
                    print(f"Warning: No sequence records found in {final_fa}")
                    reason = f"No sequence records in {direction} direction final.fa file"
                    fail_reasons.append(reason)
                    continue
                
                # Read first sequence
                record = seq_records[0]
                seq_str = str(record.seq)
                print(f"Gap {gap_name} {direction} direction, total sequence length: {len(seq_str)} bp, gap length to extract: {gap_length_from_log} bp, seed_len: {seed_len} bp")
                
                if len(seq_str) <= seed_len:
                    print(f"Warning: Sequence length {len(seq_str)} is less than or equal to seed length {seed_len}, cannot extract filled sequence")
                    reason = f"{direction} direction sequence length insufficient to extract filled portion"
                    fail_reasons.append(reason)
                    continue
                
                # Extract filled sequence
                try:
                    # Handle negative gap length (indicating overlap rather than gap)
                    if gap_length_from_log < 0:
                        print(f"Warning: Gap length is negative ({gap_length_from_log}), indicating overlap rather than gap")
                        # For overlap cases, we use an empty string as the filled sequence (since no filling is needed)
                        filled_seq = ""
                        print(f"Overlap case: Using empty string as filled sequence")
                    else:
                        if direction == 'left':
                            # Count from left to right, starting at SeedLen+1
                            filled_seq = seq_str[seed_len:seed_len+gap_length_from_log]
                            print(f"Extracting from left to right: offset={seed_len}, length={gap_length_from_log}, extracted sequence length={len(filled_seq)}")
                        else:  # direction == 'right'
                            # Count from right to left, starting at SeedLen+1
                            filled_seq = seq_str[-seed_len-gap_length_from_log:-seed_len]
                            print(f"Extracting from right to left: offset={-seed_len-gap_length_from_log} to {-seed_len}, extracted sequence length={len(filled_seq)}")
                except Exception as e:
                    print(f"Error extracting filled sequence: {str(e)}")
                    reason = f"{direction} direction error calculating sequence indexes: {str(e)}"
                    fail_reasons.append(reason)
                    continue
            except Exception as e:
                reason = f"Error filling sequence in {direction} direction: {str(e)}"
                fail_reasons.append(reason)
                print(f"{reason}")
                continue
            
            # Check extraction results
            if filled_seq is None:
                reason = f"{direction} direction unable to extract filled sequence"
                fail_reasons.append(reason)
                print(f"Gap {gap_name} {reason}")
                continue
            
            # Check extracted sequence length
            if gap_length_from_log < 0:
                # For negative gap length (overlap case), we have already set filled_seq to an empty string
                if len(filled_seq) == 0:
                    print(f"Overlap case: Using empty string as filled sequence, as expected")
                else:
                    reason = f"{direction} direction error handling overlap case (expected empty string, actual {len(filled_seq)}bp)"
                    fail_reasons.append(reason)
                    print(f"Gap {gap_name} {reason}")
                    continue
            elif len(filled_seq) != gap_length_from_log:
                print(f"Warning: Extracted sequence length ({len(filled_seq)}) does not match expected gap length ({gap_length_from_log})")
                # Try to tolerate small range of length differences to accommodate gap length and actual filled length minor discrepancies
                if abs(len(filled_seq) - gap_length_from_log) <= 5:  # Allow 5bp error
                    print(f"Error within acceptable range, continuing to use extracted sequence")
                else:
                    reason = f"{direction} direction extracted sequence length does not match expected ({gap_length_from_log}bp, actual {len(filled_seq)}bp)"
                    fail_reasons.append(reason)
                    print(f"Gap {gap_name} {reason}")
                    continue
            
            # Save filling results
            if gap_length_from_log < 0:
                print(f"Gap {gap_name} {direction} direction detected sequence overlap, no filling needed, overlap length: {abs(gap_length_from_log)} bp")
                # Save current direction filling result - overlap case
                current_fill = {
                    'direction': direction,
                    'filled_seq': filled_seq,  # Empty string
                    'gap_length': gap_length_from_log  # Keep negative value, indicates this is overlap
                }
            else:
                print(f"Gap {gap_name} {direction} direction successfully filled, length: {len(filled_seq)} bp")
                # Save current direction filling result - normal filling
                current_fill = {
                    'direction': direction,
                    'filled_seq': filled_seq,
                    'gap_length': len(filled_seq)  # Use actual extracted sequence length
                }
            
            # If no best result yet or current is left direction (prioritize left direction), update best result
            if best_fill is None or (direction == 'left' and flag == 'all'):
                best_fill = current_fill
        
        # If best filling result found, save to filled_gaps
        if best_fill:
            print(f"Gap {gap_name} using {best_fill['direction']} direction filling result")
            gap['filled'] = True
            gap['filled_seq'] = best_fill['filled_seq']
            gap['fill_direction'] = best_fill['direction']
            gap['fill_length'] = best_fill['gap_length']
            filled_gaps.append(gap)

            # Add to detailed log information
            filled_gaps_info.append({
                'chrom': gap['chrom'],
                'id': gap['id'],
                'start': gap['start'],
                'end': gap['end'],
                'gap_len': gap['gap_len'],
                'fill_direction': best_fill['direction'],
                'fill_length': best_fill['gap_length'],
                'filled_seq': best_fill['filled_seq']
            })
        else:
            gap['fail_reason'] = "; ".join(fail_reasons) if fail_reasons else "No valid filling result found"
            unfilled_gaps.append(gap)

            # Add to detailed log information
            unfilled_gaps_info.append({
                'chrom': gap['chrom'],
                'id': gap['id'],
                'start': gap['start'],
                'end': gap['end'],
                'gap_len': gap['gap_len'],
                'fail_reason': gap['fail_reason']
            })
    
    # Output statistics
    print(f"\nTotal: Among original {original_gaps_count} gaps, after merging gaps with insufficient SeedLen spacing, ignoring gaps at chromosome ends where seed sequences cannot be obtained, {len(gaps_info)} gaps can be merged using DEGAP, {len(filled_gaps)} successfully filled")

    # Process genome sequences in segments, chromosome by chromosome
    print("\n=== Starting segmented genome sequence processing ===")
    new_genome_sequences = {}

    # Group gaps by chromosome for processing
    gaps_by_chrom = {}
    for gap in filled_gaps:
        chrom = gap['chrom']
        if chrom not in gaps_by_chrom:
            gaps_by_chrom[chrom] = []
        gaps_by_chrom[chrom].append(gap)
    
    # Process chromosomes one by one
    for chrom_idx, chrom in enumerate(chromosome_ids, 1):
        print(f"\nProcessing chromosome {chrom} ({chrom_idx}/{len(chromosome_ids)})...")

        # Read single chromosome sequence
        chrom_seq = ""
        found_chrom = False
        
        with open(genome_file, 'r') as f:
            current_reading_chrom = None
            for line in f:
                if line.startswith('>'):
                    current_reading_chrom = line.strip()[1:].split()[0]
                    if current_reading_chrom == chrom:
                        found_chrom = True
                elif current_reading_chrom == chrom:
                    chrom_seq += line.strip()
                elif found_chrom:
                    # Already finished reading target chromosome, can break loop
                    break
        
        if not found_chrom:
            print(f"Warning: Chromosome {chrom} not found in genome file")
            continue

        print(f"Chromosome {chrom} sequence length: {len(chrom_seq)} bp")
        
        # Convert to list for modification
        new_seq = list(chrom_seq)

        # Find all successfully filled gaps for this chromosome
        chrom_filled_gaps = gaps_by_chrom.get(chrom, [])
        print(f"Chromosome {chrom} has {len(chrom_filled_gaps)} gaps successfully filled")

        # Replace from back to front by position (avoid position changes)
        for gap in sorted(chrom_filled_gaps, key=lambda x: x['start'], reverse=True):
            start_idx = gap['start'] - 1  # Convert to 0-based
            end_idx = gap['end'] - 1
            filled_seq = gap['filled_seq']

            print(f"  Filling gap {gap['id']} (position: {start_idx+1}-{end_idx+1}), length: {len(filled_seq)} bp")
            
            # Replace sequence
            if len(filled_seq) == 0 and gap['fill_length'] < 0:
                # Handle overlap case - first check if it's direct connection
                gap_name = f"{gap['chrom']}.{gap['id']}"
                is_direct_connection = False
                direct_connection_dir = None
                final_fa_path = None
                
                # Check result directories for both directions
                possible_dirs = []
                for dir_name in ['left', 'right']:
                    result_dir = os.path.join(degap_output_dir, f"{gap_name}.{dir_name}")
                    if os.path.exists(result_dir):
                        possible_dirs.append((dir_name, result_dir))
                
                # First determine if it's direct connection through markers in process.log
                for dir_name, result_dir in possible_dirs:
                    process_log = os.path.join(result_dir, "process.log")
                    if os.path.exists(process_log):
                        try:
                            with open(process_log, 'r') as f:
                                content = f.read()
                                if "GAP can be closed!" in content and "Direct overlap" in content:
                                    is_direct_connection = True
                                    direct_connection_dir = result_dir
                                    print(f"  Detected direct connection marker in {dir_name} direction process.log")
                                    break
                        except Exception as e:
                            print(f"  Error reading process.log file: {str(e)}")

                # If confirmed as direct connection, find corresponding final.fa file
                if is_direct_connection and direct_connection_dir:
                    # Try to find final.fa file
                    # First check standard naming
                    standard_final_fa = os.path.join(direct_connection_dir, f"{gap_name}.final.fa")
                    if os.path.exists(standard_final_fa) and os.path.getsize(standard_final_fa) > 0:
                        final_fa_path = standard_final_fa
                    
                    # If standard naming doesn't exist, check naming with direction suffix
                    if not final_fa_path:
                        dir_suffix = os.path.basename(direct_connection_dir).split('.')[-1]  # Get direction suffix
                        suffixed_final_fa = os.path.join(direct_connection_dir, f"{gap_name}.{dir_suffix}.final.fa")
                        if os.path.exists(suffixed_final_fa) and os.path.getsize(suffixed_final_fa) > 0:
                            final_fa_path = suffixed_final_fa
                    
                    # If none of the above naming exists, try to find any file containing final and fa
                    if not final_fa_path:
                        for file in os.listdir(direct_connection_dir):
                            if 'final' in file.lower() and ('.fa' in file.lower() or '.fasta' in file.lower()):
                                final_fa_path = os.path.join(direct_connection_dir, file)
                                break
                
                if is_direct_connection and final_fa_path and os.path.exists(final_fa_path) and os.path.getsize(final_fa_path) > 0:
                    print(f"  Detected direct connection file: {final_fa_path}")
                    # Read direct connection sequence
                    final_seq = ""
                    try:
                        records = list(SeqIO.parse(final_fa_path, "fasta"))
                        if records:
                            final_seq = str(records[0].seq)
                    except Exception as e:
                        print(f"  Error reading final.fa file: {str(e)}")
                        # Try simple reading
                        try:
                            with open(final_fa_path, 'r') as f:
                                for i, line in enumerate(f):
                                    if i > 0:  # Skip first line (sequence identifier)
                                        final_seq += line.strip()
                        except Exception as e2:
                            print(f"  Simple reading of final.fa file also failed: {str(e2)}")
                    
                    if final_seq:
                        print(f"  Read direct connection sequence, length: {len(final_seq)} bp")
                        
                        try:
                            # Get seed length from configuration or default value
                            seed_len = gap.get('seed_len', 1000)  # Use seed length from configuration or default value
                            print(f"  Using seed length: {seed_len} bp (read from configuration)")
                            
                            # Get seed_len bp sequence on left side of gap
                            left_seed_start = max(0, start_idx - seed_len)
                            left_seed = ''.join(new_seq[left_seed_start:start_idx])
                            
                            # Get seed_len bp sequence on right side of gap
                            right_seed_start = end_idx + 1
                            right_seed_end = min(len(new_seq), right_seed_start + seed_len)
                            right_seed = ''.join(new_seq[right_seed_start:right_seed_end])
                            
                            print(f"  Left seed sequence length: {len(left_seed)} bp")
                            print(f"  Right seed sequence length: {len(right_seed)} bp")

                            # Use shortest left and right seed sequence length for matching
                            min_seed_match = min(100, min(len(left_seed), len(right_seed)) // 2)  # Minimum length required for matching

                            # Check if left seed sequence is at the beginning of final_seq
                            match_len = min(len(left_seed), min_seed_match)
                            if match_len > 0 and final_seq.startswith(left_seed[:match_len]):
                                print(f"  final.fa sequence matches left seed sequence")
                                
                                # Check if right seed sequence is in final_seq
                                match_len_right = min(len(right_seed), min_seed_match)
                                if match_len_right > 0:
                                    # Check if right seed is in final_seq
                                    right_in_final = final_seq.find(right_seed[:match_len_right])
                                    if right_in_final > 0:
                                        print(f"  Found right seed sequence in final.fa sequence, position: {right_in_final}")

                                        # Use final_seq to replace region from left seed sequence start to right seed sequence end
                                        replace_start = left_seed_start
                                        replace_end = right_seed_end

                                        print(f"  Using direct connection sequence to replace region: {replace_start+1}-{replace_end}, range length: {replace_end-replace_start} bp")

                                        # Replace sequence
                                        new_seq[replace_start:replace_end] = list(final_seq)
                                        print(f"  Successfully used final.fa file ({len(final_seq)} bp) to replace overlap region")
                                        continue
                                    else:
                                        print(f"  Right seed sequence not found in final.fa sequence")
                            else:
                                print(f"  final.fa sequence does not match left seed sequence")
                        
                        except Exception as e:
                            print(f"  Error processing direct connection file: {str(e)}")

                # If direct connection file not found or processing failed, fall back to original handling
                print(f"  Handling overlap case: deleting {gap['gap_len']}bp N region")
                # For overlap case, we delete N region (i.e., don't add any sequence)
                new_seq[start_idx:end_idx+1] = []
            else:
                # Normal case - replace N region with filled sequence
                new_seq[start_idx:end_idx+1] = filled_seq
            
            # Update position information after filling (for detailed logging)
            for filled_gap in filled_gaps_info:
                if filled_gap['chrom'] == gap['chrom'] and filled_gap['id'] == gap['id']:
                    filled_gap['new_start'] = gap['start']
                    filled_gap['new_end'] = gap['start'] + len(filled_seq) - 1
                    break
        
        # Convert back to string and save
        new_genome_sequences[chrom] = ''.join(new_seq)
        print(f"Chromosome {chrom} processing complete, new sequence length: {len(new_genome_sequences[chrom])} bp")

        # Actively release memory
        del new_seq
        del chrom_seq
        import gc
        gc.collect()
    
    # Re-detect gaps in new sequences and generate gapDegap.log
    print("\n=== Generating post-filling gap log ===")
    degap_log_path = os.path.join(output_dir, 'gapDegap.log')

    with open(degap_log_path, 'w') as log_file:
        log_file.write("Chromosome\tId\tStart\tEnd\tGapLen\tLeftLen\tRightLen\n")

        # Process chromosomes one by one
        for chrom in chromosome_ids:
            if chrom not in new_genome_sequences:
                print(f"Warning: Skipping missing chromosome {chrom}")
                continue

            seq = new_genome_sequences[chrom]
            print(f"Detecting gaps in chromosome {chrom}...")

            # Detect all N regions (gaps), no merging
            raw_gaps = []
            in_gap = False
            gap_start = 0
            
            # Scan sequence in segments to reduce memory pressure
            segment_size = 10000000  # 10Mb段
            for i in range(0, len(seq), segment_size):
                segment = seq[i:i+segment_size]
                
                for j, c in enumerate(segment):
                    abs_pos = i + j
                    if c.upper() == 'N':
                        if not in_gap:
                            in_gap = True
                            gap_start = abs_pos
                    else:
                        if in_gap:
                            in_gap = False
                            raw_gaps.append((gap_start, abs_pos - 1))
            
            # Handle possible gap at the end
            if in_gap:
                raw_gaps.append((gap_start, len(seq) - 1))
            
            print(f"Detected {len(raw_gaps)} gaps in chromosome {chrom}")

            # Directly process each raw gap, no merging
            for gap_id, (start, end) in enumerate(raw_gaps, 1):
                # Convert chromosome coordinates (1-based)
                chr_start = start + 1
                chr_end = end + 1
                gap_len = end - start + 1
                
                # Get left sequence length - use same SeedLength limit as original gap extraction
                left_len = 0
                current_pos = start - 1
                while current_pos >= 0 and seq[current_pos].upper() != 'N':
                    left_len += 1
                    current_pos -= 1
                    # Use seed sequence length limit, consistent with extract_gap_info function
                    if left_len >= seed_len:
                        break
                
                # Get right sequence length - use same SeedLength limit as original gap extraction
                right_len = 0
                current_pos = end + 1
                while current_pos < len(seq) and seq[current_pos].upper() != 'N':
                    right_len += 1
                    current_pos += 1
                    # Use seed sequence length limit, consistent with extract_gap_info function
                    if right_len >= seed_len:
                        break
                
                # Write to log - keep even gaps on chromosome sides
                log_file.write(
                    f"{chrom}\t{gap_id}\t{chr_start}\t{chr_end}\t{gap_len}\t"
                    f"{left_len}\t{right_len}\n"
                )
    
    print(f"Post-filling gap log saved to: {degap_log_path}")

    # Return processing results
    return True, new_genome_sequences, filled_gaps_info, unfilled_gaps_info

# Add a new function to create post-processing job
def create_post_processing_job(output_dir, genome_file, flag, all_job_names):
    """
    Create a post-processing job that depends on all gap filling tasks

    Parameters:
    output_dir: Output directory path
    genome_file: Original genome file path
    flag: Fill direction
    all_job_names: List of all gap filling job names
    """
    job_scripts_dir = os.path.join(output_dir, 'jobScripts')
    post_script_path = os.path.join(job_scripts_dir, 'post_processing.sh')
    post_job_name = "degap_post"
    
    # Create post-processing script
    with open(post_script_path, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Auto-generated post-processing script\n\n")

        # Get current script directory, assume DEGAP.py is in same directory as GetJobScript.py
        f.write("# Get current script directory\n")
        f.write("SCRIPT_DIR=$(dirname \"$(readlink -f \"$0\")\")\n")
        f.write("cd \"$SCRIPT_DIR/../\"\n\n")

        # Add post-processing command
        f.write("# Post-processing command\n")
        f.write(f"python {os.path.abspath(__file__)} --post_process --reads {os.path.abspath(args.reads)} --genome {os.path.abspath(args.genome)} --out {os.path.abspath(output_dir)} --flag {flag}\n")
    
    # Set execution permissions
    os.chmod(post_script_path, 0o755)

    # Create dependency submission command
    dependency_str = " && ".join([f"done({job})" for job in all_job_names])

    # Build bsub command
    job_scripts_dir_abs = os.path.abspath(job_scripts_dir)
    output_file = os.path.join(job_scripts_dir_abs, f"{post_job_name}.out")
    error_file = os.path.join(job_scripts_dir_abs, f"{post_job_name}.err")
    
    bsub_cmd = f"bsub -J \"{post_job_name}\" -w \"{dependency_str}\" -n 1 -o {output_file} -e {error_file} -q {args.queue} \"bash {post_script_path}\""
    
    # Create submission script
    submit_script_path = os.path.join(job_scripts_dir, 'submit_post_job.sh')
    with open(submit_script_path, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write(bsub_cmd + "\n")

    # Set execution permissions
    os.chmod(submit_script_path, 0o755)

    print(f"Post-processing job script generated: {post_script_path}")
    print(f"Post-processing job submission script generated: {submit_script_path}")
    
    return submit_script_path

def check_all_gaps_finished(output_dir, flag='all'):
    """
    Check if all gap jobs are completed

    Parameters:
    output_dir: Output directory path
    flag: Fill direction, default is all (check both left and right directions)

    Returns:
    finished: Boolean value indicating whether all jobs are completed
    total_gaps: Total number of gaps
    finished_gaps: Number of completed gaps
    """
    # Read gap.log to get all gaps that need processing
    gap_log_path = os.path.join(output_dir, 'gap.log')
    if not os.path.exists(gap_log_path):
        print(f"Error: Cannot find gap log file {gap_log_path}")
        return False, 0, 0

    # Determine directions to check
    directions = []
    if flag == 'all':
        directions = ['left', 'right']
    else:
        directions = [flag]
    
    # Read all gap information
    all_gaps = []
    with open(gap_log_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                gap_id = parts[1]
                all_gaps.append((chrom, gap_id))
    
    # If all mode, count each gap only once, not each direction
    if flag == 'all':
        total_gaps = len(all_gaps)
    else:
        total_gaps = len(all_gaps) * len(directions)
        
    finished_gaps = 0
    
    # Check if each direction of each gap is completed
    degap_output_dir = os.path.join(output_dir, 'DEGAP2.0_Output')
    if not os.path.exists(degap_output_dir):
        print(f"Warning: Cannot find DEGAP2.0_Output directory: {degap_output_dir}")
        return False, total_gaps, 0
    
    # If in 'all' mode, a gap is considered complete if any direction is finished
    if flag == 'all':
        for chrom, gap_id in all_gaps:
            gap_name = f"{chrom}.{gap_id}"
            gap_finished = False

            # Check if any direction successfully closed the gap
            for direction in directions:
                result_dir = os.path.join(degap_output_dir, f"{gap_name}.{direction}")
                if not os.path.exists(result_dir):
                    continue

                # Check if process.log indicates gap has been closed
                process_log = os.path.join(result_dir, "process.log")
                if os.path.exists(process_log):
                    with open(process_log, 'r') as f:
                        content = f.read()
                        if "GAP can be closed" in content:
                            gap_finished = True
                            break

                # Check if termination marker file exists
                terminated_file = os.path.join(result_dir, "TERMINATED_BY_LEFT_SUCCESS") if direction == "right" else os.path.join(result_dir, "TERMINATED_BY_RIGHT_SUCCESS")
                if os.path.exists(terminated_file):
                    # If termination marker file exists, it means the other direction has successfully filled the gap
                    gap_finished = True
                    break

            # If not successfully closed or terminated, check if all directions have final.fa files
            if not gap_finished:
                direction_finished_count = 0
                for direction in directions:
                    result_dir = os.path.join(degap_output_dir, f"{gap_name}.{direction}")
                    if not os.path.exists(result_dir):
                        continue

                    # Check final.fa file
                    final_fa = os.path.join(result_dir, f"{gap_name}.final.fa")
                    # If standard filename doesn't exist, try using direction suffix filename format
                    if not os.path.exists(final_fa):
                        alternative_final_fa = os.path.join(result_dir, f"{gap_name}.{direction}.final.fa")
                        if os.path.exists(alternative_final_fa):
                            final_fa = alternative_final_fa
                        else:
                            # Try to find any file containing 'final' and 'fa'
                            for file in os.listdir(result_dir):
                                if 'final' in file.lower() and ('.fa' in file.lower() or '.fasta' in file.lower()):
                                    final_fa = os.path.join(result_dir, file)
                                    break

                    if os.path.exists(final_fa) and os.path.getsize(final_fa) > 0:
                        direction_finished_count += 1

                # If all directions are complete, consider this gap complete
                if direction_finished_count == len(directions):
                    gap_finished = True

            if gap_finished:
                finished_gaps += 1
    else:
        # Non-all mode, check according to original logic
        for chrom, gap_id in all_gaps:
            for direction in directions:
                gap_name = f"{chrom}.{gap_id}"
                result_dir = os.path.join(degap_output_dir, f"{gap_name}.{direction}")

                # Check if result directory exists
                if not os.path.exists(result_dir):
                    continue

                # Check if termination marker file exists
                terminated_file = os.path.join(result_dir, "TERMINATED_BY_LEFT_SUCCESS") if direction == "right" else os.path.join(result_dir, "TERMINATED_BY_RIGHT_SUCCESS")
                if os.path.exists(terminated_file):
                    # If termination marker file exists, it means the other direction has successfully filled the gap, consider this direction as completed
                    finished_gaps += 1
                    continue

                # Check if final.fa file exists and is not empty
                final_fa = os.path.join(result_dir, f"{gap_name}.final.fa")

                # If standard filename doesn't exist, try using direction suffix filename format
                if not os.path.exists(final_fa):
                    alternative_final_fa = os.path.join(result_dir, f"{gap_name}.{direction}.final.fa")
                    if os.path.exists(alternative_final_fa):
                        final_fa = alternative_final_fa
                    else:
                        # Try to find any file containing 'final' and 'fa'
                        for file in os.listdir(result_dir):
                            if 'final' in file.lower() and ('.fa' in file.lower() or '.fasta' in file.lower()):
                                final_fa = os.path.join(result_dir, file)
                                break

                if os.path.exists(final_fa) and os.path.getsize(final_fa) > 0:
                    finished_gaps += 1

    # Determine if all jobs are completed
    finished = (finished_gaps == total_gaps)

    return finished, total_gaps, finished_gaps

def main(args, stats):
    """Main process function that executes the entire gap filling workflow"""
    # Check if in post-processing mode
    if args.post_process:
        # Execute only post-processing steps
        print("\n=== Starting post-processing steps ===")
        genome_abs = os.path.abspath(args.genome)
        output_dir_abs = os.path.abspath(args.out)
        flag = args.flag if hasattr(args, 'flag') else 'all'

        print("\n=== Starting extraction of filled gap information ===")
        success, new_genome_sequences, filled_gaps_info, unfilled_gaps_info = extract_filled_gap_info(output_dir_abs, genome_abs, flag=flag)
        if not success:
            print("Failed to extract filled gap information")
            _exit_program(1)

        print("Filled gap information has been saved to:", os.path.join(output_dir_abs, 'gapDegap.log'))

        # Generate detailed process.log file
        generate_detailed_process_log(output_dir_abs, filled_gaps_info, unfilled_gaps_info)

        print("\n=== Starting to plot post-filling gap distribution ===")
        # Use the same plotting function, but with post-filling gap log
        if plot(genome_abs, os.path.join(output_dir_abs, 'gapDegap.log'), output_dir_abs, is_after_filling=True):
            print("Post-filling gap distribution plot has been saved to:", os.path.join(output_dir_abs, 'gap_plots'))
        else:
            print("Failed to generate post-filling gap distribution plot")

        # Save filled genome sequences
        print("\n=== Saving filled genome sequences ===")
        filled_genome_path = os.path.join(output_dir_abs, 'genome.filled.fasta')
        with open(filled_genome_path, 'w') as f:
            for chrom, seq in new_genome_sequences.items():
                f.write(f">{chrom}\n{seq}\n")

        print(f"Filled genome sequences have been saved to: {filled_genome_path}")
        print("\n=== Post-processing completed ===")

        # Add completion marker file
        completion_flag = os.path.join(output_dir_abs, 'DEGAP_PROCESS_COMPLETED')
        with open(completion_flag, 'w') as f:
            f.write(f"DEGAP processing completion time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        print(f"Completion marker file created: {completion_flag}")

        # If called from monitoring for post-processing, force exit; otherwise return normally
        if getattr(args, 'force_exit_after_processing', False):
            print("\n=== All steps completed, script will exit ===")
            _exit_program(0)
        return

    # The following is the pre-processing stage code
    genome_abs = os.path.abspath(args.genome)
    output_dir_abs = os.path.abspath(args.out)

    # Save configuration parameters to JSON file to ensure subsequent processing can correctly read SeedLength
    config_file = os.path.join(output_dir_abs, 'degap_config.json')
    try:
        config = {}
        if os.path.exists(config_file):
            with open(config_file, 'r') as f:
                config = json.load(f)

        # Update configuration parameters
        config['SeedLength'] = stats.get('SeedLength')

        # Write configuration file
        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)
        print(f"Configuration parameters saved to: {config_file}")
    except Exception as e:
        print(f"Error saving configuration parameters: {e}")

    # Stage 1: Extract and filter Gap information
    print("\n=== Starting Gap information extraction ===")
    extract_gap_info(genome_abs, output_dir_abs, stats['SeedLength'])
    print("Gap information has been saved to:", os.path.join(output_dir_abs, 'gapbase.log'))

    print("\n=== Starting Gap information filtering ===")
    filter_gaps(output_dir_abs, stats['SeedLength'])
    print("Filtered Gap information has been saved to:", os.path.join(output_dir_abs, 'gap.log'))
    # Breakpoint line


    print("\n=== Starting to plot gap distribution ===")
    if plot(genome_abs, os.path.join(output_dir_abs, 'gapbase.log'), output_dir_abs, is_after_filling=False):
        print("Pre-filling gap distribution plot has been saved to:", os.path.join(output_dir_abs, 'gap_plots'))
    else:
        print("Failed to generate pre-filling gap distribution plot")

    # Stage 2: Prepare filling jobs
    print(f"\n=== Starting to split HiFi reads to {os.path.join(output_dir_abs, 'reads_part')} ===")
    reads_abs = os.path.abspath(args.reads)
    reads_part_dir = split_reads(reads=reads_abs, out=output_dir_abs)

    print("\n=== Obtaining command submission script ===")
    # Build optional parameters dictionary
    optional_params = {}
    for param in ['--kmer_size', '--kmer_num', '--kmer_len', '-j', '--remove', '--edge', '--filterDepth', '--MaximumExtensionLength']:
        param_name = param.lstrip('-')
        if hasattr(args, param_name) and getattr(args, param_name) is not None:
            optional_params[param] = getattr(args, param_name)
    # Use flag and mode from command line arguments
    flag = args.flag if hasattr(args, 'flag') else 'all'
    mode = args.mode if hasattr(args, 'mode') and args.mode else 'gapfiller'
    generate_job_script(reads_abs, output_dir_abs, optional_params, flag, mode)
    
    # Stage 3: Submit and wait for jobs
    print("\n=== Generating LSF batch submission script ===")
    generate_lsf_commands(output_dir_abs)

    # Decide which submission method to use based on batch parameter
    all_job_names = []
    batch_script_path = None
    jobs_submitted = False

    # Get flag parameter
    flag = args.flag if hasattr(args, 'flag') else 'all'

    if args.batch is not None:
        # Use batch mode
        print(f"\n=== Generating batch LSF submission script (split into {args.batch} batches) ===")
        all_job_names, batch_script_path = generate_batch_lsf_commands(output_dir_abs, batch_count=args.batch, flag=flag, queue=args.queue)

        if all_job_names and len(all_job_names) > 0 and os.path.exists(batch_script_path):
            print("\n=== Executing batch LSF submission script ===")
            print(f"Executing batch LSF submission script: {batch_script_path}")
            try:
                subprocess.run(['bash', batch_script_path], check=True)
                print("LSF job submission successful")
                jobs_submitted = True
            except subprocess.CalledProcessError as e:
                print(f"LSF job submission failed: {e}")
        elif not os.path.exists(batch_script_path):
            print(f"Error: Cannot find batch script {batch_script_path}")
        else:
            print("Warning: No jobs generated, skipping batch job submission")
    else:
        # Use regular mode (all jobs in parallel, but with monitoring mechanism)
        print("\n=== Executing regular LSF batch submission script (all jobs in parallel with monitoring) ===")
        # Generate LSF command script, passing flag and queue parameters
        generate_lsf_commands(output_dir_abs, flag=flag, queue=args.queue)

        lsf_script_path = os.path.join(output_dir_abs, 'jobScripts', 'lsf_all_jobs.sh')
        if os.path.exists(lsf_script_path):
            print(f"Executing LSF batch submission script: {lsf_script_path}")
            try:
                subprocess.run(['bash', lsf_script_path], check=True)
                print("LSF job submission successful")
                jobs_submitted = True
            except subprocess.CalledProcessError as e:
                print(f"LSF job submission failed: {e}")
        else:
            print(f"Error: Cannot find batch script {lsf_script_path}")

    # If jobs were successfully submitted, start monitoring
    if jobs_submitted:
        # No longer use LSF dependency feature, but periodically check job completion
        print("\n=== Starting job completion monitoring ===")
        print("Checking completion status of all gap jobs every 30 minutes")
        print("Press Ctrl+C to interrupt monitoring at any time, post-processing can be run manually later")

        try:
            while True:
                # Check if all gap jobs are completed
                finished, total, completed = check_all_gaps_finished(output_dir_abs, flag)
                print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Completed: {completed}/{total} ({completed/total*100:.1f}%)")

                if finished:
                    print("\nAll gap jobs completed, starting post-processing steps")
                    # Directly call post-processing function
                    post_process_args = argparse.Namespace()
                    post_process_args.post_process = True
                    post_process_args.reads = reads_abs
                    post_process_args.genome = genome_abs
                    post_process_args.out = output_dir_abs
                    post_process_args.flag = flag
                    post_process_args.force_exit_after_processing = True  # Add flag to indicate forced exit after post-processing

                    # Call post-processing
                    main(post_process_args, stats)
                    # If post-processing function doesn't exit, force exit here
                    print("\n=== All steps completed, script will exit ===")
                    _exit_program(0)  # Use forced exit function

                # Wait 30 minutes
                print(f"Waiting 30 minutes before next check...")
                time.sleep(30 * 60)  # 30 minutes

        except KeyboardInterrupt:
            print("\nMonitoring interrupted. You can manually run post-processing later:")
            print(f"python {os.path.abspath(__file__)} --post_process --reads {reads_abs} --genome {genome_abs} --out {output_dir_abs} --flag {flag}")
    else:
        print("Warning: No jobs successfully submitted, skipping monitoring and post-processing")

    print("\n=== Initialization stage completed ===")
    print("You can manually run post-processing later:")
    print(f"python {os.path.abspath(__file__)} --post_process --reads {reads_abs} --genome {genome_abs} --out {output_dir_abs} --flag {flag}")

def generate_detailed_process_log(output_dir, filled_gaps_info, unfilled_gaps_info):
    """
    Generate detailed process.log file to record gap filling details

    Parameters:
    output_dir: Output directory path
    filled_gaps_info: List of filled gap information
    unfilled_gaps_info: List of unfilled gap information
    """
    log_path = os.path.join(output_dir, 'detailed_process.log')
    print(f"\n=== Generating detailed process.log file: {log_path} ===")

    with open(log_path, 'w') as f:
        # Write title
        f.write("# DEGAP2.0 Gap Filling Detailed Report\n")
        f.write("# Generated at: " + time.strftime("%Y-%m-%d %H:%M:%S") + "\n\n")

        # Write statistics
        total_gaps = len(filled_gaps_info) + len(unfilled_gaps_info)
        filled_count = len(filled_gaps_info)
        unfilled_count = len(unfilled_gaps_info)

        f.write(f"## Statistics\n")
        f.write(f"- Total Gap Count: {total_gaps}\n")
        f.write(f"- Filled Gap Count: {filled_count}\n")
        f.write(f"- Unfilled Gap Count: {unfilled_count}\n")
        f.write(f"- Fill Success Rate: {filled_count/total_gaps*100:.2f}%\n\n")

        # Write detailed information of filled gaps
        f.write("## Filled Gap Details\n")
        f.write("| No. | Chromosome | Gap ID | Original Position | Original Length | Fill Direction | Fill Length | New Position |\n")
        f.write("|-----|------------|--------|-------------------|-----------------|----------------|-------------|---------------|\n")

        for i, gap in enumerate(filled_gaps_info, 1):
            chrom = gap['chrom']
            gap_id = gap['id']
            start = gap['start']
            end = gap['end']
            gap_len = gap['gap_len']
            direction = gap['fill_direction']
            fill_length = gap['fill_length']
            new_start = gap['new_start'] if 'new_start' in gap else start
            new_end = gap['new_end'] if 'new_end' in gap else end

            f.write(f"| {i} | {chrom} | {gap_id} | {start}-{end} | {gap_len} | {direction} | {fill_length} | {new_start}-{new_end} |\n")

        f.write("\n")

        # Write unfilled gap information
        f.write("## Unfilled Gap Information\n")
        f.write("| No. | Chromosome | Gap ID | Position | Length | Failure Reason |\n")
        f.write("|-----|------------|--------|----------|--------|-----------------|\n")

        for i, gap in enumerate(unfilled_gaps_info, 1):
            chrom = gap['chrom']
            gap_id = gap['id']
            start = gap['start']
            end = gap['end']
            gap_len = gap['gap_len']
            reason = gap.get('fail_reason', 'No valid filling result found')

            f.write(f"| {i} | {chrom} | {gap_id} | {start}-{end} | {gap_len} | {reason} |\n")

        f.write("\n")

        # Write conclusion
        f.write("## Conclusion\n")
        if filled_count > 0:
            f.write(f"Successfully filled {filled_count} gaps, accounting for {filled_count/total_gaps*100:.2f}% of the total.\n")
        else:
            f.write("Failed to fill any gaps.\n")

        f.write("\n")
        f.write("# End of Report\n")

    print(f"Detailed process.log file generated: {log_path}")

# Execute main process
if __name__ == "__main__":
    try:
        # If in post-processing mode, skip index creation
        if args.post_process:
            # Read existing statistics
            stats_path = os.path.join(args.out, 'HiFi.reads.stat')
            stats = {}
            if os.path.exists(stats_path):
                with open(stats_path) as f:
                    for line in f:
                        key, value = line.strip().split('\t')
                        try:
                            stats[key] = float(value)
                        except ValueError:
                            stats[key] = value
            main(args, stats)
        else:
            # Normal workflow
            stats = process_hifi(args.reads, args.out)
            print("\nStatistics results:")
            for k, v in stats.items():
                print(f"{k}: {v}")
            main(args, stats)
    except KeyboardInterrupt:
        print("\nProgram interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nProgram execution error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

# Script usage examples:
# Basic command:
# python GetJobScript.py --reads path/to/reads --genome path/to/genome --mode gapfiller --flag left -o path/to/outputdir
# You can change default parameter values through -kn|--kmer_num  -ks|--kmer_size  -kl|--kmer_len -j --remove --edge --filterDepth --MaximumExtensionLength
# Use --batch parameter to specify the number of batches for task execution, e.g.: --batch 5 means split into 5 batches; without this parameter all jobs are submitted in parallel
# Generated scripts are located in outputdir/jobScripts/ directory:
# - all_gaps.sh: Execution script for all tasks
# - lsf_all_jobs.sh: LSF submission script for all tasks (parallel execution of all tasks)
# - lsf_all_jobs_batch_X.sh: LSF submission script for batch execution (X is the number of batches, tasks within each batch execute serially, batches execute in parallel)