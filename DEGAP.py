import re
import os
import sys
import argparse  # Replace getopt
import pysam
from pysam import AlignmentFile
import math
import Bio
from Bio import SeqIO
import GapFiller
from GapFiller import GapFiller
import CtgLinker
from CtgLinker import CtgLinker
import selectRawReads
from selectRawReads import selectRawReads
import subprocess
import glob
import telfiller
from telfiller import TelFiller

def usage():
	print ("--reads HiFi_reads.fasta")
	print ("-o | --out ./path/")
	print ("-t | --thread thread number")
	print ("--remove 1 | 2 | 3 default:2 1:only keep final result; 2: keep every round basic result ; 3 : keep all files")
	print ("--edge Edge Controller set max Edge length(missequening)")
	print ("--filterDepth num default:None. You can filtered HiFi reads by mapped depth on contig set. if num==0.3 means: mapped Hifi reads on depth>=0.3*avgdepth and depth<=(2-0.2)*avgdepth will be filtered and will not be used in whole project")
	print ("--MaximumExtensionLength num default:None. Stop Extension when reach the num")
	print ("--kmer_size | -ks num default:41. k-mer size for filtering reads")
	print ("--kmer_num | -kn num default:10. number of k-mers to use for filtering reads")
	print ("--kmer_length | -kl num default:0.1. proportion of mean read length for k-mer extraction (0.1 = 10% of mean length)")
	print ("-j num default:100. number of parallel jobs for processing reads")
	print ("--mode gapfiller | ctglinker | telfiller")
	print ("--resume num Resume from specified round (e.g.: --resume 118 continues from round118)")
	print ("--resume_auto Automatically resume from last interrupted round")
	print ("\n\ngapfiller\n")
	print ("\t--seqleft sequence before GAP")
	print ("\t--seqright sequence after GAP")
	print ("\t--flag left | right choose seqleft or seqright as seed to fill the gap(default left)")
	print ("\n\ntelfiller\n")
	print ("\t--seqleft start sequence")
	print ("\t--seqright end sequence(s)")
	print ("\t--flag left | right choose extension direction (default left)")
	print ("\n\nctglinker\n")
	print ("\t--ctgseq contig set")
	print ("-h | --help help")

def getoptions():
	parser = argparse.ArgumentParser(description='DEGAP: Dynamic Elongation of a Genome Assembly Path',
	                                add_help=False,  # Disable auto-generated help
	                                formatter_class=argparse.RawTextHelpFormatter)

	# Basic parameters
	parser.add_argument('--mode', type=str, help='Operation mode: gapfiller or ctglinker or telfiller')
	parser.add_argument('-o', '--out', type=str, help='Output directory path')
	parser.add_argument('-t', '--thread', type=str, help='Number of threads')
	parser.add_argument('-j', type=int, default=100, help='Number of parallel jobs for processing reads')
	parser.add_argument('-h', '--help', help='Show help message and exit', action='store_true')
	parser.add_argument('--reads', type=str, help='HiFi reads file path')
	parser.add_argument('--remove', type=int, default=2, 
	                   help='1: only keep final result; 2: keep every round basic result; 3: keep all files')
	parser.add_argument('--edge', type=int, default=500, help='Edge Controller set max Edge length')
	parser.add_argument('--filterDepth', type=float, help='Filter HiFi reads by mapped depth')
	parser.add_argument('--MaximumExtensionLength', type=int, help='Stop Extension when reach this length')
	
	# K-mer related parameters
	parser.add_argument('--kmer_size', '-ks', type=int, default=41, help='k-mer size for filtering reads')
	parser.add_argument('--kmer_num', '-kn', type=int, default=10, help='number of k-mers to use for filtering reads')
	parser.add_argument('--kmer_length', '-kl', type=float, default=0.1,
	                   help='proportion of mean read length for k-mer extraction (0.1 = 10%% of mean length)')

	# Resume functionality
	parser.add_argument('--resume', type=int, help='Resume from specified round')
	parser.add_argument('--resume_auto', action='store_true', help='Automatically resume from last interrupted round')

	# Parameters shared by gapfiller and telomerefiller modes
	parser.add_argument('--seqleft', type=str, help='Sequence before GAP or start sequence for telomerefiller')
	parser.add_argument('--seqright', type=str, help='Sequence after GAP or end sequence for telomerefiller')
	parser.add_argument('--flag', type=str, default='left', help='Choose seqleft or seqright as seed to fill the gap or extension direction')
	
	# ctglinker mode specific parameters
	parser.add_argument('--ctgseq', type=str, help='Contig set file path')

	# Show help if no arguments provided
	if len(sys.argv) == 1:
		print("NO PARAMETER!!!")
		print("Use '-h | --help' for same information")
		sys.exit()
	
	args = parser.parse_args()
	
	# Handle help information
	if args.help:
		usage()
		sys.exit()

	# Check required parameters
	if not args.mode:
		print("Mode must be specified (gapfiller or ctglinker or telfiller)")
		sys.exit()
	
	if not args.reads:
		print("Reads file must be specified")
		sys.exit()
	
	if not args.out:
		print("Output directory must be specified")
		sys.exit()
	
	# Create kparameters list containing kmer and parallel related parameters
	kparameters = [args.kmer_size, args.kmer_num, args.kmer_length, args.j]

	# Resume functionality parameters
	resume_params = (args.resume, args.resume_auto)

	# Build return values (maintain same return format as original function)
	if args.mode == "gapfiller":
		# Check files required for gapfiller mode
		if args.seqleft:
			if os.path.exists(args.seqleft) and os.path.getsize(args.seqleft) != 0:
				seqleft = args.seqleft
			else:
				print("seqleft file doesn't exist or is empty!")
				sys.exit()
		else:
			print("seqleft parameter is required for gapfiller mode")
			sys.exit()
			
		if args.seqright:
			if os.path.exists(args.seqright) and os.path.getsize(args.seqright) != 0:
				seqright = args.seqright
			else:
				print("seqright file doesn't exist or is empty!")
				sys.exit()
		else:
			print("seqright parameter is required for gapfiller mode")
			sys.exit()
			
		return [args.mode, args.remove, args.thread or '20', args.reads, args.out, 
		        seqleft, seqright, args.flag, args.edge, args.filterDepth, 
		        args.MaximumExtensionLength], kparameters, resume_params
	
	elif args.mode == "ctglinker":
		# Check files required for ctglinker mode
		if args.ctgseq:
			if os.path.exists(args.ctgseq) and os.path.getsize(args.ctgseq) != 0:
				seqfile = args.ctgseq
			else:
				print("contigs file doesn't exist or is empty")
				sys.exit()
		else:
			print("ctgseq parameter is required for ctglinker mode")
			sys.exit()
			
		return [args.mode, args.remove, args.thread or '20', args.reads, args.out, 
		        seqfile, args.edge, args.filterDepth, 
		        args.MaximumExtensionLength], kparameters, resume_params
	
	elif args.mode == "telfiller":
		# Check files required for telfiller mode
		if args.seqleft:
			if os.path.exists(args.seqleft) and os.path.getsize(args.seqleft) != 0:
				seqleft = args.seqleft
			else:
				print("seqleft file doesn't exist or is empty!")
				sys.exit()
		else:
			print("seqleft parameter is required for telfiller mode")
			sys.exit()
			
		if args.seqright:
			if os.path.exists(args.seqright) and os.path.getsize(args.seqright) != 0:
				seqright = args.seqright
			else:
				print("seqright file doesn't exist or is empty!")
				sys.exit()
		else:
			print("seqright parameter is required for telfiller mode")
			sys.exit()
			
		return [args.mode, args.remove, args.thread or '20', args.reads, args.out, 
		        seqleft, seqright, args.flag, args.edge, args.filterDepth, 
		        args.MaximumExtensionLength], kparameters, resume_params
	
	else:
		print("You should use gapfiller, ctglinker, or telfiller!")
		sys.exit()

parameter, kparameters, resume_params = getoptions()

#parameter: mode,remove,thread,reads,out,seqleft,seqright,flag,edge,filterDepth,MaximumExtensionLength
reads=parameter[3]
out=parameter[4]

#make outpit file
if out[-1]=="/":
	out=out[:-1]
	parameter[4]=out
if os.path.exists(out)!=True:
	os.makedirs(out)

#make reads index
print ("BUILD RAWREADS DICT")
idx_path = out+"/reads.idx"
stats_path = out+"/HiFi.reads.stat"

# Check if index file exists and is valid
rebuild_index = True
if os.path.exists(idx_path) and os.path.getsize(idx_path) > 0:
    try:
        # Try to load index file
        print("Detected existing index file, attempting to load...")
        readsdict = SeqIO.index_db(idx_path)
        # Verify if index matches current reads file
        if readsdict:
            print("Index file loaded successfully, skipping index rebuild")
            rebuild_index = False
    except Exception as e:
        print(f"Failed to load index file: {e}")
        print("Will rebuild index")

if rebuild_index:
    print("Building reads index file...")
    readsdict = SeqIO.index_db(idx_path, reads, "fasta")
    print("Index construction completed")

print ("BUILD DICT SUCCEED")
parameter.append(readsdict)

print("Splitting reads file")
# Split reads by record count (100,000 records per file)
if os.path.exists(out+"/reads_part")!=True:
	os.makedirs(out+"/reads_part")

reads_part_dir = out+"/reads_part"
# Check if directory is empty, if not empty it means already split
split_files = glob.glob(f"{reads_part_dir}/*.fa*")

if split_files:
    print(f"Found {len(split_files)} existing split files in {reads_part_dir}, skipping split step")
else:
    # Use Python subprocess module to execute seqkit command
    print("No split files found, performing reads splitting...")
    cmd = ["seqkit", "split", reads, "-O", reads_part_dir, "--force", "--by-size", "100000", "--two-pass"]
    try:
        subprocess.run(cmd, check=True)
        print(f"Splitting completed")
    except subprocess.CalledProcessError as e:
        print(f"Splitting failed: {e}")
        sys.exit(1)

# Calculate read length statistics
pwd1 = stats_path
if os.path.exists(pwd1) and os.path.getsize(pwd1) > 0:
    # Read existing statistics file
    try:
        print("Detected existing statistics file, attempting to load...")
        file1 = open(pwd1, 'r')
        mean_length = None
        reads_count = None
        total_length = None
        lenmax = None
        seedlen = None

        for i in file1:
            i1 = i.rstrip().split("\t")
            if i1[0] == "MaxLength":
                lenmax = int(i1[1])
            elif i1[0] == "SeedLength":
                seedlen = int(float(i1[1]))
            elif i1[0] == "Number":
                reads_count = int(i1[1])
            elif i1[0] == "TolalLenth" or i1[0] == "TotalLength":  # Compatible with both spellings
                total_length = int(i1[1])
        file1.close()

        # Verify if all necessary information was obtained
        if lenmax is not None and seedlen is not None and reads_count is not None and total_length is not None:
            mean_length = total_length / reads_count
            print(f"Statistics file loaded successfully: average read length {mean_length:.2f} bp, maximum read length {lenmax} bp")
        else:
            print("Statistics file format incorrect or missing key information, will recalculate")
            mean_length = None  # Trigger recalculation
    except Exception as e:
        print(f"Failed to read statistics file: {e}")
        mean_length = None  # Trigger recalculation
else:
    mean_length = None  # Trigger recalculation

# If no valid statistics information, recalculate
if mean_length is None:
    print("Creating HiFi read length statistics file...")
    file1 = open(pwd1, 'w')

    # Use created index for statistics
    n = len(readsdict)  # Directly get number of entries in index
    total_length = 0
    lenmax = 0

    print(f"Starting precise statistics calculation for all {n} reads using index...")
    count = 0
    for read_id in readsdict:
        seq_len = len(readsdict[read_id].seq)
        if seq_len > lenmax:
            lenmax = seq_len
        total_length += seq_len
        count += 1
        if count % 100000 == 0:
            print(f"Processed {count}/{n} reads ({count/n*100:.1f}%)...")

    mean_length = total_length / n if n > 0 else 0
    print(f"Statistics complete: average read length {mean_length:.2f} bp, maximum read length {lenmax} bp")
    
    l = 'Number\t' + str(n) + "\nTolalLenth\t" + str(total_length) + "\nMaxLength\t" + str(lenmax) + "\n"
    l += "MeanLength\t" + str(int(mean_length)) + "\n"
    file1.writelines(l)
    a = 10**(int(math.log(lenmax, 10)))
    b = lenmax / a + 1
    seedlen = a * b
    l = 'SeedLength\t' + str(int(seedlen)) + "\n"
    file1.writelines(l)
    file1.close()

# Calculate actual kmer_length based on average read length
print(f"Average read length: {mean_length:.2f} bp")
actual_kmer_length = int(kparameters[2] * mean_length)
print(f"Using kmer_length ratio: {kparameters[2]:.2f}, actual kmer_length: {actual_kmer_length} bp")
kparameters[2] = actual_kmer_length  # Update kmer_length in kparameters to actual value

print (lenmax)
print (seedlen)
if parameter[-3]!=None:
	selectedReads=selectRawReads(parameter,seedlen)

	parameter=parameter[:-1]

	print ("BUILD RAWREADS DICT")
	readsdict2=SeqIO.index_db(out+"/selectedReads.idx",selectedReads.readFile,"fasta")
	print ("BUILD DICT SUCCEED")
	parameter.append(readsdict2)
	parameter[3]=selectedReads.readFile

# Ensure maximum read length and seed length are added to parameter list
parameter.append(lenmax)  # Add maxReadsLen
parameter.append(seedlen) # Add seedLen

print("Calculating average and maximum read lengths...")
print(f"Average read length: {mean_length:.2f}, maximum read length: {lenmax}")
print(f"Parameter list length: {len(parameter)}")  # Print parameter length for debugging
print(parameter)

# Unpack resume parameters
resume_round, resume_auto = resume_params

# Set up resume functionality
if resume_auto:
    # Automatically search for checkpoint file
    checkpoint_file = out+"/checkpoint.info"
    if os.path.exists(checkpoint_file):
        try:
            with open(checkpoint_file, 'r') as f:
                checkpoint_data = f.read().strip()
                if checkpoint_data and checkpoint_data.startswith("round:"):
                    resume_round = int(checkpoint_data.split("round:")[1])
                    print(f"Auto resume mode: Found checkpoint information, will continue from round{resume_round}")
        except Exception as e:
            print(f"Failed to read checkpoint information: {e}")
            resume_round = None

# Execute main functionality
if parameter[0] == "gapfiller":
    if resume_round is not None:
        # If specified which round to continue from, set resume_round attribute
        print(f"Preparing to continue from round{resume_round}...")
        gapfiller = GapFiller(parameter, kparameters)
        gapfiller.resume_round = resume_round
    else:
        # Normal startup
        GapFiller(parameter, kparameters)
elif parameter[0] == "ctglinker":
    if resume_round is not None:
        print(f"Warning: ctglinker mode does not support resume functionality yet, will ignore --resume parameter")
    # Normal startup for ctglinker
    CtgLinker(parameter, kparameters)
elif parameter[0] == "telfiller":
    if resume_round is not None:
        # If specified which round to continue from, set resume_round attribute
        print(f"Preparing to continue from round{resume_round}...")
        telfiller = TelFiller(parameter, kparameters)
        telfiller.resume_round = resume_round
    else:
        # Normal startup for telfiller
        TelFiller(parameter, kparameters)
print ('welldone')


