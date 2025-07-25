import re
import Bio
import os
import sys
import re
import getopt
import pysam
from pysam import AlignmentFile
import Bio
from Bio import SeqIO
from collections import defaultdict

class FindExtensionReads(object):
	def __init__(self,roundInput,lastRoundUsedReads,usedReads,kmer_size=41,kmer_num=10,kmer_length=1000,j=70,out=None):
		# Import required modules
		import os
		import traceback
		
		self.roundInput=roundInput
		self.lastRoundUsedReads=lastRoundUsedReads
		self.usedReads=usedReads
		self.note=''
		self.kmer_size=kmer_size
		self.kmer_num=kmer_num
		self.kmer_length=kmer_length
		self.j=j
		self.out=out
		self.extensionReads=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".fa"
		self.log=self.roundInput.elongation.roundDir+"/extensionReads.log"
		
		# Ensure output directory exists
		if not os.path.exists(self.roundInput.elongation.roundDir):
			try:
				os.makedirs(self.roundInput.elongation.roundDir)
				print(f"Created output directory: {self.roundInput.elongation.roundDir}")
			except Exception as e:
				print(f"Failed to create output directory: {e}")
		
		if os.path.exists(self.extensionReads)==True and os.path.getsize(self.extensionReads)!=0:
			self.readlog()
		else:
			logfilet=open(self.log,'w')
			try:
				self.potentialExtensionReadsAln,self.minimap2Command,self.minimap2Output=self.minimap2()
				logLine='potentialExtensionReadsAln\t'+self.potentialExtensionReadsAln+"\nminimap2Command\t"+self.minimap2Command+"\nminimap2Output\t"+self.minimap2Output+"\n"
				logfilet.writelines(logLine)

				self.minimumExtensionReads()
				
				# Check if minimumExtensionReads successfully set necessary attributes
				if hasattr(self, 'minimumThresholdReadsAln'):
					logLine='minimumThresholdReadsAln\t'+self.minimumThresholdReadsAln+"\n"
					logfilet.writelines(logLine)
				
				if hasattr(self, 'minimumThresholdReadsID'):
					logLine='minimumThresholdReadsID\t'+';'.join(self.minimumThresholdReadsID)+"\n"
					logfilet.writelines(logLine)
				
				if hasattr(self, 'minimumThresholdExtensionReadsAln'):
					logLine='minimumThresholdExtensionReadsAln\t'+self.minimumThresholdExtensionReadsAln+"\n"
					logfilet.writelines(logLine)
				
				if hasattr(self, 'minimumThresholdExtensionReads'):
					logLine='minimumThresholdExtensionReads\t'+self.minimumThresholdExtensionReads+"\n"
					logfilet.writelines(logLine)
				
				if hasattr(self, 'minimumThresholdExtensionReadsID'):
					logLine='minimumThresholdExtensionReadsID\t'+';'.join(self.minimumThresholdExtensionReadsID)+"\n"
					logfilet.writelines(logLine)
			except Exception as e:
				import traceback
				print(f"Error during initialization: {e}")
				traceback.print_exc()
				self.note = 'initializationError'
				logLine='note\t'+self.note+"\n"
				logfilet.writelines(logLine)
				logfilet.close()
				return
			
			
			if self.note=='':
				self.selectMappingQuality=20
				self.selectAlignmentLength=3000
				self.selectNMAlignmentLengthratio=0.1
				
				self.selectReadsNum=0
				self.extensionReadsNum=0
				self.readsExtensionLength=1000
				self.extensionReadsEdge=10
				'''while self.extensionReadsNum<=0 and self.note=='':
					try:
						self.selectPotentialExtensionReadsAln=self.roundInput.elongation.roundDir+"/selectPotentialExtensionReads."+self.roundInput.elongation.base.tag+".bam"
						#self.selectPotentialExtensionReadsID=self.samFilter(self.potentialExtensionReadsAln,self.selectPotentialExtensionReadsAln)
						self.selectPotentialExtensionReadsID=self.samFilter(self.minimumThresholdExtensionReadsAln,self.selectPotentialExtensionReadsAln)
						self.selectReadsNum=len(self.selectPotentialExtensionReadsID)

						self.extensionReadsAln=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".bam"
						self.extensionReads=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".fa"
						self.extensionReadsID=self.extensionFinder(self.selectPotentialExtensionReadsAln,self.extensionReadsAln,self.extensionReads)
						self.extensionReadsNum=len(self.extensionReadsID)

						if self.extensionReadsNum==0:
							if self.selectMappingQuality==0:
								if self.readsExtensionLength<=10:
									if self.selectAlignmentLength<=500:
										# Fix case where normal stopping doesn't occur when reduced to minimum threshold

										# self.extensionReadsEdge=self.extensionReadsEdge+10
										# self.selectMappingQuality=20
										# self.selectAlignmentLength=3000
										# self.readsExtensionLength=1000

										# When all parameters reach most relaxed state, don't reset but set failure flag and break loop
										self.note = 'No extension reads found after relaxing all parameters'
										print(f"Warning: {self.note}")
										break # Use break statement to force exit from while loop
									else:
										self.selectAlignmentLength=self.selectAlignmentLength-100
								else:
									self.readsExtensionLength=self.readsExtensionLength-100
							else:
								self.selectMappingQuality=self.selectMappingQuality-10
					except Exception as e:
						import traceback
						print(f"Error finding extension reads: {e}")
						traceback.print_exc()
						self.note = 'extensionFinderError'
						break'''
				# ... (after initializing self.extensionReadsEdge=10)
				# Get the maximum edge value set from command line, as the upper limit for the loop
				max_edge = self.roundInput.elongation.base.edge # value is 500
				while self.extensionReadsNum <= 0 and self.note == '':
					try:
						# 1. Execute filtering
						# First call samFilter, then call extensionFinder
						self.selectPotentialExtensionReadsAln = self.roundInput.elongation.roundDir + "/selectPotentialExtensionReads." + self.roundInput.elongation.base.tag + ".bam"
						self.selectPotentialExtensionReadsID = self.samFilter(self.minimumThresholdExtensionReadsAln, self.selectPotentialExtensionReadsAln)
						
						self.extensionReadsAln = self.roundInput.elongation.roundDir + "/extensionReads." + self.roundInput.elongation.base.tag + ".bam"
						self.extensionReads = self.roundInput.elongation.roundDir + "/extensionReads." + self.roundInput.elongation.base.tag + ".fa"
						self.extensionReadsID = self.extensionFinder(self.selectPotentialExtensionReadsAln, self.extensionReadsAln, self.extensionReads)
						self.extensionReadsNum = len(self.extensionReadsID)

						# 2. If no reads found, update parameters
						if self.extensionReadsNum == 0:
							print(f"No extension reads found with params: MQ={self.selectMappingQuality}, "
								f"AlnLen={self.selectAlignmentLength}, extLen={self.readsExtensionLength}, edge={self.extensionReadsEdge}")

							# Prioritize relaxing quality and length parameters
							if self.selectMappingQuality > 0:
								self.selectMappingQuality -= 10
							elif self.readsExtensionLength > 10:
								self.readsExtensionLength -= 100
							elif self.selectAlignmentLength > 500:
								self.selectAlignmentLength -= 100
							
							# Only start relaxing edge when all other parameters are at their most relaxed
							elif self.extensionReadsEdge < max_edge:
								self.extensionReadsEdge += 10
							
							# If all parameters (including edge) have reached their most relaxed state, terminate the loop
							else:
								self.note = f'No extension reads found even after relaxing edge to {self.extensionReadsEdge}'
								print(f"Warning: {self.note}")
								break # Successfully exit the loop

					except Exception as e:
						import traceback
						print(f"Error finding extension reads: {e}")
						traceback.print_exc()
						self.note = 'extensionFinderError'
						break
				logLine='selectPotentialExtensionReadsAln\t'+self.selectPotentialExtensionReadsAln+"\n"
				logLine+='selectPotentialExtensionReadsID\t'+';'.join(self.selectPotentialExtensionReadsID)+"\n"
				logLine+='selectReadsNum\t'+str(self.selectReadsNum)+"\n"
				logLine+='extensionReadsAln\t'+self.extensionReadsAln+"\n"
				logLine+='extensionReads\t'+self.extensionReads+"\n"
				logLine+='extensionReadsID\t'+';'.join(self.extensionReadsID)+"\n"
				logLine+='extensionReadsNum\t'+str(self.extensionReadsNum)+"\n"
				logfilet.writelines(logLine)

			else:
				logLine='selectPotentialExtensionReadsAln\tNone\n'
				logLine+='selectPotentialExtensionReadsID\tNone\n'
				logLine+='selectReadsNum\t0\n'
				logLine+='extensionReadsAln\tNone\n'
				logLine+='extensionReads\tNone\n'
				logLine+='extensionReadsID\tNone\n'
				logLine+='extensionReadsNum\t0\n'
				logfilet.writelines(logLine)
			logLine='selectMappingQuality\t'+str(self.selectMappingQuality)+"\n"
			logLine+='selectAlignmentLength\t'+str(self.selectAlignmentLength)+"\n"
			logLine+='selectNMAlignmentLengthratio\t'+str(self.selectNMAlignmentLengthratio)+"\n"
			logLine+='readsExtensionLength\t'+str(self.readsExtensionLength)+"\n"
			logLine+='extensionReadsEdge\t'+str(self.extensionReadsEdge)+"\n"
			logLine+='note\t'+str(self.note)+"\n"
			logfilet.writelines(logLine)	

			logfilet.close()
	
	def minimumExtensionReads(self):
		# Import required modules
		import os
		import traceback
		
		self.selectMappingQuality=0
		self.selectAlignmentLength=500
		self.selectNMAlignmentLengthratio=0.1
		self.extensionReadsEdge=self.roundInput.elongation.base.edge
		self.readsExtensionLength=10

		self.minimumThresholdReadsAln=self.roundInput.elongation.roundDir+"/minimumThresholdReads."+self.roundInput.elongation.base.tag+".bam"
		self.minimumThresholdReadsID=self.samFilter(self.potentialExtensionReadsAln,self.minimumThresholdReadsAln)

		self.minimumThresholdExtensionReadsAln=self.roundInput.elongation.roundDir+"/minimumThresholdExtensionReads."+self.roundInput.elongation.base.tag+".bam"
		self.minimumThresholdExtensionReads=self.roundInput.elongation.roundDir+"/minimumThresholdExtensionReads."+self.roundInput.elongation.base.tag+".fa"
		
		# Check if input file exists
		if not os.path.exists(self.minimumThresholdReadsAln):
			print(f"Error: Input file does not exist: {self.minimumThresholdReadsAln}")
			self.minimumThresholdExtensionReadsID = []
			self.note = 'inputFileNotFound'
			return
		
		try:
			print("Calling extensionFinder to find minimum threshold extension reads...")
			self.minimumThresholdExtensionReadsID=self.extensionFinder(
				self.minimumThresholdReadsAln,
				self.minimumThresholdExtensionReadsAln,
				self.minimumThresholdExtensionReads
			)
			print(f"Found {len(self.minimumThresholdExtensionReadsID)} minimum threshold extension reads")
		except Exception as e:
			print(f"Error finding minimum threshold extension reads: {str(e)}")
			traceback.print_exc()
			self.minimumThresholdExtensionReadsID = []
			print("Will continue processing, but may affect result quality")

		if not hasattr(self, 'minimumThresholdExtensionReadsID') or len(self.minimumThresholdExtensionReadsID)==0:
			self.note='noExtensionReadsFoundAtMinimumThreshold'
			print(f"Setting note to: {self.note}")
		else:
			self.note=''

# Modified extensionFinder method in FindExtensionReads.py
	def extensionFinder(self, inputAln, outputAln, outputSeq):
			readslist = []
			outputSeqFile = open(outputSeq, 'w')
			inputAlnFile = AlignmentFile(inputAln, "rb", check_sq=False)
			outputAlnFile = AlignmentFile(outputAln, "wb", template=inputAlnFile)
			
			# ===== Core modification: restore 1.0's direct memory access method =====
			for r in inputAlnFile:
					queryID = r.query_name
					
					# Directly get sequence from memory dictionary (restore 1.0 logic)
					try:
							query = self.roundInput.elongation.base.readsDict[queryID]
					except KeyError:
							print(f"Warning: Cannot find ID {queryID} in readsDict, skipping this record")
							continue
					
					# Restore 1.0's sequence processing logic
					if r.is_reverse:
							queryseq = query.seq.reverse_complement()
					else:
							queryseq = query.seq
					
					# Calculate boundary distance (maintain 1.0 logic)
					queryDistance, refDistance, extensionLength, extensionReadsSeq = self.calculateBoundDistance(
							queryID, query, r
					)
					
					# Condition judgment (maintain 1.0 logic)
					if self.readsExtensionLength <= extensionLength:
							if queryDistance <= self.extensionReadsEdge and refDistance <= self.extensionReadsEdge:
									if queryID not in self.usedReads:
											if queryID not in readslist:
													l = '>' + queryID + "\n" + str(queryseq) + "\n"
													outputSeqFile.writelines(l)
													readslist.append(queryID)
											outputAlnFile.write(r)
			# ===== End of modification =====
			
			outputSeqFile.close()
			inputAlnFile.close()
			outputAlnFile.close()
			return readslist	
	def calculateBoundDistance(self,queryID,query,r):
		if r.is_reverse:
			queryseq=query.seq.reverse_complement()
		else:
			queryseq=query.seq
		
		qs,qe,rs,re=self.findAlnPosition(queryseq,self.roundInput.inputSeedSequence,r)
		if self.roundInput.elongation.base.flag=='left':
			queryDistance=qs
			extensionReadsSeq=queryseq
			refDistance=len(self.roundInput.inputSeedSequence.seq)-re
			extensionLength=len(queryseq)-qe
		else:
			queryDistance=len(queryseq)-qe
			extensionReadsSeq=queryseq
			refDistance=rs
			extensionLength=qs
		if qs<0 or len(queryseq)<qe or qe<0 or rs<0 or re<0 :#or len(self.roundInput.inputSeedSequence.seq)<re:
			print ('wrong sequence end',qs,qe,rs,re)
			print (queryDistance,refDistance,extensionLength,extensionReadsSeq)
			print (qs<0)
			print (len(queryseq)<qe)
			print (qe<0)
			print (rs<0)
			print (re<0)
			print (len(self.roundInput.inputSeedSequence.seq)<re)
			print (len(self.roundInput.inputSeedSequence.seq),re)
			sys.exit()
		return queryDistance,refDistance,extensionLength,extensionReadsSeq
	
	def findAlnPosition(self,queryseq,refseq,r):
		ct=r.cigartuples
		scpstart=ct[0]
		if scpstart[0]==4:
			scps=scpstart[1]
		else:
			scps=0
	
		scpend=ct[-1]
		if scpend[0]==4:
			scpe=scpend[1]
		else:
			scpe=0

		if queryseq==r.query_sequence or r.query_alignment_sequence==None:
			qs=r.query_alignment_start
			qe=r.query_alignment_end
		else:
			qaln=r.query_alignment_sequence
			qs=str(queryseq).index(qaln)
			qe=qs+len(r.query_alignment_sequence)
		rs=r.reference_start
		re=r.reference_end
		if qs<0:
			print ('wrong query start',qs)
		if qe>len(queryseq):
			print ('wrong query end',qe,len(queryseq))
		return qs,qe,rs,re

		

	def samFilter(self,inputAln,outputAln):
		samFile=AlignmentFile(inputAln,"rb",check_sq=False)
		outputBamFile=AlignmentFile(outputAln,"wb", template=samFile)
		
		readslist=[]
		print (inputAln)
		sMQ=self.selectMappingQuality
		sAlignmentLength=self.selectAlignmentLength
		sNMAlignmentLengthr=self.selectNMAlignmentLengthratio

		for r in samFile:
			if r.is_unmapped==False and r.mapping_quality>=sMQ:
				for i in r.tags:
					if i[0]=='NM':
						NM=i[1]
				AlignmentLength=len(r.query_alignment_sequence)
				if  AlignmentLength>=sAlignmentLength and float(NM)/AlignmentLength<=sNMAlignmentLengthr:
					outputBamFile.write(r)
					if r.query_name not in readslist:
						readslist.append(r.query_name)
		samFile.close()
		outputBamFile.close()
		return readslist
	
	def jellyfish(self, kmerseq, kmer_size, kmer_num, kmer_length=1000):
		"""
		Extract kmers from sequence based on flag parameter
		First extract a subsequence of length kmer_length from one end of the sequence:
		- If flag=right, extract from the left end of the sequence
		- If flag=left, extract from the right end of the sequence
		Then sample kmer_num kmers from the subsequence

		Parameters:
			kmerseq: Sequence to extract k-mers from
			kmer_size: Length of k-mer
			kmer_num: Number of k-mers to extract
			kmer_length: Length of subsequence for k-mer extraction, usually a proportion of average read length

		Returns:
			File path containing k-mers
		"""
		# Import required modules
		import os
		
		# Create output directory (if it doesn't exist)
		output_dir = self.roundInput.elongation.roundDir
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)

		# Record original sequence length
		seq_length = len(kmerseq)
		print(f"Original sequence length: {seq_length} bp")

		# If sequence is too short, cannot extract k-mer
		if seq_length < kmer_size:
			print(f"Warning: Sequence length ({seq_length}) is smaller than k-mer size ({kmer_size}), cannot extract k-mer")
			return None
		
		# Extract subsequence
		sub_seq = ""
		if self.roundInput.elongation.base.flag == 'right':
			# Extract subsequence of length kmer_length from left end
			sub_seq_length = min(kmer_length, seq_length)
			sub_seq = kmerseq[:sub_seq_length]
			print(f"Extracted subsequence of length {sub_seq_length}bp from left end for k-mer extraction")
		else:
			# Extract subsequence of length kmer_length from right end
			sub_seq_length = min(kmer_length, seq_length)
			sub_seq = kmerseq[-sub_seq_length:]
			print(f"Extracted subsequence of length {sub_seq_length}bp from right end for k-mer extraction")

		# Ensure subsequence length is sufficient for k-mer extraction
		if len(sub_seq) < kmer_size:
			print(f"Warning: Subsequence length ({len(sub_seq)}) is smaller than k-mer size ({kmer_size}), cannot extract k-mer")
			return None

		# If subsequence length is smaller than expected kmer_length, adjust sampling strategy
		if len(sub_seq) < kmer_length:
			print(f"Note: Subsequence length ({len(sub_seq)}bp) is smaller than expected kmer_length ({kmer_length}bp), will use entire sequence and adjust sampling strategy")
		
		# Generate sampling points
		if self.roundInput.elongation.base.flag == 'left':
			# Applicable for starting from right end
			def generate_points(n):
				u = [i / n for i in range(1, n+1)]  # Generate u_i = 1/n, 2/n, ..., 1
				x = [num**0.5 for num in u]         # Convert to x_i = sqrt(u_i)
				x.reverse()                         # Reverse order, first point is right endpoint 1
				return x
		else:
			# Applicable for starting from left end
			def generate_points(n):
				u = [i / n for i in range(1, n+1)]  # Generate u_i = 1/n, 2/n, ..., 1
				x = [num**0.5 for num in u]
				y = []
				for i in x:
					i = 1-i
					y.append(i)
				y.reverse()
				return y
		
		# Generate sampling point positions
		points = generate_points(kmer_num)

		# Extract k-mers
		kmers = []

		# Adjust sampling strategy to ensure sufficient k-mer extraction even when sequence is short
		available_positions = len(sub_seq) - kmer_size + 1
		if available_positions < kmer_num:
			print(f"Warning: Available positions ({available_positions}) is less than requested k-mer count ({kmer_num}), will use uniform sampling")
			# Uniform sampling
			if available_positions <= 1:
				# Only one position available, use directly
				kmers.append((sub_seq[:kmer_size], 0))
			else:
				# Uniformly distribute sampling points
				step = available_positions / min(available_positions, kmer_num)
				for i in range(min(available_positions, kmer_num)):
					pos = int(i * step)
					kmer = sub_seq[pos:pos + kmer_size]
					kmers.append((kmer, pos))
		else:
			# Normal sampling
			for point in points:
				# Map points in 0-1 interval to subsequence length
				pos = int(point * (len(sub_seq) - kmer_size))

				# Ensure position is valid
				pos = max(0, min(pos, len(sub_seq) - kmer_size))

				if self.roundInput.elongation.base.flag == 'left':
					# Extract kmer_size length leftward from position
					start = max(0, pos - kmer_size + 1)
					kmer = sub_seq[start:start + kmer_size]
				else:
					# Extract kmer_size length rightward from position
					kmer = sub_seq[pos:pos + kmer_size]

				kmers.append((kmer, pos))

		print(f"Generated a total of {len(kmers)} k-mers")
		
		# Write to output file, ensuring k-mer sequences are uppercase
		kmer_output = os.path.join(output_dir, "kmers.txt")
		with open(kmer_output, 'w') as f:
			for kmer, pos in kmers:
				# Convert k-mer to uppercase before writing to file
				uppercase_kmer = str(kmer).upper()
				f.write(f"{uppercase_kmer}\n")

				# Calculate relative position for output
				rel_pos = pos / len(sub_seq)
				rel_pos_percent = rel_pos * 100

				print(f"Selected k-mer: {uppercase_kmer} (position: {pos}, relative position: {rel_pos_percent:.1f}%)")

		# If extracted k-mer count is insufficient, output warning
		if len(kmers) < kmer_num:
			print(f"Warning: Only extracted {len(kmers)} k-mers, less than requested {kmer_num}")
		
		return str(kmer_output)

	def seqkit(self, kmer_list, reads):
		"""
		Call seqkit grep command to extract sequences containing specified k-mers from kmer_list in reads
		Use GNU parallel to process split reads files in parallel, improving processing speed

		Parameters:
			kmer_list: File path containing k-mers
			reads: Reads file path to search

		Returns:
			Output file path
		"""
		# Import required modules
		import subprocess
		import os
		import glob
		import time
		import shutil
		
		# Ensure using absolute paths
		kmer_list = os.path.abspath(kmer_list)
		reads = os.path.abspath(reads)

		print(f"seqkit function using absolute paths: kmer_list={kmer_list}, reads={reads}")

		# Create output directory (if it doesn't exist)
		output_dir = self.roundInput.elongation.roundDir
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)

		# Create temporary output directory
		seqkit_tmp_dir = os.path.join(output_dir, "seqkit_tmp")
		if not os.path.exists(seqkit_tmp_dir):
			os.makedirs(seqkit_tmp_dir)

		# Set final output file path
		output_file = os.path.join(output_dir, "seqkitOutput.fa")

		# Get split reads file paths
		# Use out parameter passed during initialization to get DEGAP.py output directory
		if self.out:
			# If out parameter was provided during initialization, use it first
			base_out_dir = os.path.abspath(self.out)
		else:
			# Try to get from other attributes
			base_out_dir = None
			if hasattr(self.roundInput.elongation.base, "out"):
				base_out_dir = os.path.abspath(self.roundInput.elongation.base.out)
			elif hasattr(self.roundInput.elongation, "out"):
				base_out_dir = os.path.abspath(self.roundInput.elongation.out)

			# If still unable to get, use current working directory as fallback
			if not base_out_dir:
				base_out_dir = os.path.abspath(os.getcwd())
				print(f"Warning: Could not determine base output directory, using current directory: {base_out_dir}")

		# Build reads_part directory path
		reads_part_dir = os.path.join(base_out_dir, "reads_part")
		print(f"Using reads_part directory: {reads_part_dir}")

		# Find split files
		reads_files = glob.glob(f"{reads_part_dir}/*.fasta")
		
		if not reads_files:
			reads_files = glob.glob(f"{reads_part_dir}/*.fa*")

		if not reads_files:
			# If split reads files are not found, use original method
			print(f"Warning: No split read files found in {reads_part_dir}, using original read file")
			cmd = f"seqkit grep -f {kmer_list} -s {reads} > {output_file}"

			try:
				process = subprocess.run(cmd, shell=True, check=True)
				print(f"Seqkit completed with single file")
				return str(output_file)
			except subprocess.CalledProcessError as e:
				print(f"Seqkit command failed: {e}")
				raise
		
		print(f"Found {len(reads_files)} split read files in {reads_part_dir}")

		# Clean up possibly existing old output files
		for old_file in glob.glob(f"{seqkit_tmp_dir}/*.output.fasta"):
			try:
				os.remove(old_file)
			except:
				pass

		# Create log file for parallel tasks
		parallel_log = os.path.join(output_dir, "parallel_seqkit.log")

		# 1. Use GNU parallel to process split reads files in parallel, capturing output and exit status
		print("Starting parallel processing of split read files...")
		try:
			# Use --joblog parameter to record job status
			parallel_cmd = f"parallel --joblog {output_dir}/parallel_jobs.log -j {self.j} 'seqkit grep -f {kmer_list} -s {{}} -o {seqkit_tmp_dir}/{{/.}}.output.fasta' ::: {reads_part_dir}/*.fa*"
			
			print(f"Running parallel seqkit with command: {parallel_cmd}")
			# Capture standard output and error output
			result = subprocess.run(
				parallel_cmd,
				shell=True,
				check=True,
				stdout=subprocess.PIPE,
				stderr=subprocess.PIPE,
				text=True
			)

			# Print command output information
			if result.stdout:
				print("Parallel command output:")
				print(result.stdout[:1000] + "..." if len(result.stdout) > 1000 else result.stdout)

			if result.stderr:
				print("Parallel command warnings/errors:")
				print(result.stderr[:1000] + "..." if len(result.stderr) > 1000 else result.stderr)

			print(f"GNU parallel command completed with return code {result.returncode}")

			# Check parallel-generated joblog file to confirm all tasks are completed
			if os.path.exists(f"{output_dir}/parallel_jobs.log"):
				with open(f"{output_dir}/parallel_jobs.log", 'r') as log_file:
					lines = log_file.readlines()
					# Skip header line
					data_lines = [line for line in lines if line.strip() and not line.startswith("#")]
					completed_jobs = len(data_lines)
					print(f"Parallel jobs completed: {completed_jobs}/{len(reads_files)}")

					# Modify logic for determining task failure: exitval=0 means success
					failed_jobs = []
					for line in data_lines:
						cols = line.split("\t")
						if len(cols) >= 7:  # Ensure at least 7 columns
							exitval = cols[6].strip()
							if exitval != "0":  # Only consider failure when exitval is not 0
								failed_jobs.append(line)

					if failed_jobs:
						print(f"Warning: {len(failed_jobs)} parallel jobs failed with non-zero exit codes")
						# Can add more detailed failure reason analysis
						for i, job in enumerate(failed_jobs[:3]):  # Only show first 3 failed tasks
							print(f"  Failed job {i+1}: {job.strip()}")
						if len(failed_jobs) > 3:
							print(f"  ... and {len(failed_jobs)-3} more failures")
			
		except subprocess.CalledProcessError as e:
			print(f"Parallel seqkit command failed with return code {e.returncode}")
			print(f"Error output: {e.stderr if hasattr(e, 'stderr') else 'No error output available'}")
			if os.path.exists(f"{output_dir}/parallel_jobs.log"):
				print("Checking job log for partial results...")
				# Even if command fails, we try to process possible partial results
			else:
				print("No job log found, parallel execution may have failed completely")
				raise

		# Wait a short time to ensure all files are written completely
		time.sleep(5)

		# 2. Merge results
		output_files = glob.glob(f"{seqkit_tmp_dir}/*.output.fasta")

		if not output_files:
			print("Warning: No output files were generated from parallel seqkit")
			# Possibly no matches, create an empty file
			open(output_file, 'w').close()
			return str(output_file)

		print(f"Merging {len(output_files)} output files...")
		
		# Use cat command to merge files
		try:
			merge_cmd = f"cat {seqkit_tmp_dir}/*.output.fasta > {output_file}"
			print(f"Merging results with command: {merge_cmd}")
			process = subprocess.run(merge_cmd, shell=True, check=True)
			print(f"Successfully merged {len(output_files)} output files")
		except subprocess.CalledProcessError as e:
			print(f"Merge command failed: {e}")

			# If cat command fails, try using Python to merge files
			try:
				print("Attempting to merge files using Python...")
				with open(output_file, 'w') as outf:
					for f in output_files:
						with open(f, 'r') as inf:
							shutil.copyfileobj(inf, outf)
				print("Python file merge successful")
			except Exception as e2:
				print(f"Python merge also failed: {e2}")
				raise

		# 3. Clean up intermediate files
		try:
			subprocess.run(f"rm -r {seqkit_tmp_dir}", shell=True, check=True)
			print(f"Cleaned up temporary files in {seqkit_tmp_dir}")
		except subprocess.CalledProcessError as e:
			print(f"Warning: Failed to clean up temporary files: {e}")

		# Check output file
		if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
			print(f"Seqkit completed successfully, output file size: {os.path.getsize(output_file)} bytes")
		else:
			print(f"Warning: Seqkit output file is empty or doesn't exist")
		
		return str(output_file)

	def seqkitfliter(self, inputfasta, inputkmer):
		"""
		Filter reads in large files, keeping reads containing low-frequency kmers

		Parameters:
			inputfasta: Input fasta file path
			inputkmer: File path containing kmers

		Returns:
			Filtered file path or original file path (if file is small)
		"""
		# Import required modules
		import subprocess
		import os
		from Bio import SeqIO
		import collections

		# Ensure using absolute paths
		inputfasta = os.path.abspath(inputfasta)
		inputkmer = os.path.abspath(inputkmer)

		print(f"seqkitfliter function using absolute paths: inputfasta={inputfasta}, inputkmer={inputkmer}")

		# Create output directory (if it doesn't exist)
		output_dir = self.roundInput.elongation.roundDir
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)

		# Set output file path
		output_file = os.path.join(output_dir, "seqkitOutputFliter.fa")

		# 1. Check inputfasta file size
		fasta_size = os.path.getsize(inputfasta) / (1024 * 1024 * 1024)  # Convert to GB
		print(f"Input file size: {fasta_size:.2f} GB")

		# If file is smaller than 1GB, return original file directly
		if fasta_size < 1.0:
			print(f"File size is smaller than 1GB, no filtering needed, returning original file directly")
			return inputfasta
		
		# 2. Read kmer file
		kmers = []
		with open(inputkmer, 'r') as f:
			for line in f:
				kmer = line.strip()
				if kmer:
					kmers.append(kmer)

		print(f"Read {len(kmers)} kmers")

		# Count occurrences of each kmer in fasta file
		kmer_counts = collections.defaultdict(int)

		# Use seqkit locate command to find kmer occurrence positions and counts
		for kmer in kmers:
			try:
				cmd = f"seqkit locate -p {kmer} {inputfasta} -m 0"
				result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

				# Count result lines (minus header line)
				lines = result.stdout.strip().split('\n')
				count = len(lines) - 1 if lines and lines[0] else 0
				kmer_counts[kmer] = count

				print(f"kmer: {kmer} occurrence count: {count}")

			except subprocess.CalledProcessError as e:
				print(f"Error finding kmer {kmer}: {e}")
				kmer_counts[kmer] = 0
		
		# Dynamic threshold adjustment logic - performed after all kmer statistics are completed
		thresholds = [50, 100]

		# First round check for kmers ≤50 times
		current_threshold = thresholds[0]
		low_freq_kmers = [kmer for kmer, count in kmer_counts.items() if count <= current_threshold and count > 0]

		print(f"Found {len(low_freq_kmers)} low-frequency kmers (occurrence count ≤ {current_threshold})")

		# Check if quantity requirement is met
		if len(low_freq_kmers) < 0.5 * len(kmers):
			print(f"Found low-frequency kmers insufficient {0.5 * len(kmers)} ({len(low_freq_kmers)}), relaxing threshold to {thresholds[1]}")
			current_threshold = thresholds[1]
			low_freq_kmers = [kmer for kmer, count in kmer_counts.items() if count <= current_threshold and count > 0]
			print(f"Found {len(low_freq_kmers)} low-frequency kmers (occurrence count ≤ {current_threshold})")

		if not low_freq_kmers:
			print("No qualifying low-frequency kmers found, returning original file")
			return inputfasta

		# 3. Create temporary file to save low-frequency kmers
		low_freq_kmer_file = os.path.join(output_dir, "low_freq_kmers.txt")
		with open(low_freq_kmer_file, 'w') as f:
			for kmer in low_freq_kmers:
				f.write(f"{kmer}\n")

		# Use seqkit grep to extract reads containing low-frequency kmers
		try:
			cmd = f"seqkit grep -f {low_freq_kmer_file} -s {inputfasta} -o {output_file}"
			print(f"Executing command: {cmd}")
			subprocess.run(cmd, shell=True, check=True)

			# Check output file
			if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
				print(f"Filtering completed, output file size: {os.path.getsize(output_file) / (1024*1024):.2f} MB")
				# 4. Return output file path
				return str(output_file)
			else:
				print(f"No valid results after filtering, returning original file")
				return inputfasta

		except Exception as e:
			print(f"Error during filtering process: {e}")
			return inputfasta

		# Clean up temporary files
		try:
			os.remove(low_freq_kmer_file)
		except:
			pass

	def blast(self, kmer_list, reads):
		"""
		Use blastn command to align kmer_list with reads and extract matching read sequences

		Parameters:
			kmer_list: File path containing k-mers
			reads: Reads file path to search

		Returns:
			Output file path
		"""
		# Import required modules
		import subprocess
		import os
		from Bio import SeqIO
		import glob
		import time
		import shutil

		# Ensure using absolute paths
		kmer_list = os.path.abspath(kmer_list)
		reads = os.path.abspath(reads)

		print(f"blast function using absolute paths: kmer_list={kmer_list}, reads={reads}")

		# Create output directory (if it doesn't exist)
		output_dir = self.roundInput.elongation.roundDir
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)

		# Set output file paths
		blast_out = os.path.join(output_dir, "blastOutput.txt")
		output_file = os.path.join(output_dir, "blastOutput.fa")

		# Create FASTA file containing kmers (blastn requires FASTA format input)
		kmer_fasta = os.path.join(output_dir, "kmers.fasta")
		with open(kmer_list, 'r') as kmer_in, open(kmer_fasta, 'w') as kmer_out:
			kmer_count = 0
			for line in kmer_in:
				kmer = line.strip()
				if kmer:
					kmer_out.write(f">kmer_{kmer_count}\n{kmer}\n")
					kmer_count += 1

		# Check input file size
		reads_size = os.path.getsize(reads) / (1024 * 1024 * 1024)  # Convert to GB
		print(f"Input file size: {reads_size:.2f} GB")
		
		# If file is larger than 1GB, need to split for processing
		if reads_size > 1.0:
			print(f"Detected large file ({reads_size:.2f} GB), will perform split processing")

			# Create temporary directory to save split files
			split_dir = os.path.join(output_dir, "blast_split_tmp")
			if not os.path.exists(split_dir):
				os.makedirs(split_dir)

			# Use seqkit split to split file
			try:
				# 10000 sequences per file
				split_cmd = f"seqkit split -s 10000 {reads} -O {split_dir} --force"
				print(f"Executing split command: {split_cmd}")
				subprocess.run(split_cmd, shell=True, check=True)

				# Get list of split files
				split_files = glob.glob(f"{split_dir}/*.fa*")
				print(f"Successfully split into {len(split_files)} small files")

				# Create temporary directory to save blast results
				blast_tmp_dir = os.path.join(output_dir, "blast_results_tmp")
				if not os.path.exists(blast_tmp_dir):
					os.makedirs(blast_tmp_dir)

				# Run blastn on each small file
				matched_ids = set()
				for i, split_file in enumerate(split_files):
					split_base = os.path.basename(split_file)
					blast_result = os.path.join(blast_tmp_dir, f"{split_base}.blast.txt")

					print(f"Processing split file {i+1}/{len(split_files)}: {split_base}")
					try:
						# Execute blastn command
						blast_cmd = f"blastn -query {kmer_fasta} -subject {split_file} -outfmt 6 -out {blast_result}"
						subprocess.run(blast_cmd, shell=True, check=True)

						# Collect matched IDs
						if os.path.exists(blast_result) and os.path.getsize(blast_result) > 0:
							with open(blast_result, 'r') as f:
								for line in f:
									fields = line.strip().split('\t')
									if len(fields) >= 2:
										matched_ids.add(fields[1])  # Second column is subject ID

					except subprocess.CalledProcessError as e:
						print(f"Warning: blastn failed when processing file {split_file}: {e}")
						# Continue processing other files

				print(f"All split files processed, found {len(matched_ids)} matching reads IDs")

				# Extract matching sequences from original reads file
				with open(output_file, 'w') as out_f:
					for record in SeqIO.parse(reads, "fasta"):
						if record.id in matched_ids:
							SeqIO.write(record, out_f, "fasta")

				# Clean up temporary files
				try:
					import shutil
					shutil.rmtree(split_dir)
					shutil.rmtree(blast_tmp_dir)
					print("Temporary file cleanup completed")
				except Exception as e:
					print(f"Warning: Failed to clean up temporary files: {e}")

			except Exception as e:
				print(f"Split processing error: {e}")
				# Return seqkit result on failure
				print("Due to blastn processing failure, will use seqkit result directly")
				return reads
		else:
			# File is small, process directly
			try:
				# Execute blastn command
				cmd = f"blastn -query {kmer_fasta} -subject {reads} -outfmt 6 -out {blast_out}"
				print(f"Executing blastn command: {cmd}")
				subprocess.run(cmd, shell=True, check=True)

				# Extract matching reads IDs from blast results
				matched_ids = set()
				with open(blast_out, 'r') as f:
					for line in f:
						fields = line.strip().split('\t')
						if len(fields) >= 2:
							matched_ids.add(fields[1])  # Second column is subject ID

				# Extract matching sequences from reads file
				with open(output_file, 'w') as out_f:
					for record in SeqIO.parse(reads, "fasta"):
						if record.id in matched_ids:
							SeqIO.write(record, out_f, "fasta")

			except Exception as e:
				print(f"blastn processing error: {e}")
				# Return seqkit result on failure
				print("Due to blastn processing failure, will use seqkit result directly")
				return reads

		# Clean up temporary files
		try:
			os.remove(kmer_fasta)
		except:
			pass

		# Check output file
		if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
			print(f"blastn filtering completed, output file size: {os.path.getsize(output_file) / (1024*1024):.2f} MB")
			return str(output_file)
		else:
			print(f"blastn did not generate valid results, will use seqkit result")
			return reads

	def kmerfilter(self, inputSeq, kmer_size, kmer_num):
		"""
		Use k-mer filtering method to screen HiFi reads

		Parameters:
			inputSeq: Input sequence file
			kmer_size: Size of k-mer
			kmer_num: Number of k-mers needed

		Returns:
			Filtered reads file path
		"""
		# Import required modules
		import os
		import subprocess
		import math
		from Bio import SeqIO

		print(f"Executing k-mer filtering, parameters: kmer_size={kmer_size}, kmer_num={kmer_num}, kmer_length={self.kmer_length}")

		# Get DEGAP.py parameters
		is_right_flag = self.roundInput.elongation.base.flag == 'right'

		# Extract sequence based on flag parameter
		kmerseq = ""
		for seq_record in SeqIO.parse(inputSeq, "fasta"):
			seq = str(seq_record.seq)

			if is_right_flag:
				# Extract sequence from left end
				kmerseq = seq
			else:
				# Extract sequence from right end
				kmerseq = seq

			# Only process first sequence
			break

		print(f"Extracted sequence length: {len(kmerseq)} bp, used for k-mer generation")

		# Use jellyfish to generate k-mer list - ensure parameters are passed correctly
		# Important: pass the same kmer_size, kmer_num and kmer_length parameters
		kmer_list = self.jellyfish(
			kmerseq=kmerseq,
			kmer_size=kmer_size,
			kmer_num=kmer_num,
			kmer_length=self.kmer_length
		)

		if not kmer_list:
			print("Warning: jellyfish failed to generate valid k-mer list")
			return None

		# Ensure using absolute path
		reads_path = os.path.abspath(self.roundInput.elongation.base.reads)
		print(f"Using absolute path to access reads file: {reads_path}")

		# Use seqkit to screen reads containing these k-mers
		try:
			seqkitOutput = self.seqkit(kmer_list, reads_path)

			# Check seqkitOutput size, if larger than 1GB call blast function for secondary filtering
			if os.path.exists(seqkitOutput):
				output_size = os.path.getsize(seqkitOutput) / (1024 * 1024 * 1024)  # Convert to GB
				if output_size > 1.0:
					print(f"seqkitOutput size is {output_size:.2f} GB, larger than 1GB, calling blast function for secondary filtering")
					blastOutput = self.blast(kmer_list, seqkitOutput)

					# Check blastOutput size, if larger than 1GB call new blast_filter function for tertiary filtering
					if os.path.exists(blastOutput):
						blast_size = os.path.getsize(blastOutput) / (1024 * 1024 * 1024)  # Convert to GB
						if blast_size > 1.0:
							print(f"blastOutput size is {blast_size:.2f} GB, larger than 1GB, calling blast_filter for tertiary filtering")
							filteredOutput = self.blast_filter(blastOutput, kmer_list)
							if os.path.exists(filteredOutput) and os.path.getsize(filteredOutput) > 0:
								print(f"Using blast_filter filtered result: {filteredOutput}")
								return filteredOutput
						else:
							print(f"blastOutput size is {blast_size:.2f} GB, smaller than 1GB, no further filtering needed")

					if os.path.exists(blastOutput) and os.path.getsize(blastOutput) > 0:
						print(f"Using blast filtered result: {blastOutput}")
						return blastOutput
				else:
					print(f"seqkitOutput size is {output_size:.2f} GB, smaller than 1GB, no further filtering needed")

			return seqkitOutput
		except Exception as e:
			print(f"Error occurred during k-mer filtering process: {e}")
			return None

	def minimap2(self):
		# Import required modules
		import os
		import subprocess
		import time

		alnname=self.roundInput.elongation.roundDir+"/potentialExtensionReads."+self.roundInput.elongation.base.tag+".bam"
		alnname1=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".fa"

		# If target file already exists, delete it first to avoid appending
		if os.path.exists(alnname):
			try:
				os.remove(alnname)
			except Exception as e:
				print(f"Warning: Unable to delete existing file {alnname}: {e}")

		# Add filtering layer
		print("Applying k-mer filtering...")
		# Use instance variables kmer_size and kmer_num
		filtered_reads = self.kmerfilter(self.roundInput.inputSeq, kmer_size=self.kmer_size, kmer_num=self.kmer_num)
		# Ensure filtering result is valid
		if filtered_reads and os.path.exists(filtered_reads) and os.path.getsize(filtered_reads) > 0:
			reads_to_use = os.path.abspath(filtered_reads)
			print(f"Using filtered reads file: {reads_to_use}, kmer_size={self.kmer_size}, kmer_num={self.kmer_num}")
		else:
			# Ensure using absolute path
			reads_to_use = os.path.abspath(self.roundInput.elongation.base.reads)
			print(f"k-mer filtering did not generate valid file, using original reads absolute path: {reads_to_use}")

		# Ensure inputSeq also uses absolute path
		input_seq_abs = os.path.abspath(self.roundInput.inputSeq)
		print(f"minimap2 using absolute paths: input_seq={input_seq_abs}, reads={reads_to_use}")
		
		commandline="minimap2 -t "+self.roundInput.elongation.base.thread+" -Y -ax asm20 "+input_seq_abs+" "+reads_to_use+" | samtools view -bS >"+alnname

		# If file already exists and is valid, return directly
		if os.path.exists(alnname)==True and os.path.getsize(alnname)!=0 and os.path.exists(alnname1)==True and os.path.getsize(alnname1)!=0:
			return alnname,commandline,str(0)
		else:
			# Use subprocess library instead of os.system, add timeout control
			try:
				print(f"Executing command: {commandline}")
				start_time = time.time()
				# Set timeout to 30 minutes
				result = subprocess.run(commandline, shell=True, timeout=1800, capture_output=True)

				# Check return code
				if result.returncode != 0:
					print(f"minimap2 command execution failed, return code: {result.returncode}")
					print(f"Error output: {result.stderr.decode('utf-8', errors='ignore')}")

					# Retry up to 2 times
					retries = 0
					while result.returncode != 0 and retries < 2:
						retries += 1
						print(f"Retry {retries} of minimap2 command...")
						result = subprocess.run(commandline, shell=True, timeout=1800, capture_output=True)

						if result.returncode == 0:
							break

						if retries == 2:
							print("minimap2 multiple retries failed, unable to execute correctly!")
							sys.exit(1)

				# Check output file size
				if os.path.exists(alnname):
					file_size = os.path.getsize(alnname)
					if file_size == 0:
						print(f"Warning: minimap2 generated BAM file size is 0 bytes")
					else:
						print(f"minimap2 generated BAM file size: {file_size} bytes")

				return alnname, commandline, str(result.returncode)

			except subprocess.TimeoutExpired:
				print("minimap2 command execution timeout, possibly due to processing large files")
				# If timeout, try to terminate process and clean up possible partial files
				if os.path.exists(alnname):
					os.remove(alnname)
				sys.exit(1)
			except Exception as e:
				print(f"Error occurred while executing minimap2 command: {e}")
				sys.exit(1)

	def readlog(self):
		logfilet=open(self.log,'r')
		for row in logfilet:
			row1=row.rstrip().split('\t')
			if row1[0]=='potentialExtensionReadsAln':
				self.potentialExtensionReadsAln=row1[1]
			elif row1[0]=='minimap2Command':
				self.minimap2Command=row1[1]
			elif row1[0]=='minimap2Output':
				self.minimap2Output=row1[1]
			elif row1[0]=='minimumThresholdReadsAln':
				self.minimumThresholdReadsAln=row1[1]
			elif row1[0]=='minimumThresholdReadsID':
				self.minimumThresholdReadsID=row1[1].split(';')
			elif row1[0]=='minimumThresholdExtensionReadsAln':
				self.minimumThresholdExtensionReadsAln=row1[1]
			elif row1[0]=='minimumThresholdExtensionReads':
				self.minimumThresholdExtensionReads=row1[1]
			elif row1[0]=='minimumThresholdExtensionReadsID':
				self.minimumThresholdExtensionReadsID=row1[1].split(';')
			elif row1[0]=='selectPotentialExtensionReadsAln':
				self.selectPotentialExtensionReadsAln=row1[1]

			elif row1[0]=='selectPotentialExtensionReadsID':
				if row1[1]!='None':
					self.selectPotentialExtensionReadsID=row1[1].split(';')
				else:
					self.selectPotentialExtensionReadsID=[]
			elif row1[0]=='selectReadsNum':
				self.selectReadsNum=int(row1[1])
			elif row1[0]=='extensionReadsAln':
				self.extensionReadsAln=row1[1]
			elif row1[0]=='extensionReads':
				self.extensionReads=row1[1]
			elif row1[0]=='extensionReadsID':
				if row1[1]!='None':
					self.extensionReadsID=row1[1].split(';')
				else:
					self.extensionReadsID=[]
			elif row1[0]=='extensionReadsNum':
				self.extensionReadsNum=int(row1[1])
			elif row1[0]=='selectMappingQuality':
				self.selectMappingQuality=int(row1[1])
			elif row1[0]=='selectAlignmentLength':
				self.selectAlignmentLength=int(row1[1])
			elif row1[0]=='selectNMAlignmentLengthratio':
				self.selectNMAlignmentLengthratio=float(row1[1])
			elif row1[0]=='readsExtensionLength':
				self.readsExtensionLength=int(row1[1])
			elif row1[0]=='extensionReadsEdge':
				self.extensionReadsEdge=int(row1[1])
			elif row1[0]=='note':
				if len(row1)==1:
					self.note=''
				else:
					self.note=row1[1]
		logfilet.close()

	def blast_filter(self, blast_result, kmer_list):
		"""
		Further filter blast results based on kmer frequency

		Parameters:
			blast_result: FASTA file path after blast filtering
			kmer_list: File path containing k-mers

		Returns:
			Further filtered FASTA file path
		"""
		# Import required modules
		import os
		import subprocess
		from Bio import SeqIO
		import collections
		import time

		# Ensure using absolute paths
		blast_result = os.path.abspath(blast_result)
		kmer_list = os.path.abspath(kmer_list)

		print(f"blast_filter function starting processing: blast_result={blast_result}, kmer_list={kmer_list}")

		# Create output directory (if it doesn't exist)
		output_dir = self.roundInput.elongation.roundDir
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)

		# Set output file path
		output_file = os.path.join(output_dir, "blast_filtered.fa")

		# Check blast_result file size
		blast_size = os.path.getsize(blast_result) / (1024 * 1024 * 1024)  # Convert to GB
		print(f"Blast result file size: {blast_size:.2f} GB")

		# If file is smaller than 1GB, return original file directly
		if blast_size < 1.0:
			print(f"Blast result file is smaller than 1GB, no further filtering needed")
			return blast_result
		
		# Read kmer file
		kmers = []
		with open(kmer_list, 'r') as f:
			for line in f:
				kmer = line.strip()
				if kmer:
					kmers.append(kmer)

		kmer_num = len(kmers)
		print(f"Read {kmer_num} kmers")

		# Count occurrences of each kmer in blast result file
		kmer_counts = collections.defaultdict(int)
		print(f"Starting to count occurrences of each kmer in blast results...")

		start_time = time.time()
		for kmer in kmers:
			try:
				cmd = f"seqkit locate -p {kmer} {blast_result} -m 0"
				result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

				# Count result lines (minus header line)
				lines = result.stdout.strip().split('\n')
				count = len(lines) - 1 if lines and lines[0] else 0
				kmer_counts[kmer] = count

				print(f"kmer: {kmer} occurrence count: {count}")

			except subprocess.CalledProcessError as e:
				print(f"Error finding kmer {kmer}: {e}")
				kmer_counts[kmer] = 0

		elapsed_time = time.time() - start_time
		print(f"Statistics completed, time elapsed: {elapsed_time:.2f} seconds")
		
		# Dynamic threshold adjustment logic
		thresholds = [50, 100, 150, 200]
		selected_kmers = []

		for threshold in thresholds:
			# Find kmers with occurrence count below threshold
			low_freq_kmers = [kmer for kmer, count in kmer_counts.items() if count <= threshold and count > 0]
			print(f"Number of kmers with occurrence count ≤ {threshold}: {len(low_freq_kmers)}")

			# Check if quantity requirement is met (at least half of total kmers)
			if len(low_freq_kmers) >= 0.5 * kmer_num:
				print(f"Found sufficient low-frequency kmers (threshold={threshold}), will use these kmers for filtering")
				selected_kmers = low_freq_kmers
				break

		# If all thresholds don't meet conditions, select relatively low-frequency kmers
		if not selected_kmers:
			print(f"All preset thresholds don't meet conditions, will select relatively low-frequency {int(0.5 * kmer_num)} kmers")

			# Sort all kmers by occurrence frequency
			sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1])

			# Select half of the kmers with lowest occurrence count
			selected_kmers = [kmer for kmer, _ in sorted_kmers[:int(0.5 * kmer_num)]]
			print(f"Selected {len(selected_kmers)} relatively low-frequency kmers")

		if not selected_kmers:
			print("No suitable low-frequency kmers found, returning original blast result")
			return blast_result

		# Create temporary file to save selected kmers
		selected_kmer_file = os.path.join(output_dir, "selected_kmers.txt")
		with open(selected_kmer_file, 'w') as f:
			for kmer in selected_kmers:
				f.write(f"{kmer}\n")

		print(f"Will use {len(selected_kmers)} kmers for final filtering")

		# Use seqkit grep to extract reads containing selected kmers
		try:
			cmd = f"seqkit grep -f {selected_kmer_file} -s {blast_result} -o {output_file}"
			print(f"Executing command: {cmd}")
			subprocess.run(cmd, shell=True, check=True)

			# Check output file
			if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
				output_size = os.path.getsize(output_file) / (1024*1024)  # MB
				print(f"Final filtering completed, output file size: {output_size:.2f} MB")
				return str(output_file)
			else:
				print(f"No valid results after filtering, returning original blast result")
				return blast_result

		except Exception as e:
			print(f"Error during final filtering process: {e}")
			return blast_result

		# Clean up temporary files
		try:
			os.remove(selected_kmer_file)
		except:
			pass
