#!/usr/bin/env python3
# -*-coding:utf-8 -*-

#**********************************************************************************
#filename:GAP_Filler.py
#*********************************************************************************

import re
import os
import sys
import getopt
import pysam
from pysam import AlignmentFile
import Bio
from Bio import SeqIO
import GapFillerClass
from GapFillerClass import GapFillerClass,mummer

class GapFiller(object):
	def __init__(self,parameterlist, kparameters):
		self.mode,self.remove,self.thread,self.reads,self.out,self.seqleft,self.seqright,self.flag,self.edge,self.filterDepth,self.MaximunExtensionLength,self.readsDict,self.maxReadsLen,self.seedLen=parameterlist
		self.kmer_size, self.kmer_num, self.kmer_length, self.j = kparameters
		self.resume_round = None  # Default: do not resume from checkpoint
		
		out=self.out
		self.log=out+"/process.log"
		if not os.path.exists(self.log):
			open(self.log,'w').close()
		self.summary=out+"/process.summary"
		if not os.path.exists(self.summary):
			open(self.summary,'w').close()
		self.agp=out+"/process.agp"
		self.usedReads=[]
		if not os.path.exists(self.agp):
			open(self.agp,'w').close()
		self.name=out.split('/')[-1]
		self.outfile=out+"/process"
		if not os.path.exists(self.outfile):
			os.makedirs(self.outfile)

		# Check if checkpoint file exists
		self.checkpoint_file = out+"/checkpoint.info"
		if os.path.exists(self.checkpoint_file):
			self.load_checkpoint()

		if self.flag=='left':
			self.initialSeq=self.seqleft
			self.terminalSeq=self.seqright
			self.tag='left'
		else:
			self.initialSeq=self.seqright
			self.terminalSeq=self.seqleft
			self.tag='right'

		# Add seed sequence length check logic
		self.check_and_adjust_seed_sequences()

		self.Elongation=Elongation(self)
	
	def load_checkpoint(self):
		"""Read checkpoint information and decide which round to continue from"""
		try:
			with open(self.checkpoint_file, 'r') as f:
				checkpoint_data = f.read().strip()
				if checkpoint_data and checkpoint_data.startswith("round:"):
					self.resume_round = int(checkpoint_data.split("round:")[1])
					print(f"Found checkpoint information, will continue from round{self.resume_round}")
		except Exception as e:
			print(f"Failed to read checkpoint information: {e}")
			self.resume_round = None

	def save_checkpoint(self, round_num):
		"""Save current execution round"""
		try:
			with open(self.checkpoint_file, 'w') as f:
				f.write(f"round:{round_num}")
		except Exception as e:
			print(f"Failed to save checkpoint information: {e}")
	
	def check_and_adjust_seed_sequences(self):
		"""Check seed sequence length, use entire sequence if shorter than seedLen"""
		try:
			# Check initial seed sequences
			for seq_file in [self.initialSeq, self.terminalSeq]:
				if not os.path.exists(seq_file):
					continue

				# Read sequence file
				original_seqs = []
				for record in SeqIO.parse(seq_file, "fasta"):
					original_seqs.append((record.id, record.seq))

				# Skip if no sequences
				if not original_seqs:
					continue

				# Check each sequence length
				need_update = False
				for idx, (seq_id, seq) in enumerate(original_seqs):
					if len(seq) < self.seedLen:
						print(f"Warning: Sequence {seq_id} length ({len(seq)}bp) is shorter than set seed length ({self.seedLen}bp), will use entire sequence.")
						need_update = True

				# Continue to next file if no update needed
				if not need_update:
					continue

				# Create temporary file name
				temp_file = seq_file + ".temp"

				# Write new file
				with open(temp_file, 'w') as f:
					for seq_id, seq in original_seqs:
						f.write(f">{seq_id}\n{seq}\n")

				# Replace original file
				os.replace(temp_file, seq_file)
				print(f"Updated sequence file: {seq_file}")

		except Exception as e:
			print(f"Error checking and adjusting seed sequences: {e}")
			import traceback
			traceback.print_exc()

class Elongation(object):
	def __init__(self,base):
		self.base=base
		self.roundNum=1
		self.usedReads=[]
		self.extensionLen=0
		self.endSignal=False
		self.out=base.out
		self.resume_data = {}  # Store data needed for recovery
		self.has_direct_connection = False  # New: flag for direct connection existence
		self.direct_connection_result = None  # New: save direct connection result

		logfile=open(self.base.log,'w')
		summaryfile=open(self.base.summary,'w')
		
		# Check if seqleft and seqright can be directly connected before starting extension
		self.check_direct_connection(logfile)

		# If need to continue from specific round
		if self.base.resume_round is not None:
			self.resume_from_checkpoint(self.base.resume_round)
		
		while self.endSignal==False:
			self.ElongationInit(logfile)
			
			# Save checkpoint information for current round
			self.base.save_checkpoint(self.roundNum)
			
			try:
				self.lastRoundUsedReads=self.ElongateSeq(logfile,summaryfile,self.lastRoundUsedReads,self.extensionLen)
				if self.roundResult.ExtensionReads.note=='' and 'No extension contigs or reads found' not in self.roundResult.ExtensionContigs.selectContigNote:
					self.extensionLen=self.extensionLen+self.roundResult.ExtensionContigs.extensionLength
				if self.roundNum==1:
					if self.roundResult.ExtensionReads.note=='' and 'No extension contigs or reads found' not in self.roundResult.ExtensionContigs.selectContigNote:
						atgInitial=self.roundResult.ExtensionContigs.selectExtensionContigsAln[0][0]
					else:
						atgInitial=''
			except Exception as e:
				# Record checkpoint when exception occurs
				print(f"Round {self.roundNum} execution error: {str(e)}")
				self.base.save_checkpoint(self.roundNum)
				raise  # Re-raise exception for external handling
			
			self.roundNum+=1
			self.removeFile()
		summaryfile.close()

		# Delete checkpoint file after completing all rounds
		if os.path.exists(self.base.checkpoint_file):
			os.remove(self.base.checkpoint_file)

		# Check if direct connection result should be used before generating final result
		use_direct_connection = False

		# Check if direct connection result exists and should be used
		if self.has_direct_connection and self.direct_connection_result:
			logfile.write("\n========== Check if direct connection result should be used ==========\n")
			
			# Check if extension result can close gap
			extension_can_close_gap = False
			if hasattr(self, 'roundResult') and hasattr(self.roundResult, 'roundOutput'):
				if os.path.exists(self.roundDir+"/linkedSequence.fasta") and os.path.getsize(self.roundDir+"/linkedSequence.fasta") != 0:
					extension_can_close_gap = True

					# Further check if gap length is negative
					process_log_negative = False
					try:
						if os.path.exists(self.base.log):
							with open(self.base.log, 'r') as f:
								log_content = f.read()
								# Search for latest GAP Length record in log
								gap_length_matches = re.findall(r'GAP Length: ([-\d]+)', log_content)
								if gap_length_matches:
									last_gap_length = int(gap_length_matches[-1])
									if last_gap_length < 0:
										process_log_negative = True
										logfile.write(f"Detected negative GAP Length in extension phase: {last_gap_length}\n")
					except Exception as e:
						logfile.write(f"Error checking process.log: {str(e)}\n")
			
			# Decide whether to use direct connection result
			if not extension_can_close_gap or process_log_negative:
				use_direct_connection = True
				logfile.write("Decided to use direct connection result because:\n")
				if not extension_can_close_gap:
					logfile.write("- Extension phase failed to close gap\n")
				if process_log_negative:
					logfile.write("- Extension phase gap length is negative\n")

				# Use direct connection result
				logfile.write(f"Using direct connection result file: {self.direct_connection_result['final_seq_file']}\n")
				self.finalSeq = self.direct_connection_result['final_seq_file']

				# Create marker file indicating direct connection result was used
				with open(os.path.join(self.out, "USED_DIRECT_CONNECTION"), 'w') as f:
					f.write(f"Used direct connection result instead of extension result\n")
					f.write(f"Direct connection file: {self.direct_connection_result['final_seq_file']}\n")
					f.write(f"GAP Length: {self.direct_connection_result['gap_length']}\n")
					f.write(f"Overlap length: {self.direct_connection_result['overlap_length']}\n")
					f.write(f"Identity: {self.direct_connection_result['identity']}%\n")
			else:
				logfile.write("Extension phase successfully closed gap and gap length is not negative, using extension result\n")
			
			logfile.write("============================\n\n")
		
		self.roundNum=self.roundNum-1
		
		# If not using direct connection result, generate final result in original way
		if not use_direct_connection:
			self.finalSeq=self.base.out+'/'+self.base.name+'.final.fa'
			if os.path.exists(self.roundDir+"/linkedSequence.fasta")==True and  os.path.getsize(self.roundDir+"/linkedSequence.fasta")!=0:
				fileofs=open(self.finalSeq,'w')
				l='>'+self.base.name+"\n"
				fileofs.writelines(l)
				for gseq in SeqIO.parse(self.roundResult.roundOutput.linkedSequence,'fasta'):
					l0=gseq.description
					l1=l0.split('\t')
					l=gseq.seq+'\n'
					fileofs.writelines(l)
					for l11 in l1:
						l2=l11.split(':')
						if l2[0]=="Aln":
							atgTerminal=l2[1]
				fileofs.close()
			elif self.roundResult.ExtensionReads.note!='' or 'No extension contigs or reads found' in self.roundResult.ExtensionContigs.selectContigNote or "Reach the maximum Length" in self.roundResult.ExtensionContigs.selectContigNote:
				fileofs=open(self.finalSeq,'w')
				if self.roundResult.ExtensionReads.note=='':
					if 'Reach the maximum Length' in self.roundResult.ExtensionContigs.selectContigNote:
						l='>'+self.base.name+"_reachMaximumLength\n"
					else:
						l='>'+self.base.name+"_noExtensionContigsorReads\n"
				else:
					l='>'+self.base.name+"_noExtensionContigsorReads\n"
				fileofs.writelines(l)
				atgTerminal=''
				for gseq in SeqIO.parse(self.roundInputSeq,'fasta'):
					l=gseq.seq+'\n'
					fileofs.writelines(l)
				fileofs.close()
			else:
				fileofs=open(self.finalSeq,'w')
				l='>'+self.base.name+"_noNewExtensionReads\n"
				fileofs.writelines(l)
				atgTerminal=''
				for gseq in SeqIO.parse(self.roundResult.roundOutput.outputSequence,'fasta'):
					l=gseq.seq+'\n'
					fileofs.writelines(l)
				fileofs.close()
		
		# If --remove=1 is set, ensure final result file won't be deleted
		if self.base.remove == 1:
			# Ensure final result file has been generated and has content
			if os.path.exists(self.finalSeq) and os.path.getsize(self.finalSeq) > 0:
				print(f"Generated final result file: {self.finalSeq}")

				# Copy important log files to safe location
				import shutil
				log_backup = self.base.out + '/process.log.backup'
				summary_backup = self.base.out + '/process.summary.backup'
				agp_backup = self.base.out + '/process.agp.backup'

				if os.path.exists(self.base.log):
					shutil.copy2(self.base.log, log_backup)
				if os.path.exists(self.base.summary):
					shutil.copy2(self.base.summary, summary_backup)
				if os.path.exists(self.base.agp):
					shutil.copy2(self.base.agp, agp_backup)

				print("Backed up important log files")
	
		for fgseq in SeqIO.parse(self.finalSeq,'fasta'):
			l='Final ExtensionSequence: '+fgseq.id+"\n"
			l=l+'Final EXtendFile: '+self.finalSeq+"\n"
			logfile.writelines(l)
		logfile.close()
		
		agpfile=open(self.base.agp,'w')
		for gseq in SeqIO.parse(self.base.initialSeq,'fasta'):
			initialSequence=gseq
		if atgInitial!='':
			if self.base.flag=='left':
				atgInitial1=atgInitial.split('\t')
				st=len(initialSequence.seq)-(int(atgInitial1[7])-int(atgInitial1[1]))
				if atgTerminal=='':
					l=self.base.name+"\t1"+"\t"+str(st)+"\t1\tw\t"+initialSequence.id+"\t1\t"+str(st)+"\t+\n"
					agpfile.writelines(l)
					l=self.base.name+"\t"+str(st+1)+"\t"+str(len(fgseq.seq))+"\t2\tw\t"+fgseq.id+"\t"+str(st+1)+"\t"+str(len(fgseq.seq))+"\t+\n"
					agpfile.writelines(l)
					self.TerminalSeq=None
				else:
					l=self.base.name+"\t1"+"\t"+str(st)+"\t1\tw\t"+initialSequence.id+"\t1\t"+str(st)+"\t+\n"
					agpfile.writelines(l)
					atgTerminal1=atgTerminal.split(';')
					eT=len(fgseq.seq)-(int(atgTerminal1[7])-int(atgTerminal1[4]))
					l=self.base.name+"\t"+str(st+1)+"\t"+str(eT)+"\t2\tw\t"+fgseq.id+"\t"+str(st+1)+"\t"+str(eT)+"\t+\n"
					agpfile.writelines(l)
					Terminalname1=atgTerminal1[-2]
					l=self.base.name+"\t"+str(int(eT)+1)+"\t"+str(len(fgseq.seq))+"\t3\tw\t"+Terminalname1+"\t"+str(int(atgTerminal1[1])+1)+"\t"+atgTerminal1[7]+"\t+\n"
					agpfile.writelines(l)
					self.TerminalSeq=Terminalname1
		
			else:
				atgInitial1=atgInitial.split('\t')
				st=len(fgseq.seq)-len(initialSequence.seq)+int(atgInitial1[0])-1
				if atgTerminal=='':
					atgInitial1=atgInitial.split('\t')
					l=self.base.name+"\t1\t"+str(st)+"\t1\tw\t"+fgseq.id+"\t1\t"+str(st)+"\t+\n"
					agpfile.writelines(l)
					l=self.base.name+"\t"+str(st+1)+"\t"+str(len(fgseq))+"\t2\tw\t"+initialSequence.id+"\t"+atgInitial1[0]+"\t"+str(len(initialSequence.seq))+"\t+\n"
					agpfile.writelines(l)
					self.TerminalSeq=None
					
				else:
					atgTerminal1=atgTerminal.split(';')
					atgInitial1=atgInitial.split('\t')
					sT=int(atgTerminal1[1])
					l=self.base.name+"\t1\t"+str(sT)+"\t1\tw\t"+atgTerminal1[-2]+"\t1\t"+atgTerminal1[1]+"\t+\n"
					agpfile.writelines(l)
					l=self.base.name+"\t"+str(sT+1)+"\t"+str(st)+"\t2\tw\t"+fgseq.id+"\t"+str(sT+1)+"\t"+str(st)+"\t+\n"
					agpfile.writelines(l)
					l=self.base.name+"\t"+str(st+1)+"\t"+str(len(fgseq))+"\t3\tw\t"+initialSequence.id+"\t1\t"+atgTerminal1[1]+"\t+\n"
					agpfile.writelines(l)
					self.TerminalSeq=atgTerminal1[-2]
		else:
			l=self.base.name+"\t1"+"\t"+str(len(initialSequence.seq))+"\t1\tw\t"+initialSequence.id+"\t1\t"+str(len(initialSequence.seq))+"\t+\n"
			agpfile.writelines(l)
		agpfile.close()
		self.removeFile()
				
	def removeFile(self):
		if os.path.exists(self.roundDir+"/hifiasm")==True or os.path.exists(self.roundDir+"/potentialExtensionReads.left.sam")==True:
			if self.base.remove==2 or self.base.remove==1:
				commondline='rm '+self.roundDir+"/*.bam"
				os.system(commondline)
				commondline='rm -rf '+self.roundDir+"/hifiasm"
				os.system(commondline)
		if self.endSignal==True:
			if os.path.exists(self.base.outfile)==True:
				if self.base.remove==1:
					# Modified: Don't delete entire output directory immediately, clean up after all processing is complete
					# First save final result file
					final_seq = self.base.out+'/'+self.base.name+'.final.fa'
					if os.path.exists(final_seq):
						# Ensure final result file has been generated
						print(f"Keeping final result file: {final_seq}")
					else:
						print(f"Warning: Final result file not yet generated, skipping output directory cleanup")
						return

					# Only delete intermediate process files, keep final results
					for root, dirs, files in os.walk(self.base.outfile):
						for file in files:
							file_path = os.path.join(root, file)
							os.remove(file_path)

					# Don't delete directory structure, only delete files
					print(f"Cleaned up intermediate process files, kept final results and directory structure")

	def ElongateSeq(self,logfile,summaryfile,lastRoundUsedReads,extensionLen):
		self.lastRoundUsedReads=lastRoundUsedReads
		roundLog=open(self.roundLog,'w')
		roundSummary=open(self.roundSummary,'w')
		self.roundResult=GapFillerClass(self,out=self.out)
		print ("MaximunExtensionLength",self.base.MaximunExtensionLength,"TotalExtensionLength",extensionLen)
		#sys.exit()
		if self.base.MaximunExtensionLength!=None:
			if extensionLen>self.base.MaximunExtensionLength: #Max
				#print (self.roundResult.ExtensionContigs.extensionLength,self.base.MaximunExtensionLength,"MaximunExtensionLength")
				self.roundResult.ExtensionContigs.selectContigNote=self.roundResult.ExtensionContigs.selectContigNote+"Reach the maximum Length\n"
				print (self.roundResult.ExtensionContigs.selectContigNote)
				#sys.exit()
		logLine,summeryLine=self.writelog(extensionLen)
		print (logLine)
		roundLog.writelines(logLine)
		roundSummary.writelines(summeryLine)
		logfile.writelines(logLine)
		summaryfile.writelines(summeryLine)
		roundLog.close()
		roundSummary.close()
		if self.roundResult.ExtensionReads.note=='' and 'No extension contigs or reads found' not in self.roundResult.ExtensionContigs.selectContigNote and "Reach the maximum Length" not in self.roundResult.ExtensionContigs.selectContigNote:
			return self.roundResult.roundOutput.ExtensionUsedReads
		else:
			return []
	
	def writelog(self,extensionLen):
		logLine='\n\n*****************\n\n'
		logLine+='\toutputPath: '+str(self.roundDir)+"\n"
		logLine+="\tseedSequenceLength: "+str(self.base.seedLen)+"\n"
		logLine+="\tinitialSequenceFile: "+str(self.base.initialSeq)+"\n"
		for gseq in SeqIO.parse(self.roundInputSeq,'fasta'):
			inputSeq=gseq
		logLine+="\t\tinitialSeqnenceID: "+str(inputSeq.id)+"\n\t\tinitialSeqnenceLength: "+str(len(inputSeq.seq))+"\n"
		logLine+="\tterminalSequenceFile: "+str(self.base.terminalSeq)+"\n"
		logLine+="\tseedSequenceFile: "+str(self.roundResult.roundInput.inputSeq)+"\n"
		logLine+="\t\tseedSeqnenceID: "+str(self.roundResult.roundInput.inputSeedSequence.id)+"\n\t\tseedSeqnenceLength: "+str(len(self.roundResult.roundInput.inputSeedSequence.seq))+"\n\n"
		if self.roundResult.ExtensionReads.note=='' and 'No extension contigs or reads found' not in self.roundResult.ExtensionContigs.selectContigNote and "Reach the maximum Length" not in self.roundResult.ExtensionContigs.selectContigNote:
			logLine+='minimap2Commond: '+str(self.roundResult.ExtensionReads.minimap2Command)+"\n"
			logLine+='\textensionReads: \n\tselectReadsNum: '+str(self.roundResult.ExtensionReads.selectReadsNum)+"\n"
			logLine+="\t\tselectReadsAln: "+str(self.roundResult.ExtensionReads.selectPotentialExtensionReadsAln)+"\n"
			logLine+="\t\tselectMappingQuality: "+str(self.roundResult.ExtensionReads.selectMappingQuality)+"\n"
			logLine+="\t\tselectAlignmentLength: "+str(self.roundResult.ExtensionReads.selectAlignmentLength)+"\n"
			logLine+="\t\tselectNMAlignmentLengthratio: "+str(self.roundResult.ExtensionReads.selectNMAlignmentLengthratio)+"\n\n"
			logLine+="\textensionReadsNum: "+str(self.roundResult.ExtensionReads.extensionReadsNum)+"\n"
			logLine+="\t\textensionReadsAln: "+str(self.roundResult.ExtensionReads.extensionReadsAln)+"\n"
			logLine+="\t\textensionReadsFile: "+str(self.roundResult.ExtensionReads.extensionReads)+"\n"
			logLine+="\t\textensionReadsMinimumExtensionLength: "+str(self.roundResult.ExtensionReads.readsExtensionLength)+"\n"
			logLine+="\t\textensionReadsMaximumEdge: "+str(self.roundResult.ExtensionReads.extensionReadsEdge)+"\n\n"

			logLine+="\textensionSequnece: "+str(self.roundResult.ExtensionContigs.extensionContigs)+"\n"
			logLine+="\textensionSequneceNote: "+str(self.roundResult.ExtensionContigs.extensionSeqNote)+"\n"+str(self.roundResult.ExtensionContigs.selectContigNote)+"\n"
			logLine+="\t\textensionSequneceIdentity: "+str(self.roundResult.ExtensionContigs.selectContigIdentity)+"\n"
			logLine+="\t\textensionSequneceMinimumExtensionLength: "+str(self.roundResult.ExtensionContigs.selectContigAlnLength)+"\n"
			logLine+="\t\textensionSequneceMaximumEdge: "+str(self.roundResult.ExtensionContigs.selectContigDistance)+"\n"
			logLine+="\t\textensionSequneceAlnMerge: "+str(self.roundResult.ExtensionContigs.contigAlnMerge)+"\n"
			logLine+="\t\textensionSequneceAlnMergeIdentity: "+str(self.roundResult.ExtensionContigs.contigAlnMergeIdentity)+"\n"
			logLine+="\textensionSeedSequenceFile: "+str(self.roundResult.ExtensionContigs.extensionSequence)+"\n"
			logLine+="\textensionLength: "+str(self.roundResult.ExtensionContigs.extensionLength)+"\n"
			logLine+="\t\textensionContigOrReadsID:\n\t\t\t"+'\n\t\t\t'.join(self.roundResult.ExtensionContigs.extensionContigID)+"\n"
			newReads,note=self.updateUsedReads()
			logLine+=note
			
			# Check if self.roundResult has roundOutput attribute
			if hasattr(self.roundResult, 'roundOutput'):
				logLine+="\t\tusedReadsNum: "+str(len(self.roundResult.roundOutput.ExtensionUsedReads))+"\n"
				logLine+="\t\tusedReads:\n\t\t\t"+"\n\t\t\t".join(self.roundResult.roundOutput.ExtensionUsedReads)+"\n"
				logLine+="\t\tusedNewReads:\n\t\t\t"+"\n\t\t\t".join(newReads)+"\n\n"
				logLine+="\toutputFile: "+str(self.roundResult.roundOutput.outputSequence)+"\n"
				logLine+="\t\toutputSequenceLength: "+str(self.roundResult.roundOutput.totalOutputSequenceLength)+"\n\n"
				logLine+="\t\ttotalExtensionLength: "+str(extensionLen)+"\n\n"
				
				if self.roundResult.roundOutput.linkedSequenceNote!='':
					logLine+='\tGAP can be closed!\n'+str(self.roundResult.roundOutput.linkedSequenceNote)+"\nLinkedSequence File: "+str(self.roundResult.roundOutput.linkedSequence)+"\n"
					logLine+='Endloop!\t'+str(self.roundResult.roundOutput.linkedSequence)+"\n"
					for gseq in SeqIO.parse(self.roundResult.roundOutput.linkedSequence,'fasta'):
						l0=gseq.description
						l1=l0.split('\t')
						for l11 in l1:
							l2=l11.split(':')
							if l2[0]=="Aln":
								closectg=l2[1].split(';')
								logLine+="Linked ctg:\t"+closectg[-2]+"\n"
								logLine+='\t'.join(closectg)+"\n"
							self.endSignal=True
				else:
					logLine+='\tGAP still not closed!\n'
			else:
				# If no roundOutput attribute, add appropriate information
				logLine+="\t\tNo output sequence generated in this round.\n"
				logLine+="\t\ttotalExtensionLength: "+str(extensionLen)+"\n\n"
				logLine+='\tGAP could not be processed further - missing output information.\n'
				self.endSignal=True
			
			logLine+='\n\n*****************\n\n'
			summeryLine='round'+str(self.roundNum)+"\t"+str(len(inputSeq.seq))+"\t"
			
			# Also check roundOutput attribute to correctly generate summary line
			if hasattr(self.roundResult, 'roundOutput'):
				summeryLine+=str(self.roundResult.roundOutput.totalOutputSequenceLength)+"\t"+str(self.roundResult.ExtensionContigs.extensionLength)+"\t"+"-ovl-".join(self.roundResult.ExtensionContigs.extensionContigID)+"\t"+str(len(newReads))+"\t"+";".join(newReads)+"\t"+str(len(self.roundResult.roundOutput.ExtensionUsedReads))+"\t"+";".join(self.roundResult.roundOutput.ExtensionUsedReads)+"\n"
			else:
				summeryLine+="0\t0\tno_extension\t0\tnone\t0\tnone\n"
			
			return logLine,summeryLine
		else:
			if self.roundResult.ExtensionReads.note!='':
				logLine+='No ExtensionReads or ExtensionContig Found\n'
				logLine+="\n\tPossible reason 1: Your original data can only support assembly up to the current length. Assembly results can be found at: "+self.base.out+'/'+self.base.name+'.final.fa'+"\n\n"
				logLine+="\tPossible reason 2: There may be an issue with the kmer filtering parameters. The current default parameters for kmer filtering are: --kmer_size 41, --kmer_num 4. "
				logLine+="This parameter combination is based on rice HiFi sequencing data. "
				logLine+="If you are using it to assemble other species' genomes, you can set appropriate parameters through --kmer_size and --kmer_num.\n\n"
			else:
				if "Reach the maximum Length" not in self.roundResult.ExtensionContigs.selectContigNote:
					logLine+='No ExtensionReads or ExtensionContig Found\n'
				else:
					logLine+='Reach the maximum Length\n'	
			logLine+="Endloop!\t"+self.roundResult.ExtensionReads.note+"\n"
			logLine+="\t\ttotalExtensionLength: "+str(extensionLen)+"\n\n"
			self.endSignal=True
			summeryLine=''
			return logLine,summeryLine

	def updateUsedReads(self):
		newReads=[]
		sameWithLastRound=[]
		note=''
		for i1 in self.roundResult.roundOutput.ExtensionUsedReads:
			if i1 not in self.usedReads:
				self.usedReads.append(i1)
				newReads.append(i1)
			if i1 in self.lastRoundUsedReads:
				sameWithLastRound.append(i1)
		if len(newReads)==0 :
			note='No New ExtensionReads Found,\t'
			if len(sameWithLastRound)==0:
				note+='Not same ExtensionReads with last round ExtensionReads, end up a loop!!\n'
				self.endSignal=True
			else:
				note+='Same ExtensionReads with the last round,continune the loop!\n'
				if self.roundResult.ExtensionContigs.extensionLength==0:
					note+='However,ExtensionLength==0,end up a loop!!\n'
					self.endSignal=True
		return newReads,note

	
	def ElongationInit(self,logfile):
		l='\n\n*****************\n\nExtensionRound '+str(self.roundNum)+'\n'
		logfile.writelines(l)
		print (l)
		self.roundDir=self.base.outfile+'/round'+str(self.roundNum)
		if not os.path.exists(self.roundDir):
			os.makedirs(self.roundDir)
		self.lastRoundDir=self.base.outfile+'/round'+str(self.roundNum-1)
		self.roundLog=self.base.outfile+'/round'+str(self.roundNum)+"/log"
		self.roundSummary=self.base.outfile+'/round'+str(self.roundNum)+"/summary"

		if self.roundNum!=1:
			self.roundInputSeq=self.lastRoundDir+'/outputExtensionSequence.fasta'
			# Check updated seed sequence
			self.check_seed_sequence_length()
		else:
			self.roundInputSeq=self.base.initialSeq
			self.lastRoundUsedReads=[]

	def check_seed_sequence_length(self):
		"""Check seed sequence length in each iteration, use entire sequence if shorter than seedLen"""
		try:
			if not os.path.exists(self.roundInputSeq):
				print(f"Warning: Cannot find seed sequence file {self.roundInputSeq}")
				return

			# Read sequence file
			original_seqs = []
			for record in SeqIO.parse(self.roundInputSeq, "fasta"):
				original_seqs.append((record.id, record.seq))

			# Return if no sequences
			if not original_seqs:
				print(f"Warning: No sequences in seed sequence file {self.roundInputSeq}")
				return

			# Check each sequence length
			need_update = False
			for idx, (seq_id, seq) in enumerate(original_seqs):
				if len(seq) < self.base.seedLen:
					print(f"Round {self.roundNum}: Sequence {seq_id} length ({len(seq)}bp) is shorter than set seed length ({self.base.seedLen}bp), will use entire sequence.")
					need_update = True

			# Return directly if no update needed
			if not need_update:
				return

			# Create temporary file name
			temp_file = self.roundInputSeq + ".temp"

			# Write new file
			with open(temp_file, 'w') as f:
				for seq_id, seq in original_seqs:
					f.write(f">{seq_id}\n{seq}\n")

			# Replace original file
			os.replace(temp_file, self.roundInputSeq)
			print(f"Updated seed sequence file: {self.roundInputSeq}")

		except Exception as e:
			print(f"Error checking seed sequence length: {e}")
			import traceback
			traceback.print_exc()

	def resume_from_checkpoint(self, target_round):
		"""Resume execution state from specified round"""
		print(f"Preparing to resume execution from round{target_round}...")

		# Set initial round to target round
		self.roundNum = target_round

		# Determine previous round directory
		last_round = target_round - 1
		if last_round > 0:
			self.lastRoundDir = self.base.outfile+'/round'+str(last_round)
			
			# Read previous round's output sequence as input for this round
			self.roundInputSeq = self.lastRoundDir+'/outputExtensionSequence.fasta'

			# If previous round's output sequence doesn't exist, use original sequence
			if not os.path.exists(self.roundInputSeq):
				print(f"Warning: Cannot find previous round's output sequence {self.roundInputSeq}, will use initial sequence")
				self.roundInputSeq = self.base.initialSeq
				self.lastRoundUsedReads = []
			else:
				# Try to recover used reads information from previous round's summary file
				summary_file = self.lastRoundDir+"/summary"
				self.lastRoundUsedReads = []
				self.usedReads = []

				if os.path.exists(summary_file):
					try:
						with open(summary_file, 'r') as f:
							last_line = f.readlines()[-1].strip()
							if last_line:
								fields = last_line.split('\t')
								if len(fields) >= 9:  # Ensure sufficient fields
									self.lastRoundUsedReads = fields[8].split(';') if fields[8] != 'none' else []
									self.usedReads = self.lastRoundUsedReads.copy()
					except Exception as e:
						print(f"Failed to recover reads information from summary file: {e}")
						self.lastRoundUsedReads = []
						
				# Try to read previous round's extension length
				logs_file = self.lastRoundDir+"/log"
				if os.path.exists(logs_file):
					try:
						with open(logs_file, 'r') as f:
							content = f.read()
							extension_matches = re.findall(r'totalExtensionLength:\s+(\d+)', content)
							if extension_matches:
								self.extensionLen = int(extension_matches[-1])
								print(f"Recovered current extension length from log file: {self.extensionLen}")
					except Exception as e:
						print(f"Failed to recover extension length from log file: {e}")
						# If unable to recover extension length, try to calculate
						try:
							initial_seq_len = 0
							for gseq in SeqIO.parse(self.base.initialSeq,'fasta'):
								initial_seq_len = len(gseq.seq)

							current_seq_len = 0
							for gseq in SeqIO.parse(self.roundInputSeq,'fasta'):
								current_seq_len = len(gseq.seq)

							self.extensionLen = current_seq_len - initial_seq_len
							if self.extensionLen < 0:
								self.extensionLen = 0
							print(f"Calculated extension length: {self.extensionLen}")
						except Exception as calc_e:
							print(f"Failed to calculate extension length: {calc_e}")
							self.extensionLen = 0
		else:
			# If first round, use initial settings
			self.roundInputSeq = self.base.initialSeq
			self.lastRoundUsedReads = []

		print(f"Ready to continue execution from round{target_round}, current extension length: {self.extensionLen}")

	def check_direct_connection(self, logfile):
		"""Check if seqleft and seqright can be directly connected (including tandem repeats and normal overlaps)"""
		import tempfile
		import os

		logfile.write("\n========== Direct Connection Check ==========\n")

		# Check if files exist
		if not os.path.exists(self.base.seqleft) or not os.path.exists(self.base.seqright):
			logfile.write("Cannot check direct connection: left or right sequence file does not exist\n")
			logfile.write("============================\n\n")
			return

		# Read left and right sequences
		left_seq = None
		right_seq = None

		for record in SeqIO.parse(self.base.seqleft, "fasta"):
			left_seq = record.seq
			left_id = record.id
			break

		for record in SeqIO.parse(self.base.seqright, "fasta"):
			right_seq = record.seq
			right_id = record.id
			break

		if left_seq is None or right_seq is None:
			logfile.write("Cannot check direct connection: left or right sequence is empty\n")
			logfile.write("============================\n\n")
			return

		logfile.write(f"Left sequence length: {len(left_seq)}bp\n")
		logfile.write(f"Right sequence length: {len(right_seq)}bp\n")
		
		# Create temporary directory for mummer analysis
		temp_dir = os.path.join(self.out, "direct_connection_check")
		os.makedirs(temp_dir, exist_ok=True)

		# Write sequences to temporary files
		left_file = os.path.join(temp_dir, "left.fasta")
		right_file = os.path.join(temp_dir, "right.fasta")

		with open(left_file, "w") as f:
			f.write(f">{left_id}\n{left_seq}\n")

		with open(right_file, "w") as f:
			f.write(f">{right_id}\n{right_seq}\n")

		# Use mummer for alignment
		mummer_prefix = os.path.join(temp_dir, "direct_connection")
		mummer_output = mummer(left_file, right_file, mummer_prefix)

		logfile.write(f"Running mummer alignment: {mummer_output}\n")
		
		# Analyze alignment results
		if os.path.exists(mummer_output) and os.path.getsize(mummer_output) > 0:
			# Collect all possible overlaps
			alignments = []
			# Collect simple overlaps (non-tandem repeats)
			simple_overlaps = []
			
			with open(mummer_output, "r") as f:
				for line in f:
					if line.startswith("#"):
						continue
					parts = line.strip().split("\t")
					if len(parts) < 10:
						continue
					
					# Parse alignment information
					try:
						s1, e1, s2, e2 = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
						length = int(parts[4])
						identity = float(parts[6])

						# Check various types of overlaps
						len_left = len(left_seq)
						len_right = len(right_seq)

						# Condition 1: Simple overlap - left sequence end matches right sequence start
						is_simple_overlap = (e1 > len_left * 0.8 and s2 < len_right * 0.2)

						# Condition 2: Tandem repeat - matching region also appears in middle part of left sequence
						is_tandem_repeat = (e1 > len_left * 0.5 and s2 < len_right * 0.5 and length > 20)

						# Collect by different types of overlaps
						if is_simple_overlap and length >= 10 and identity >= 80:
							simple_overlaps.append((s1, e1, s2, e2, length, identity, parts))
							logfile.write(f"Found simple overlap: left {s1}-{e1}, right {s2}-{e2}, length {length}bp, identity {identity}%\n")

						if is_tandem_repeat and length >= 20 and identity >= 80:
							alignments.append((s1, e1, s2, e2, length, identity, parts))
							logfile.write(f"Found tandem repeat match: left {s1}-{e1}, right {s2}-{e2}, length {length}bp, identity {identity}%\n")
					except ValueError:
						continue
			
			# Prioritize simple overlaps if they exist
			if simple_overlaps:
				# Sort by overlap length and identity, select the best one
				best_overlap = sorted(simple_overlaps, key=lambda x: (x[4], x[5]), reverse=True)[0]
				s1, e1, s2, e2, length, identity, parts = best_overlap

				# Calculate gap length (negative value indicates overlap)
				gap_length = -(length)  # Overlap length as negative gap length

				# Build final sequence (remove overlapping part)
				final_seq = left_seq + right_seq[e2:]

				# Create final result file
				direct_connection_final = os.path.join(self.base.out, f"{self.base.name}.direct.final.fa")

				with open(direct_connection_final, "w") as f:
					f.write(f">{self.base.name}_direct_overlap\n{final_seq}\n")

				# Record to log
				logfile.write(f"Found direct overlap between left and right sequences!\n")
				logfile.write(f"Overlap region: left sequence {s1}-{e1}, right sequence {s2}-{e2}\n")
				logfile.write(f"Overlap length: {length}bp, identity: {identity}%\n")
				logfile.write(f"GAP Length: {gap_length}\n")  # Negative value indicates overlap
				logfile.write(f"Direct connection sequence saved to: {direct_connection_final}\n")
				logfile.write("GAP can be closed!\n")  # Important marker for extract_filled_gap_info detection
				logfile.write("Direct overlap\n")  # Mark this as direct overlap for subsequent identification
				logfile.write("But will enter extension phase first, if extension phase fails or result gap length is negative, then use direct connection result\n")
				logfile.write("============================\n\n")

				# Set direct connection flag and save result, but don't set endSignal (allow entering extension phase)
				self.has_direct_connection = True
				self.direct_connection_result = {
					"final_seq_file": direct_connection_final,
					"gap_length": gap_length,
					"overlap_length": length,
					"identity": identity
				}
				
				return True
		
			# If found qualifying tandem repeat alignments
			if alignments:
				# Sort by alignment length and identity, take best result
				best_alignment = sorted(alignments, key=lambda x: (x[4], x[5]), reverse=True)[0]
				s1, e1, s2, e2, length, identity, parts = best_alignment

				# Calculate gap length (negative value indicates overlap)
				gap_length = -(length)  # Overlap length as negative gap length

				# Extract more complete repeat unit
				# We need to find complete repeat unit, not just matching region
				# Repeat unit length should be distance between consecutive matching regions

				# Look for second match if it exists
				second_alignments = []
				with open(mummer_output, "r") as f:
					for line in f:
						if line.startswith("#"):
							continue
						parts = line.strip().split("\t")
						if len(parts) < 10:
							continue
						
						try:
							ss1, se1, ss2, se2 = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
							s_length = int(parts[4])
							s_identity = float(parts[6])
							
							# Find matches that are different from the best match but have similar patterns
							if s_length >= 20 and s_identity >= 80:
								# Ensure this is not the same match
								if not (ss1 == s1 and se1 == e1 and ss2 == s2 and se2 == e2):
									second_alignments.append((ss1, se1, ss2, se2, s_length, s_identity))
						except ValueError:
							continue

				# Find possible repeat unit length from second match
				repeat_unit_length = None

				if second_alignments:
					# Sort by match length and similarity
					second_best = sorted(second_alignments, key=lambda x: (x[4], x[5]), reverse=True)[0]
					ss1, se1, ss2, se2, s_length, s_identity = second_best

					# Calculate distance between two matches
					if s1 < ss1:  # First match is on the left
						distance = ss1 - s1
					else:  # Second match is on the left
						distance = s1 - ss1

					if distance > 0:
						repeat_unit_length = distance
						logfile.write(f"Estimated repeat unit length based on distance between two matches: {repeat_unit_length}bp\n")

				# If no second match found, try using match length as repeat unit length
				if repeat_unit_length is None:
					repeat_unit_length = length
					logfile.write(f"No second match found, using match length as repeat unit length: {repeat_unit_length}bp\n")
				
				# Create final sequence file
				direct_connection_final = os.path.join(self.base.out, f"{self.base.name}.direct.final.fa")

				# Build final sequence based on repeat unit length
				# We will keep all left-end sequences, then add a complete repeat unit at the end
				repeat_unit = None

				# Try to extract repeat unit
				if repeat_unit_length > 0:
					if s1 + repeat_unit_length <= len(left_seq):
						repeat_unit = left_seq[s1:s1+repeat_unit_length]
					elif s2 + repeat_unit_length <= len(right_seq):
						repeat_unit = right_seq[s2:s2+repeat_unit_length]

				# If unable to extract valid repeat unit, use matching region as repeat unit
				if repeat_unit is None or len(repeat_unit) < 20:  # Prevent repeat unit from being too short
					repeat_unit = left_seq[s1:e1]
					repeat_unit_length = len(repeat_unit)
					logfile.write(f"Using matching region as repeat unit, length: {repeat_unit_length}bp\n")

				# Build final sequence
				final_seq = left_seq + repeat_unit + right_seq[e2:]
				
				with open(direct_connection_final, "w") as f:
					f.write(f">{self.base.name}_tandem_repeat\n{final_seq}\n")
				
				# Record to log
				logfile.write(f"Found tandem repeats between left and right sequences!\n")
				logfile.write(f"Matching region: left sequence {s1}-{e1}, right sequence {s2}-{e2}\n")
				logfile.write(f"Match length: {length}bp, similarity: {identity}%\n")
				logfile.write(f"Repeat unit length: {repeat_unit_length}bp\n")
				logfile.write(f"GAP Length: {gap_length}\n")  # Negative value indicates overlap
				logfile.write(f"Direct connection sequence saved to: {direct_connection_final}\n")
				logfile.write("GAP can be closed!\n")  # Important marker for extract_filled_gap_info detection
				logfile.write("Tandem repeat\n")  # Mark this as tandem repeat
				logfile.write("But will enter extension phase first, if extension phase fails or result gap length is negative, then use direct connection result\n")
				logfile.write("============================\n\n")
				
				# Set direct connection flag and save result, but don't set endSignal (allow entering extension phase)
				self.has_direct_connection = True
				self.direct_connection_result = {
					"final_seq_file": direct_connection_final,
					"gap_length": gap_length,
					"overlap_length": length,
					"identity": identity,
					"is_tandem_repeat": True,
					"repeat_unit_length": repeat_unit_length
				}
				
				return True
		
		# If no valid direct connection found
		logfile.write("No valid direct connection found\n")
		logfile.write("============================\n\n")
		return False

