import re
import os
import sys
import math
import Bio
from Bio import SeqIO
import shutil
import GapFiller
from GapFiller import GapFiller
import selectRawReads
from selectRawReads import selectRawReads
import subprocess
import glob
import FindExtensionReads
from FindExtensionReads import FindExtensionReads
import FindExtensionContigs
from FindExtensionContigs import FindExtensionContigs
import GapFillerClass
from GapFillerClass import GapFillerClass, OutputSequence, InputSequence, mummer

class TelFiller:
    """
    TelFiller class - A tool for extending sequences toward a set of potential terminators,
    specifically designed for filling telomeric regions.
    """
    
    def __init__(self, parameter, kparameters):
        """
        Initializes the TelFiller object.
        
        Args:
            parameter (list): A list of operational parameters.
            kparameters (list): A list of k-mer related parameters.
        """
        # Parse and set parameters
        self._parse_parameters(parameter, kparameters)

        # Create a 'base' object to mimic the structure GapFillerClass expects
        self._create_base_context()
        
        # Create work directories
        if not os.path.exists(self.workDir):
            os.makedirs(self.workDir)
        if not os.path.exists(self.base.outfile):
            os.makedirs(self.base.outfile)

        # Print initialization summary
        self._print_summary()

        # Main workflow
        self.check_direct_connection()
        
        if not self.has_direct_connection:
            self.run_telfiller()
        else:
            print("TelFiller process finished: Direct connection found.")

    def _parse_parameters(self, parameter, kparameters):
        """Helper to parse input parameter lists."""
        self.mode = parameter[0]
        self.remove = parameter[1]
        self.thread = parameter[2]
        self.reads = parameter[3]
        self.workDir = parameter[4]
        self.seqleft = parameter[5]
        self.seqright = parameter[6]
        self.flag = parameter[7]
        self.edge = parameter[8]
        self.filterDepth = parameter[9]
        self.MaximumExtensionLength = parameter[10]
        self.readsdict = parameter[11] if len(parameter) > 11 else None
        self.maxReadsLen = parameter[12] if len(parameter) > 12 else None
        self.seedLen = parameter[13] if len(parameter) > 13 else None
        
        self.kmer_size, self.kmer_num, self.kmer_length, self.j = kparameters
        self.resume_round = None # Placeholder for future resume functionality

    def _create_base_context(self):
        """Creates a context object to maintain compatibility with GapFiller components."""
        class Base:
            pass
        self.base = Base()
        self.base.mode, self.base.remove, self.base.thread, self.base.reads = self.mode, self.remove, self.thread, self.reads
        self.base.out = self.workDir
        self.base.seqleft, self.base.seqright, self.base.flag, self.base.edge = self.seqleft, self.seqright, self.flag, self.edge
        self.base.tag = self.flag
        self.base.filterDepth, self.base.MaximunExtensionLength = self.filterDepth, self.MaximumExtensionLength
        self.base.readsDict, self.base.maxReadsLen, self.base.seedLen = self.readsdict, self.maxReadsLen, self.seedLen
        self.base.kmer_size, self.base.kmer_num, self.base.kmer_length, self.base.j = self.kmer_size, self.kmer_num, self.kmer_length, self.j
        self.base.name = os.path.basename(self.workDir)
        self.base.outfile = os.path.join(self.workDir, "process")

    def _print_summary(self):
        """Prints a summary of the initialized parameters."""
        print("--- TelFiller Initialized ---")
        print(f"  Work Directory: {self.workDir}")
        print(f"  Reads File: {self.reads}")
        print(f"  Start Sequence: {self.seqleft if self.flag == 'left' else self.seqright}")
        print(f"  End Sequence(s): {self.seqright if self.flag == 'left' else self.seqleft}")
        print(f"  Extension Direction: {self.flag}")
        print(f"  K-mer Params: size={self.kmer_size}, num={self.kmer_num}, length={self.kmer_length}")
        print("-----------------------------")

    def check_direct_connection(self):
        """
        Checks if the start sequence can be directly connected to any of the end sequences.
        This handles a 1-vs-N comparison.
        """
        print("Checking for direct connections...")
        self.has_direct_connection = False
        direct_conn_dir = os.path.join(self.workDir, "direct_connection")
        os.makedirs(direct_conn_dir, exist_ok=True)
        
        log_file = os.path.join(self.workDir, "process.log")
        with open(log_file, 'a') as logfile:
            logfile.write("\n========== Direct Connection Check ==========\n")
            
            source_file = self.base.seqleft if self.flag == 'left' else self.base.seqright
            target_file = self.base.seqright if self.flag == 'left' else self.base.seqleft

            try:
                source_record = next(SeqIO.parse(source_file, "fasta"))
                target_records = list(SeqIO.parse(target_file, "fasta"))
            except (StopIteration, FileNotFoundError) as e:
                logfile.write(f"Error reading sequence files: {e}\n")
                return

            logfile.write(f"Source: {source_record.id} ({len(source_record.seq)} bp)\n")
            logfile.write(f"Targets: {len(target_records)} sequences\n")

            connected_pairs = []
            for target_record in target_records:
                temp_source_file = os.path.join(direct_conn_dir, "temp_source.fa")
                temp_target_file = os.path.join(direct_conn_dir, "temp_target.fa")
                SeqIO.write(source_record, temp_source_file, "fasta")
                SeqIO.write(target_record, temp_target_file, "fasta")

                mummer_prefix = os.path.join(direct_conn_dir, f"{source_record.id}_vs_{target_record.id}")
                mummer_output = mummer(temp_source_file, temp_target_file, mummer_prefix)
                
                # Analyze mummer results for a valid link at the sequence ends
                if os.path.exists(mummer_output) and os.path.getsize(mummer_output) > 0:
                    with open(mummer_output, 'r') as f:
                        for row in f:
                            parts = row.strip().split('\t')
                            if len(parts) < 8 or not parts[0].strip().isdigit(): continue
                            
                            e1, s2, length, identity = int(parts[1]), int(parts[2]), int(parts[4]), float(parts[6])
                            len1, len2 = int(parts[7]), int(parts[8])
                            
                            # Check for overlap at the correct ends with sufficient quality
                            if e1 >= len1 - self.base.edge and s2 <= self.base.edge and length >= 100 and identity >= 95:
                                connected_pairs.append({'row': row, 'source': source_record, 'target': target_record})
                                logfile.write(f"  Found potential link: {source_record.id} -> {target_record.id}, len={length}, id={identity}%\n")

            if connected_pairs:
                # Select the best connection based on length and identity
                best_pair = sorted(connected_pairs, key=lambda p: (int(p['row'].split('\t')[4]), float(p['row'].split('\t')[6])), reverse=True)[0]
                self.has_direct_connection = True
                self._finalize(status="direct_connection", link_info=best_pair)
                logfile.write(f"Success! Best direct connection found: {best_pair['source'].id} -> {best_pair['target'].id}\n")
            else:
                logfile.write("No direct connections found.\n")

    def run_telfiller(self):
        """Executes the main iterative extension workflow."""
        print("No direct connection found. Starting iterative extension...")
        
        # Initialize loop state
        self.roundNum = 1
        self.usedReads, self.lastRoundUsedReads = [], []
        self.extensionLen = 0
        self.endSignal = False
        
        self.roundInputSeq = self.base.seqleft if self.flag == 'left' else self.base.seqright
        self.terminalSeqFile = self.base.seqright if self.flag == 'left' else self.base.seqleft

        while not self.endSignal and self.roundNum <= 100:
            self._setup_round()
            
            # Create a mock Elongation object for this round
            class MockElongation: pass
            elongation_context = MockElongation()
            elongation_context.base = self.base
            elongation_context.roundDir = self.roundDir
            elongation_context.roundInputSeq = self.roundInputSeq
            elongation_context.usedReads = self.usedReads
            elongation_context.lastRoundUsedReads = self.lastRoundUsedReads
            elongation_context.extensionLen = self.extensionLen
            elongation_context.endSignal = False

            # --- Core Extension Steps ---
            input_seq_obj = InputSequence(elongation_context)
            ext_reads = FindExtensionReads(input_seq_obj, self.lastRoundUsedReads, self.usedReads, self.base.kmer_size, self.base.kmer_num, self.base.kmer_length, self.base.j, self.base.out)

            if ext_reads.note or ext_reads.extensionReadsNum == 0:
                self._handle_stop_condition("stalled_no_reads", self.roundInputSeq)
                break

            ext_contigs = FindExtensionContigs(ext_reads)
            if 'No extension contigs or reads found' in ext_contigs.selectContigNote:
                self._handle_stop_condition("stalled_no_contigs", self.roundInputSeq)
                break

            output_seq_obj = OutputSequence(ext_contigs, elongation_context)
            
            # Check if the newly extended sequence can link to a terminal sequence
            is_linked, link_info = self._check_for_link(output_seq_obj.outputSequence)
            if is_linked:
                self._finalize(status="success", link_info=link_info)
                self.endSignal = True
            else:
                # Prepare for the next round
                self.roundInputSeq = output_seq_obj.outputSequence
                self.lastRoundUsedReads = output_seq_obj.ExtensionUsedReads
                self.usedReads.extend([r for r in self.lastRoundUsedReads if r not in self.usedReads])
                self.extensionLen += ext_contigs.extensionLength
                self.roundNum += 1
                
                # Check for other stop conditions
                if self.base.MaximunExtensionLength and self.extensionLen >= self.base.MaximunExtensionLength:
                    self._handle_stop_condition("max_length", self.roundInputSeq)
    
        if self.roundNum > 100:
            self._handle_stop_condition("max_rounds", self.roundInputSeq)

    def _setup_round(self):
        """Initializes directories and logs for the current extension round."""
        print(f"\n--- Starting Round {self.roundNum} (Total Extension: {self.extensionLen} bp) ---")
        self.roundDir = os.path.join(self.base.outfile, f"round{self.roundNum}")
        os.makedirs(self.roundDir, exist_ok=True)

    def _check_for_link(self, extended_seq_file):
        """Checks if the extended sequence can link to any of the terminal sequences."""
        print(f"Checking for link between {os.path.basename(extended_seq_file)} and terminal sequences...")
        best_link = None
        
        for term_record in SeqIO.parse(self.terminalSeqFile, "fasta"):
            temp_term_file = os.path.join(self.roundDir, "temp_terminal.fa")
            SeqIO.write(term_record, temp_term_file, "fasta")
            
            mummern = os.path.join(self.roundDir, f"link_check.mummer.{term_record.id}")
            mummerout = mummer(extended_seq_file, temp_term_file, mummern)
            
            if os.path.exists(mummerout) and os.path.getsize(mummerout) > 0:
                with open(mummerout, 'r') as f:
                    for row in f:
                        parts = row.strip().split('\t')
                        if len(parts) < 8 or not parts[0].strip().isdigit(): continue
                        e_ref, s_query, length, identity = int(parts[1]), int(parts[2]), int(parts[4]), float(parts[6])
                        len_ref, len_query = int(parts[7]), int(parts[8])
                        
                        if e_ref >= len_ref - self.base.edge and s_query <= self.base.edge and length >= 100 and identity >= 95:
                            current_link = {'row': row, 'extended_seq_file': extended_seq_file, 'target': term_record}
                            if best_link is None or float(parts[6]) > float(best_link['row'].split('\t')[6]):
                                best_link = current_link
        
        if best_link:
            return True, best_link
        return False, None

    def _handle_stop_condition(self, reason, final_file):
        """Centralized handler for loop termination."""
        stop_messages = {
            "stalled_no_reads": "No new extension reads found.",
            "stalled_no_contigs": "Could not assemble valid extension contigs.",
            "max_length": "Reached maximum extension length.",
            "max_rounds": "Reached maximum number of rounds (100)."
        }
        print(f"Stopping extension: {stop_messages.get(reason, 'Unknown reason.')}")
        self.endSignal = True
        self._finalize(status=reason, final_sequence_file=final_file)

    def _finalize(self, status, link_info=None, final_sequence_file=None):
        """Creates the final output sequence file."""
        print(f"Finalizing process with status: {status.upper()}")
        final_fa = os.path.join(self.workDir, "final.fa")
        
        final_seq = None
        header = f">telomere_filled_{status}"

        if status in ["success", "direct_connection"] and link_info:
            parts = link_info['row'].strip().split('\t')
            # In direct_connection, source is extended, target is terminal.
            # In success, extended_seq_file is extended, target is terminal.
            extended_seq = link_info['source'].seq if status == "direct_connection" else next(SeqIO.parse(link_info['extended_seq_file'], 'fasta')).seq
            terminal_seq = link_info['target'].seq
            
            # Alignment: extended_seq (ref) vs terminal_seq (query)
            ref_end, query_start, query_end = int(parts[1]), int(parts[2]), int(parts[3])

            if query_start < query_end: # Forward match
                final_seq = extended_seq[:ref_end] + terminal_seq[query_end:]
            else: # Reverse-complement match
                terminal_revcomp = terminal_seq.reverse_complement()
                final_seq = extended_seq[:ref_end] + terminal_revcomp[query_start:]
            
            header += f"_linked_to_{link_info['target'].id}"

        elif final_sequence_file and os.path.exists(final_sequence_file):
            final_seq = next(SeqIO.parse(final_sequence_file, "fasta")).seq
            print(f"Process stopped. The longest achieved sequence is being saved.")

        if final_seq:
            with open(final_fa, 'w') as f:
                f.write(f"{header}\n{str(final_seq)}\n")
            print(f"Final sequence saved to: {final_fa}")
        else:
            print("Error: Could not generate a final sequence file.")