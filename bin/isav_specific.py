import os
import sys
import re
import time
from pathlib import Path
from collections import defaultdict

# Add package directory to Python path
# home = str(Path.home())
# sys.path.append(f'{home}/git/gitlab/vsnp3/bin')
from download_fasta_by_acc import Downloader


class ISAV_Specific():

    def order_segments_by_coverage(self, virus_dict):
        """
        Orders virus segments based on segment/VP numbers and selects highest coverage when duplicates exist.
        
        Args:
            virus_dict (dict): Dictionary containing virus sequences with headers and coverage info
            
        Returns:
            dict: Ordered dictionary with segments arranged by genome position

        Template"
            When using this file as a template update the 
                segment_order dictionary
                if 1 <= potential_segment <= 10:
                for segment_num in range(1, 11):
                
            - Keeping the quotes around the segment names helps with specificity
            - Note that the FASTA headers are being converted to upper case so keep segment names in upper case
            to match the FASTA headers
        """
        # Segment to genome position mapping
        segment_order = {
            '(PB2)': 1, '(PB1)': 2, '(NP)': 3, '(PA)': 4, '(F)': 5,
            '(HE)': 6, '(NS1)': 7, '(NEP)': 7, '(M1)': 8, 
            'SEGMENT 1': 1, 'SEGMENT 2': 2, 'SEGMENT 3': 3, 'SEGMENT 4': 4,
            'SEGMENT 5': 5, 'SEGMENT 6': 6, 'SEGMENT 7': 7, 'SEGMENT 8': 8,
            'SEG 1': 1, 'SEG 2': 2, 'SEG 3': 3, 'SEG 4': 4, 'SEG 5': 5,
            'SEG 6': 6, 'SEG 7': 7, 'SEG 8': 8
        }
        
        # Initialize dictionary to store sequences by segment number
        segments = defaultdict(list)
        
        # First pass: Classify sequences by segment number
        for seq_id, seq_info in virus_dict.items():
            header = seq_info['header'].upper()
            segment_num = None
            
            # Check for VP/NS notation
            for seg_name, seg_num in segment_order.items():
                if seg_name in header:
                    segment_num = seg_num
                    break
            
            # Check for direct segment number notation if no VP/NS found
            if segment_num is None:
                # Look for patterns like "segment 1" or just "1" with word boundaries
                segment_matches = re.findall(r'\b(?:SEGMENT\s*)?(\d{1,2})\b', header)
                if segment_matches:
                    try:
                        potential_segment = int(segment_matches[0])
                        if 1 <= potential_segment <= 8:  # Valid segment range
                            segment_num = potential_segment
                    except ValueError:
                        continue
            
            if segment_num:
                segments[segment_num].append((seq_id, seq_info))
        
        # Create ordered dictionary selecting highest coverage for each segment
        ordered_dict = {}
        
        # Process segments in order (1-8)
        for segment_num in range(1, 9):
            if segment_num in segments:
                # Sort by percent coverage and take the highest
                best_sequence = max(segments[segment_num], 
                                key=lambda x: x[1]['percent_covered'])
                
                # Add segment information to the sequence info
                seq_info = best_sequence[1].copy()
                seq_info['segment_number'] = segment_num
                
                # Try to determine segment name from header
                header = best_sequence[1]['header'].upper()
                segment_name = f"Segment {segment_num}"
                for seg_name, seg_num in segment_order.items():
                    if seg_num == segment_num and seg_name in header:
                        segment_name = seg_name
                        break
                        
                seq_info['segment_name'] = segment_name
                ordered_dict[best_sequence[0]] = seq_info
        
        return ordered_dict
    
    def concatenate_fasta_files(self, ordered_dict, output_file="concatenated_specific.fasta"):
        """
        Concatenates downloaded FASTA files based on the ordered dictionary.
        Preserves segment order and adds segment information to headers.
        
        Args:
            ordered_dict (dict): Ordered dictionary containing sequence information
            output_file (str): Name of the output file (default: concatenated_sequences.fasta)
        """
        with open(output_file, 'w') as outfile:
            for seq_id, seq_info in ordered_dict.items():
                # Construct input filename - assuming it matches the accession ID
                input_file = f"{seq_id}.fasta"
                
                try:
                    with open(input_file, 'r') as infile:
                        # Read the content of the input file
                        content = infile.read().strip()
                        
                        # Split into header and sequence
                        parts = content.split('\n', 1)
                        if len(parts) != 2:
                            print(f"Warning: Unexpected format in {input_file}")
                            continue
                            
                        header, sequence = parts
                        
                        # Add segment information to header if not present
                        segment_info = f"segment_{seq_info['segment_number']}"
                        if 'segment_name' in seq_info:
                            segment_info = f"{segment_info}_{seq_info['segment_name']}"
                        
                        # Write modified header and sequence
                        outfile.write(f"{header} {segment_info}\n")
                        outfile.write(f"{sequence}\n")
                        
                except FileNotFoundError:
                    print(f"Warning: Could not find file for {seq_id}")
                    continue
                    
        print(f"Concatenation complete. Output written to {output_file}")

    def run(self, alignment_stats=None):
        dict_ordered = self.order_segments_by_coverage(alignment_stats)
        
        print(f"Downloading {len(dict_ordered)} ISAV segments...")
        downloaded_count = 0
        failed_downloads = []
        
        for i, (seq_id, seq_inf) in enumerate(dict_ordered.items(), 1):
            print(f"Downloading segment {i}/{len(dict_ordered)}: {seq_id}")
            
            # Add rate limiting - wait 2 seconds between downloads (except for the first)
            if i > 1:
                print("  Waiting 2 seconds (NCBI rate limiting)...")
                time.sleep(2)
            
            downloader = Downloader(seq_id)
            max_retries = 3
            success = False
            
            # Retry logic with exponential backoff
            for attempt in range(max_retries):
                try:
                    downloader.fasta()
                    success = True
                    downloaded_count += 1
                    print(f"  ✓ Successfully downloaded {seq_id}")
                    break
                except Exception as e:
                    print(f"  ⚠ Attempt {attempt + 1} failed for {seq_id}: {e}")
                    if attempt < max_retries - 1:
                        wait_time = (attempt + 1) * 5  # 5, 10, 15 seconds
                        print(f"    Waiting {wait_time} seconds before retry...")
                        time.sleep(wait_time)
            
            if not success:
                failed_downloads.append(seq_id)
                print(f"  ✗ Failed to download {seq_id} after {max_retries} attempts")
        
        print(f"\nDownload summary:")
        print(f"  Successfully downloaded: {downloaded_count}/{len(dict_ordered)}")
        if failed_downloads:
            print(f"  Failed downloads: {len(failed_downloads)}")
            print(f"  Failed accessions: {', '.join(failed_downloads)}")
        
        # Only proceed with concatenation if we have some successful downloads
        if downloaded_count > 0:
            print("Proceeding with concatenation of successfully downloaded segments...")
            # Filter out failed downloads from the ordered dict
            successful_dict = {k: v for k, v in dict_ordered.items() if k not in failed_downloads}
            self.concatenate_fasta_files(successful_dict)
            
            # Clean up individual FASTA files
            for seq_id, seq_inf in successful_dict.items():
                try:
                    os.remove(f'{seq_id}.fasta')
                except FileNotFoundError:
                    pass  # File already removed or never existed
        else:
            print("Warning: No segments were successfully downloaded")
            # Create an empty output file to prevent downstream errors
            with open("concatenated_specific.fasta", 'w') as f:
                f.write(">No_Segments_Available\nNNNNNNNNNN\n")