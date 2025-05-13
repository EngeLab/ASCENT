import pysam
import sys
import argparse
import os

def calculate_metrics(bam_input):
    """
    Calculate fragment length and overlap for each read pair.
    Works with both filename or stdin ("-")
    """
    bam = pysam.AlignmentFile(bam_input, "rb")
    
    current_reads = []
    current_name = None
    
    try:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
                
            if current_name != read.query_name:
                if len(current_reads) == 2:
                    r1, r2 = sorted(current_reads, key=lambda x: x.reference_start if x.reference_start is not None else float('inf'))
                    
                    # Skip if either read is unmapped or missing positions
                    if (r1.reference_name == r2.reference_name and 
                        not r1.is_unmapped and not r2.is_unmapped and
                        r1.reference_start is not None and r1.reference_end is not None and
                        r2.reference_start is not None and r2.reference_end is not None):
                        
                        frag_length = r2.reference_end - r1.reference_start
                        overlap = max(0, min(r1.reference_end, r2.reference_end) - 
                                     max(r1.reference_start, r2.reference_start))
                        
                        print(f"{frag_length}\t{overlap}")
                
                current_reads = [read]
                current_name = read.query_name
            else:
                current_reads.append(read)
        
        # Process last pair
        if len(current_reads) == 2:
            r1, r2 = sorted(current_reads, key=lambda x: x.reference_start if x.reference_start is not None else float('inf'))
            if (r1.reference_name == r2.reference_name and 
                not r1.is_unmapped and not r2.is_unmapped and
                r1.reference_start is not None and r1.reference_end is not None and
                r2.reference_start is not None and r2.reference_end is not None):
                
                frag_length = r2.reference_end - r1.reference_start
                overlap = max(0, min(r1.reference_end, r2.reference_end) - 
                             max(r1.reference_start, r2.reference_start))
                print(f"{frag_length}\t{overlap}")
                
    except BrokenPipeError:
        sys.stderr.close()
        sys.stdout = open(os.devnull, 'w')
        
    finally:
        bam.close()

def main():
    parser = argparse.ArgumentParser(description='Calculate fragment length and overlap')
    parser.add_argument('bam', nargs='?', default='-', 
                      help='Input BAM file (name-sorted) or - for stdin')
    
    args = parser.parse_args()
    calculate_metrics(args.bam)

if __name__ == "__main__":
    main()
