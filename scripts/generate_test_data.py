#!/usr/bin/env python3
"""
Generate synthetic FASTQ data for benchmarking rustmasker vs bbmask.sh

Creates reads with varying complexity:
- Low complexity (homopolymers, dinucleotide repeats)
- Medium complexity (simple repeats)
- High complexity (random sequences)

Usage:
    ./generate_test_data.py --num-reads 10000 --read-length 1000 -o test.fastq
    ./generate_test_data.py --num-reads 1000 --read-length 100000 -o long_reads.fastq  # ONT-like
"""

import argparse
import random
import sys
from typing import TextIO


def generate_low_complexity(length: int) -> str:
    """Generate low-complexity sequences (homopolymers and dinucleotide repeats)"""
    pattern_type = random.choice(['homopolymer', 'dinucleotide', 'trinucleotide'])

    if pattern_type == 'homopolymer':
        base = random.choice('ACGT')
        return base * length
    elif pattern_type == 'dinucleotide':
        bases = random.choice(['AT', 'CG', 'AC', 'GT', 'AG', 'CT'])
        repeat = (bases * (length // 2 + 1))[:length]
        return repeat
    else:  # trinucleotide
        bases = ''.join(random.choices('ACGT', k=3))
        repeat = (bases * (length // 3 + 1))[:length]
        return repeat


def generate_medium_complexity(length: int) -> str:
    """Generate medium-complexity sequences (short repeats)"""
    repeat_unit = ''.join(random.choices('ACGT', k=random.randint(4, 8)))
    return (repeat_unit * (length // len(repeat_unit) + 1))[:length]


def generate_high_complexity(length: int) -> str:
    """Generate high-complexity random sequences"""
    return ''.join(random.choices('ACGT', k=length))


def generate_quality_string(length: int, quality: int = 30) -> str:
    """Generate quality scores (Phred+33)"""
    qual_char = chr(33 + quality)  # Default Q30
    return qual_char * length


def write_fastq_record(fh: TextIO, read_id: int, sequence: str, quality: str):
    """Write a single FASTQ record"""
    fh.write(f"@read_{read_id}\n")
    fh.write(f"{sequence}\n")
    fh.write("+\n")
    fh.write(f"{quality}\n")


def generate_fastq(num_reads: int, read_length: int, output_file: str,
                   low_pct: float = 0.3, medium_pct: float = 0.3):
    """
    Generate FASTQ file with mixed complexity reads

    Args:
        num_reads: Total number of reads to generate
        read_length: Length of each read (can vary for ONT-like data)
        output_file: Output FASTQ filename
        low_pct: Percentage of low-complexity reads
        medium_pct: Percentage of medium-complexity reads
        (remainder will be high complexity)
    """

    # Calculate how many reads of each type
    num_low = int(num_reads * low_pct)
    num_medium = int(num_reads * medium_pct)
    num_high = num_reads - num_low - num_medium

    print(f"Generating {num_reads} reads of length {read_length}:")
    print(f"  - Low complexity: {num_low} ({low_pct*100:.1f}%)")
    print(f"  - Medium complexity: {num_medium} ({medium_pct*100:.1f}%)")
    print(f"  - High complexity: {num_high} ({(1-low_pct-medium_pct)*100:.1f}%)")

    with open(output_file, 'w') as fh:
        read_id = 0

        # Generate low complexity reads
        for _ in range(num_low):
            # Vary length slightly for realism (±10%)
            length = int(read_length * random.uniform(0.9, 1.1))
            sequence = generate_low_complexity(length)
            quality = generate_quality_string(length)
            write_fastq_record(fh, read_id, sequence, quality)
            read_id += 1

        # Generate medium complexity reads
        for _ in range(num_medium):
            length = int(read_length * random.uniform(0.9, 1.1))
            sequence = generate_medium_complexity(length)
            quality = generate_quality_string(length)
            write_fastq_record(fh, read_id, sequence, quality)
            read_id += 1

        # Generate high complexity reads
        for _ in range(num_high):
            length = int(read_length * random.uniform(0.9, 1.1))
            sequence = generate_high_complexity(length)
            quality = generate_quality_string(length)
            write_fastq_record(fh, read_id, sequence, quality)
            read_id += 1

    # Calculate file size
    import os
    size_mb = os.path.getsize(output_file) / (1024 * 1024)
    print(f"\nGenerated {output_file} ({size_mb:.2f} MB)")
    print(f"Total bases: {num_reads * read_length:,}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate synthetic FASTQ data for benchmarking',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate 10K short reads (Illumina-like)
  %(prog)s --num-reads 10000 --read-length 150 -o illumina_test.fastq

  # Generate 1K long reads (ONT-like)
  %(prog)s --num-reads 1000 --read-length 50000 -o ont_test.fastq

  # Generate 100 very long reads (stress test)
  %(prog)s --num-reads 100 --read-length 200000 -o ont_long_test.fastq
        """
    )

    parser.add_argument('-n', '--num-reads', type=int, default=10000,
                        help='Number of reads to generate (default: 10000)')
    parser.add_argument('-l', '--read-length', type=int, default=1000,
                        help='Average read length in bases (default: 1000)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output FASTQ file')
    parser.add_argument('--low-complexity', type=float, default=0.3,
                        help='Fraction of low-complexity reads (default: 0.3)')
    parser.add_argument('--medium-complexity', type=float, default=0.3,
                        help='Fraction of medium-complexity reads (default: 0.3)')

    args = parser.parse_args()

    # Validate percentages
    if args.low_complexity + args.medium_complexity > 1.0:
        print("Error: low_complexity + medium_complexity must be <= 1.0", file=sys.stderr)
        sys.exit(1)

    generate_fastq(
        num_reads=args.num_reads,
        read_length=args.read_length,
        output_file=args.output,
        low_pct=args.low_complexity,
        medium_pct=args.medium_complexity
    )


if __name__ == '__main__':
    main()
