import argparse
from bam_tool_rep.bam_tools import load_read_bam_multi
import pickle


def main():
    parser = argparse.ArgumentParser(
        description="Filter and load reads from a BAM file with customizable thresholds."
    )
    parser.add_argument(
        "bam",
        type=str,
        help="Path to the BAM file."
    )
    parser.add_argument(
        "--remove-shorter-than",
        type=int,
        default=5000,
        help="Remove reads shorter than this length (in bp). Default: 5000."
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        default=None,
        help="Maximum number of reads to process. Default: None (process all)."
    )
    parser.add_argument(
        "--brdu-threshold",
        type=float,
        default=0.05,
        help="Lower mean BrdU threshold. Remove reads with lower BrdU content. Default: 0.05."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=3,
        help="Number of threads to use. Default: 3."
    )
    parser.add_argument(
        "--res",
        type=int,
        default=100,
        help="Resolution parameter for load_read_bam_multi. Default: 100."
    )
    parser.add_argument(
        "--check-size",
        action="store_true",
        help="Print the number of reads returned."
    )

    args = parser.parse_args()

    r = load_read_bam_multi(
        args.bam,
        res=args.res,
        remove_less_than={"b": args.brdu_threshold},
        remove_shorter_than=args.remove_shorter_than,
        maxi=args.max_reads,
        threads=args.threads,
    )

    pickle_file = args.bam + ".pickle"
    with open(pickle_file, "wb") as f:
        pickle.dump(r, f)


    if args.check_size:
        print(f"Number of reads returned: {len(r)}")
        print(f"Saved results to {pickle_file}")
    else:
        print(f"Reads loaded successfully and saved to {pickle_file}.")


if __name__ == "__main__":
    main()
