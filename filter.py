import pysam

from FamilyRecord import FamilyRecord


def filterBAM(inFile: pysam.AlignmentFile, outFile: pysam.AlignmentFile, args):
    pass

def getArgs(parser):

    # SNV arguments
    snvArgs = parser.add_argument_group("SNV Argumentss")

    snvArgs.add_argument(
        "-tb", "--target_bed",
        required=False,
        help="A tab-delinated file listing regions on which variant calling will be restricted to"
    )
    snvArgs.add_argument(
        "-vaft", "--variant_allele_fraction_threshold",
        default=0.01,
        type=float,
        help="Minimum variant frequency threshold for each strand [Default: %(default)s]"
    )
    snvArgs.add_argument(
        "-mo", "--min_molecules",
        default=40,
        type=int,
        help="Number of total molecules (supporting or otherwise) required to call a variant at this position. Reduce this if you are running only on positions you expect to be mutated [Default: %(default)s]"
    )
    snvArgs.add_argument(
        "-mum", "--mutant_molecules",
        default=3,
        required=False,
        type=int,
        help="Number of TOTAL molecules (i.e. on the forward and reverse strand) required to call a variant as real (set to 0 if you are forcing variant calling at known sites) [Default: %(default)s]"
    )
    snvArgs.add_argument(
        "-mrpu", "--min_reads_per_uid",
        default=2,
        type=int,
        help="Bases with support between MRPU and SSBT will be classified as a weak supported base [Default: %(default)s]"
    )
    snvArgs.add_argument(
        "-ssbt", "--strong_supported_base_threshold",
        default=3,
        type=int,
        help="Bases with support equal to or greater then SSBT, will be classified as a strong supported base [Default: %(default)s]"
    )

    snvArgs.add_argument(
        "-eds", "--enforce_dual_strand",
        action='store_true',
        help="require at least one molecule to be read in each of the two directions"
    )
    snvArgs.add_argument(
        "-mq", "--min_qual",
        default=30,
        type=int,
        help="Minimum base quality threshold, below which a WEAK base will be ignored'")

    # Filter Args
    filterArgs = parser.add_argument_group("Filter Arguments")
    filterArgs.add_argument("-ss", "--allow_single_stranded", action="store_true", default=False,
                            help="Allow variants with only single stranded support [Default: %(default)s]")
    filterArgs.add_argument("-sb", "--strand_bias_threshold", default=0.05, type=float,
                            help="Strand bias p-value threshold, below which vairants will be discarded [Default: %(default)s]")
    filterArgs.add_argument("-st", "--strong_base_threshold", default=1, type=int,
                            help="Strong supported base count threshold [Default: %(default)s]")
    filterArgs.add_argument("-wt", "--weak_base_threshold", default=2, type=int,
                            help="Weak supported base count theshold [Default: %(default)s]")
    filterArgs.add_argument("-md", "--min_depth", type=int, default=2,
                            help="Minimum depth threshold [Default: %(default)s]")