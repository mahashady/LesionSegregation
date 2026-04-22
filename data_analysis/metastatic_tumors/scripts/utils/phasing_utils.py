import os

def build_bam_path(sample: str, subset: str) -> str:
    return os.path.join(
        BAM_BASE_DIRS[subset],
        f"{sample}T_dedup.realigned.minisorted.bam",
    )


def build_intersect1_path(sample: str, chrom: str, pos_multi: int, pos_germline: int, subset: str) -> str:
    return os.path.join(
        BAM_BASE_DIRS[subset],
        f"{sample}T_{chrom}:{pos_multi}-{pos_germline}_dedup.realigned.mini.intersect1.bam",
    )


def build_intersect2_path(sample: str, chrom: str, pos_multi: int, pos_germline: int, subset: str) -> str:
    return os.path.join(
        BAM_BASE_DIRS[subset],
        f"{sample}T_{chrom}:{pos_multi}-{pos_germline}_dedup.realigned.mini.intersect2.bam",
    )


def classify_phasing(germline_counts: list[int]) -> str:
    """
    germline_counts order:
      [n_ref_ref, n_alt_ref, n_ref_alt1, n_alt_alt1, n_ref_alt2, n_alt_alt2]
    """
    if (
        germline_counts[2] != 0
        and germline_counts[3] == 0
        and germline_counts[4] != 0
        and germline_counts[5] == 0
    ):
        return "Phased"
    elif (
        germline_counts[2] == 0
        and germline_counts[3] != 0
        and germline_counts[4] == 0
        and germline_counts[5] != 0
    ):
        return "Phased"
    elif (
        germline_counts[2] == 0
        and germline_counts[3] != 0
        and germline_counts[4] != 0
        and germline_counts[5] == 0
    ):
        return "Anti"
    elif (
        germline_counts[2] != 0
        and germline_counts[3] == 0
        and germline_counts[4] == 0
        and germline_counts[5] != 0
    ):
        return "Anti"
    else:
        return "No_phasing"


def phase(
    sample: str,
    chrom: str,
    pos_multi: int,
    ref_multi: str,
    alt1_multi: str,
    alt2_multi: str,
    pos_germline: int,
    ref_germline: str,
    alt_germline: str,
    subset: str,
) -> list[int]:
    bam_path = build_bam_path(sample, subset)
    intersect1_path = build_intersect1_path(sample, chrom, pos_multi, pos_germline, subset)
    intersect2_path = build_intersect2_path(sample, chrom, pos_multi, pos_germline, subset)

    samfile = pysam.AlignmentFile(bam_path, "rb")
    output1 = pysam.AlignmentFile(intersect1_path, "wb", template=samfile)
    output2 = pysam.AlignmentFile(intersect2_path, "wb", template=samfile)

    for read in samfile.fetch(chrom, pos_multi - 1, pos_multi):
        output1.write(read)
    output1.close()

    pysam.index(intersect1_path)

    output1 = pysam.AlignmentFile(intersect1_path, "rb")
    for read in output1.fetch(chrom, pos_germline - 1, pos_germline):
        output2.write(read)
    output1.close()

    samfile.close()
    output2.close()

    pysam.index(intersect2_path)

    samfile = pysam.AlignmentFile(intersect2_path, "rb")

    read_bases = {}
    if pos_multi < pos_germline:
        start = pos_multi - 1
        end = pos_germline - 1
    else:
        start = pos_germline - 1
        end = pos_multi - 1

    for pileupcolumn in samfile.pileup(chrom, start, end):
        if pileupcolumn.pos == (pos_multi - 1) or pileupcolumn.pos == (pos_germline - 1):
            print(pileupcolumn.n)
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del or pileupread.is_refskip or pileupread.query_position is None:
                    continue

                read_name = pileupread.alignment.query_name
                if read_name not in read_bases:
                    read_bases[read_name] = ["N", "N"]

                base = pileupread.alignment.query_sequence[pileupread.query_position]
                if pileupcolumn.pos == (pos_multi - 1):
                    read_bases[read_name][0] = base
                if pileupcolumn.pos == (pos_germline - 1):
                    read_bases[read_name][1] = base

    samfile.close()

    read_bases_filtered = {}
    for read_name, bases in read_bases.items():
        if bases[0] != "N" and bases[1] != "N":
            read_bases_filtered[read_name] = bases

    print("Len dict = " + str(len(read_bases)))
    print("Len dict = " + str(len(read_bases_filtered)))

    ref_list = []
    alt1_list = []
    alt2_list = []

    for read_name in read_bases_filtered:
        if read_bases_filtered[read_name][0] == alt1_multi:
            alt1_list.append(read_bases_filtered[read_name][1])
        elif read_bases_filtered[read_name][0] == alt2_multi:
            alt2_list.append(read_bases_filtered[read_name][1])
        elif read_bases_filtered[read_name][0] == ref_multi:
            ref_list.append(read_bases_filtered[read_name][1])

    n_ref_ref_list = ref_list.count(ref_germline)
    n_alt_ref_list = ref_list.count(alt_germline)
    n_ref_alt1_list = alt1_list.count(ref_germline)
    n_alt_alt1_list = alt1_list.count(alt_germline)
    n_ref_alt2_list = alt2_list.count(ref_germline)
    n_alt_alt2_list = alt2_list.count(alt_germline)

    return [
        n_ref_ref_list,
        n_alt_ref_list,
        n_ref_alt1_list,
        n_alt_alt1_list,
        n_ref_alt2_list,
        n_alt_alt2_list,
    ]
