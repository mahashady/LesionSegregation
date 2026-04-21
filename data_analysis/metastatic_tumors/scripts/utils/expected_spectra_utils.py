from collections import defaultdict
import itertools
import json
import os
import csv

def load_json(filepath: str):
    with open(filepath, "r") as f:
        return json.load(f)

def ensure_parent_dir(filepath: str) -> None:
    parent = os.path.dirname(filepath)
    if parent:
        os.makedirs(parent, exist_ok=True)

def load_sample_mutation_rates(filepath: str) -> dict[str, float]:
    mutation_rates = {}
    with open(filepath, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            context = row["mutation"].strip()
            rate = float(row["Mutrate"])
            mutation_rates[context] = rate
    return mutation_rates

def group_by_trinuc(mutation_rates):
    """
    Group mutation-rate entries by trinucleotide context.
    """    
    grouped = defaultdict(list)
    for label, rate in mutation_rates.items():
        trinuc = label[0:3]
        grouped[trinuc].append((label, rate))
    return grouped

def compute_expected_triallelic(mutation_rates, genome_context):
    """
    Compute expected triallelic counts from a dictionary of per-context mutation rates.

    For each trinucleotide context, all unordered pairs of distinct mutation
    channels are considered. The expected number of sites with both events is:

        genome_count(trinuc) * rate1 * rate2

    and this expected count is added to both component labels.
    """    
    grouped = group_by_trinuc(mutation_rates)
    expected = defaultdict(float)

    for trinuc, mutations in grouped.items():
        context_count = genome_context.get(trinuc)
        if context_count is None:
            continue

        for (label1, rate1), (label2, rate2) in itertools.combinations(mutations, 2):
            joint_prob = rate1 * rate2
            expected_count = context_count * joint_prob
            expected[label1] += expected_count
            expected[label2] += expected_count

    return expected

def parse_signature_context(raw: str) -> str:
    """
    Convert signature notation like A[C>A]A into ACA>A.
    Adjust this parser if your input format differs.
    """
    raw = str(raw)
    if len(raw) < 7:
        raise ValueError(f"Unexpected mutation context format: {raw}")
    return raw[0] + raw[2] + raw[6] + ">" + raw[4]

