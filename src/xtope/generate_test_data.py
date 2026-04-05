#!/usr/bin/env python3
"""
Generate a synthetic antigen CSV for testing.
Creates random antigen sequences with some intentional similarities planted
for validation. Sequences are generated without any N-terminal tag — N-tag
stripping is opt-in and not assumed for test data.
"""

import random
import csv

random.seed(42)

AA = "ACDEFGHIKLMNPQRSTVWY"


def random_seq(length):
    return "".join(random.choices(AA, k=length))


def mutate(seq, rate=0.1):
    """Introduce random point mutations."""
    return "".join(
        c if random.random() > rate else random.choice(AA)
        for c in seq
    )


records = []

# Plant 5 "families" of similar antigens (should cluster together)
families = []
for i in range(5):
    base = random_seq(random.randint(40, 80))
    families.append(base)
    for j in range(10):
        variant = mutate(base, rate=0.05 + random.random() * 0.1)
        records.append({
            "antigen_id": f"FAM{i:02d}_VAR{j:02d}",
            "sequence":   variant,
        })

# Fill remaining with random sequences
for i in range(len(records), 500):
    length = random.randint(40, 120)
    records.append({
        "antigen_id": f"RND_{i:04d}",
        "sequence":   random_seq(length),
    })

random.shuffle(records)

output = "data/test_antigens.csv"
with open(output, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["antigen_id", "sequence"])
    writer.writeheader()
    writer.writerows(records)

print(f"Generated {len(records)} test antigens → {output}")
print(f"Planted {len(families)} similarity families (10 variants each)")
print(f"Note: sequences contain no N-terminal tag (N-tag stripping is opt-in via --tag)")
