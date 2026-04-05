#!/usr/bin/env python3
"""
Generate a synthetic antigen CSV for testing.
Creates sequences with a His6-ABP N-terminal tag + random antigen region,
with some intentional similarities planted for validation.
"""

import random
import csv

random.seed(42)

AA = "ACDEFGHIKLMNPQRSTVWY"

TAG = "MHHHHHHGSSG"  # His6-ABP placeholder — replace with your real tag

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
        full_seq = TAG + variant
        records.append({
            "antigen_id": f"FAM{i:02d}_VAR{j:02d}",
            "sequence":   full_seq,
        })

# Fill remaining with random sequences
for i in range(len(records), 500):
    length = random.randint(40, 120)
    full_seq = TAG + random_seq(length)
    records.append({
        "antigen_id": f"RND_{i:04d}",
        "sequence":   full_seq,
    })

random.shuffle(records)

output = "data/test_antigens.csv"
with open(output, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["antigen_id", "sequence"])
    writer.writeheader()
    writer.writerows(records)

print(f"Generated {len(records)} test antigens → {output}")
print(f"Tag used: {TAG}")
print(f"Planted {len(families)} similarity families (10 variants each)")
print(f"\nUpdate DEFAULT_TAG in tag_stripper.py to match: {TAG}")
