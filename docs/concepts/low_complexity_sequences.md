# Low-complexity sequences

**What it is:** Stretches of DNA made of very few distinct patterns, e.g. `AAAAAAAA` or `ATATATAT`.

**Biology / intuition:** Real genomes contain repeats (microsatellites, poly-A tails), but in sequencing data low-complexity regions cause trouble: reads align to many places equally well, assemblers get confused, and BLAST gives misleading hits. They can also be artefacts (adapter dimers, homopolymer errors). Tools like DUST and RepeatMasker flag them.

**How we detect it:** Slide a window of k = 4 along the sequence, collect every 4-mer, and compute (number of unique 4-mers) ÷ (total 4-mers). If that ratio is under 0.30, flag the sequence.

**When it's wrong to use:**
- **The ratio depends on length, and this is a real bug.** There are only 4⁴ = 256 possible 4-mers, so the unique count can never exceed 256. Once a sequence is longer than about 860 bp (256 ÷ 0.30), the ratio is *mathematically forced* below 0.30, so **every sequence over roughly 860 bp is flagged low-complexity**, however diverse. It only looks fine on our sample data because those sequences are short (median 140 bp). Verify this yourself in W2.
- Ratio of unique to total k-mers is a heuristic, not a standard method like DUST (which scores triplet frequency in windows).
- It flags a whole sequence, not the low-complexity *region* inside it.

**Where in our code:** `_check_low_complexity` in [fasta_qc.py](../../src/bioagent/pipelines/fasta_qc.py). `LOW_COMPLEXITY_THRESHOLD = 0.30`.

**Likely interview questions**
1. *Why do we care about low-complexity regions?* — They make alignment and assembly ambiguous and can inflate hits, so QC flags them and tools often mask them.
2. *What's a weakness of your detector?* — The unique/total 4-mer ratio is capped at 256 ÷ length, so any sequence over about 860 bp is flagged whatever its content. A windowed method like DUST is the right approach. *(Earn in W2.)*
