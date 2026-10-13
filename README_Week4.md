# PRA2003 – Week 4: Bacteria per event in 5 million events

Counting bacterial strains per event in a sample of 5M events, with
statistical uncertainties from the sub-sampling method.

## Contents

1. [Overview](#1-overview)
2. [Files](#2-files)
3. [How to run](#3-how-to-run)
4. [The data](#4-the-data)
5. [Method](#5-method)
6. [Results](#6-results)
7. [Notes and limitations](#7-notes-and-limitations)

---

## 1. Overview

The goal this week is to measure the **average number of each bacterial
strain per event**, together with its **statistical uncertainty**, using
the entire sample of 5M events (10 files of 500K events each). The
uncertainties are calculated with the **sub-sampling method** from the
lecture.

Each strain appears as a normal (wild-type) form with a positive ID and a
mutant form with the same ID and a minus sign (for example E. coli WT = 211
and E. coli mutant = -211).

## 2. Files

| File | Purpose |
|---|---|
| `README.md` | Current file |
| `Week4Sub_Sampling.R` | Counts every ID in each of the 10 data files, creating all of the sub-samples, which are outputted in a csv file as a list |
| `Week4Analysis.R` | Reads the list and computes the final averages and uncertainties |
| `sub_sample_results.csv` | Output of "Week4Sub_Sampling.R": count and average per ID for each file |
| `final_results.csv` | Output of "Week4Analysis.R": final average ± uncertainty per ID |

The raw data files (`output-Set1.txt` … `output-Set10.txt`, about 800 MB
each) are **not** included in the repository due to their size. They
can be downloaded here --> https://surfdrive.surf.nl/files/index.php/s/7udCnWTk4yMUASD.

## 3. How to run

1. Put the 10 raw data files in the same folder as the scripts.
2. Open `Week4Sub_Sampling.R` in VS Code and press Run Code (Ctrl + Alt + N).
   This reads the raw data and writes `sub_sample_results.csv`.
3. Open `Week4Analysis.R` in VS Code and press Run Code (Ctrl + Alt + N).
   This reads `sub_sample_results.csv` and writes `final_results.csv`.

Part 2 only needs the CSV from part 1, so the final results can be
reproduced without the raw data files.

## 4. The data

- 10 files, each containing about 500K events (5M events in total).
- Each event starts with a header line (event number and number of
  bacteria), followed by one line per bacterium: three momentum components
  and the bacterial ID.
- The IDs of the 12 strains of interest are:
- Empty events are excluded from the average

| ID | Strain | ID | Strain |
|---:|---|---:|---|
| 211 | E. coli WT | -211 | E. coli mutant |
| 321 | Bacillus subtilis WT | -321 | Bacillus subtilis mutant |
| 2212 | Pseudomonas aeruginosa WT | -2212 | Pseudomonas aeruginosa antibiotic-resistant |
| 3122 | Streptococcus pneumoniae | -3122 | Capsule-deficient S. pneumoniae |
| 3312 | Mycobacterium tuberculosis | -3312 | Drug-resistant M. tuberculosis |
| 3334 | Salmonella enterica | -3334 | Salmonella mutant |

The data contain further IDs that are not one of these 12 strains. They
appear as "Unknown" in `final_results.csv` and are not part of the result
below.

## 5. Method

The sub-sampling technique:

1. **Split** the sample into 10 sub-samples of similar size. Each data file
   (about 500K events) is one sub-sample.
2. **Analyse each sub-sample separately:** for every ID, compute the
   average number per event in that file.
3. **Central value:** the inverse-variance weighted average of the 10
   sub-sample results (weight = 1 / uncertainty² per sub-sample).
4. **Statistical uncertainty:** the standard deviation of the 10
   sub-sample results.

The standard deviation describes how much the result varies from one
sub-sample of about 500K events to another. Only the statistical
uncertainty is quoted. No systematic uncertainty was evaluated, and since
the result is a number per event, it has no unit.

**Asymmetry test (Section 6.2/6.3):** for each pair X vs -X, the
difference is computed **within each of the 10 sub-samples first**
(differenceᵢ = averageᵢ(X) − averageᵢ(−X)), and only then combined: the
reported difference is the mean of the 10 per-file differences, and its
uncertainty is the standard deviation of those 10 differences. This is
different from combining the two strains' already-finished uncertainties
by quadrature (√(σ_X² + σ_-X²)), which assumes X and -X vary
independently from file to file. They do not — both are counted in the
same events, so a file with more bacteria overall tends to have more of
both forms together. Taking the difference within each file first cancels
that shared fluctuation before it is combined, which is why this method is
used instead of quadrature.

## 6. Results

### 6.1 Average count per event (5M events, 10 sub-samples)

| ID | Strain | Total count | Average / event | Uncertainty |
|---:|---|---:|---:|---:|
| 211 | E. coli WT | 92,126,688 | 19.9495 | ± 0.0327 |
| -211 | E. coli mutant | 91,977,542 | 19.9172 | ± 0.0319 |
| 321 | Bacillus subtilis WT | 11,587,227 | 2.50915 | ± 0.00477 |
| -321 | Bacillus subtilis mutant | 11,560,946 | 2.50346 | ± 0.00550 |
| 2212 | Pseudomonas aeruginosa WT | 5,578,693 | 1.20803 | ± 0.00190 |
| -2212 | P. aeruginosa antibiotic-resistant | 5,468,447 | 1.18416 | ± 0.00241 |
| 3122 | Streptococcus pneumoniae | 1,277,330 | 0.276599 | ± 0.00107 |
| -3122 | Capsule-deficient S. pneumoniae | 1,254,690 | 0.271696 | ± 0.000985 |
| 3312 | Mycobacterium tuberculosis | 182,139 | 0.0394412 | ± 0.000284 |
| -3312 | Drug-resistant M. tuberculosis | 180,104 | 0.0390005 | ± 0.000402 |
| 3334 | Salmonella enterica | 5,482 | 0.00118710 | ± 0.0000417 |
| -3334 | Salmonella mutant | 5,318 | 0.00115158 | ± 0.0000508 |

The data also contain 26 other IDs that are not one of the 12 strains. They
are listed as "Unknown" in `final_results.csv`.

### 6.2 Asymmetry: wild type vs mutant

| Pair | X - (-X) | n σ | Asymmetry A | Asymmetric? |
|---|---:|---:|---:|:---:|
| 211 vs -211 (E. coli) | 0.0323 ± 0.0045 | **7.15** | **(0.08 ± 0.01)%** | **yes** |
| 321 vs -321 (B. subtilis) | 0.00569 ± 0.00329 | 1.73 | (0.11 ± 0.07)% | no |
| 2212 vs -2212 (P. aeruginosa) | 0.02387 ± 0.00237 | **10.06** | **(1.00 ± 0.10)%** | **yes** |
| 3122 vs -3122 (S. pneumoniae) | 0.00490 ± 0.00058 | **8.41** | **(0.89 ± 0.11)%** | **yes** |
| 3312 vs -3312 (M. tuberculosis) | 0.00044 ± 0.00049 | 0.90 | (0.56 ± 0.62)% | no |
| 3334 vs -3334 (Salmonella) | 0.000036 ± 0.000069 | 0.52 | (1.5 ± 2.9)% | no |

The averages of X and -X are in table 6.1.

#### The answer for each pair, and why

- **211 vs -211, E. coli: yes, a significant asymmetry.**
  - The WT is 0.0323 per event more common, and once the difference is
    taken per file (accounting for the strong file-to-file correlation
    between the two counts), the uncertainty shrinks to 0.0045 — a
    **7.15σ** effect.
  - A = **(0.08 ± 0.01)%**. Small in absolute terms, but measured
    precisely because E. coli is so abundant.
- **321 vs -321, Bacillus subtilis: no asymmetry.**
  - The difference is **1.73σ**, and A = (0.11 ± 0.07)% is consistent with
    zero.
  - Any real asymmetry is smaller than about 0.3%.
- **2212 vs -2212, Pseudomonas aeruginosa: yes, a clear asymmetry.**
  - The wild type is 0.0239 per event more common than the
    antibiotic-resistant strain.
  - That is **10.06σ**, far above the 3σ threshold.
  - A = **(1.00 ± 0.10)%**.
- **3122 vs -3122, Streptococcus pneumoniae: yes, an asymmetry.**
  - The wild type is more common than the capsule-deficient strain by
    **8.41σ**.
  - A = **(0.89 ± 0.11)%**.
- **3312 vs -3312, Mycobacterium tuberculosis: no asymmetry.**
  - The difference is **0.90σ**, and A = (0.56 ± 0.62)% is consistent with
    zero.
  - With about 180,000 counts per strain, the data can only rule out
    asymmetries larger than about 2.4%.
- **3334 vs -3334, Salmonella: no asymmetry.**
  - The difference is **0.52σ**.
  - With only about 5,400 counts per strain, A = (1.5 ± 2.9)% is very
    imprecise, so an asymmetry smaller than about 10% cannot be ruled out.

### 6.3 Asymmetry: the other 12 particle/antiparticle pairs

The same test applied to the 12 other pairs of IDs in the data:

| Pair | X - (-X) | n σ | Asymmetry A | Asymmetric? |
|---|---:|---:|---:|:---:|
| 3212 vs -3212 | 0.00173 ± 0.00076 | 2.26 | (0.57 ± 0.25)% | no |
| 3222 vs -3222 | 0.00173 ± 0.00059 | 2.93 | (0.58 ± 0.20)% | no |
| 3112 vs -3112 | 0.00202 ± 0.00046 | **4.43** | **(0.68 ± 0.15)%** | **yes** |
| 3322 vs -3322 | 0.00026 ± 0.00030 | 0.88 | (0.34 ± 0.39)% | no |
| 431 vs -431 | -0.00005 ± 0.00017 | 0.32 | (-0.23 ± 0.74)% | no |
| 531 vs -531 | 0.00002 ± 0.00008 | 0.32 | (1.3 ± 3.9)% | no |
| 4232 vs -4232 | 0.000003 ± 0.000033 | 0.09 | (0.2 ± 2.7)% | no |
| 4132 vs -4132 | 0.000008 ± 0.000041 | 0.19 | (0.7 ± 3.5)% | no |
| 5132 vs -5132 | 0.000014 ± 0.000020 | 0.70 | (14 ± 19)% | no |
| 5232 vs -5232 | -0.000002 ± 0.000016 | 0.14 | (-2 ± 15)% | no |
| 4332 vs -4332 | 0.000004 ± 0.000012 | 0.37 | (13 ± 36)% | no |
| 5332 vs -5332 | 0.0000006 ± 0.0000035 | 0.18 | (18 ± 86)% | no |

- **1 of the 12 pairs is above 3σ:** 3112 vs -3112, at 4.43σ, with
  A = (0.68 ± 0.15)%. 3222 (2.93σ) is close but does not cross the
  threshold.
- **The rarest pairs** (5132, 5232, 4332, 5332, with 7–260 counts each)
  have such large uncertainties that they cannot show anything either way.

---

## 7. Notes and limitations

- **Meaning of the uncertainty.** The quoted uncertainty is the spread of
  one 500K-event sub-sample. The uncertainty of the average over all 10
  sub-samples is smaller by a factor of about √10 (the standard error of
  the mean).
- **Wild type vs mutant.** For each strain the wild type and the mutant
  have similar averages (for example E. coli: 19.949 vs 19.917). The counts
  of a wild type and its mutant come from the same events and can be
  correlated, so the uncertainty on their difference is not obtained by
  simply adding the two uncertainties in quadrature. A difference is
  computed within each sub-sample and the spread taken over those values
  (see Section 5).
- **Empty events.** Events with 0 bacteria were excluded from the event
  count used in the average (though they were still scanned to check
  their header). Across the 10 sub-samples, 4,617,993 of the nominal
  5,000,000 events contained at least one bacterium and were used in the
  averages; the remaining 382,007 events (about 7.6%) were empty and
  excluded.
- **Statistical only.** No systematic uncertainty was estimated.
