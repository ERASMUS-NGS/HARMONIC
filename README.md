# HARMONIC: Hierarchy-based Active-tree Redundancy-Minimized ONtology Integrative Comparator
HARMONIC is an R package for Gene Ontology (GO) Biological Process enrichment analysis that reduces redundancy in top-ranked results by restructuring enriched terms into ontology-guided term trees based on the GO DAG. It defines and filters active trees using coverage/consistency criteria to suppress structurally driven false positives arising from hierarchical dependencies. HARMONIC supports single- and multiple- case study designs by integrating enrichment outputs across conditions, providing scoring options and visualization utilities to summarize common and condition-specific biological processes. Optional features include ORA-style soft filtering and fold-change-aware weighting for tree and term prioritization.

![Alt text](./FIG/HARMONIC_F1.jpg "HARMONIC")
>- **Problem:** GO Biological Process enrichment results are often dominated by biologically similar top terms, and GO’s hierarchical DAG + gene sharing can create **structural false positives**, reducing interpretability and confidence.
>- **Solution (HARMONIC):** An analytical framework that reorganizes GO outputs into **ontology-informed term trees** and applies **structural filtering** based on **within-tree consistency**.
>- **Key idea 1 — Tree organizing:** HARMONIC groups similar terms into trees to reduce redundancy and compress results into interpretable units.
>- **Key idea 2 — Active-tree filtering:** HARMONIC suppresses isolated, single-term–driven signals using an **active-tree criterion**.
>- **Multiple-case support:** HARMONIC aligns results from multiple comparisons into a **shared tree structure** and classifies **common vs condition-specific** trees using **HWES** and relative signal strength across conditions.
>- **Parameter recommendations:** A parameter search provides **input-size–specific recommended settings**.
>- **Benchmarking:** Evaluated across **31 single-case analyses** against **ORA** and **GSEA** using **adjusted p-values**, **fraction of terms assigned**, and **rich factors** to assess interpretability.
>- **Validation:** Demonstrated reproduction of reported signals across public datasets spanning **knockout**, **cancer**, **differentiation**, **dose–response**, and **cohort** designs.
>- **Usability:** Provides **R-based visualization functions** for both single- and multiple-case settings to aid interpretation.
>- **Take-home:** HARMONIC reduces redundancy and structural false positives driven by the GO DAG, enabling **reproducible term-pattern summarization** for multi-comparison studies.

## Installation

### Linux / macOS (Using conda OR micromamba at terminal)

```bash
micromamba create -n HARMONIC
micromamba activate HARMONIC

micromamba install -c conda-forge -c r -c bioconda -y \
jupyter_core jupyter_client jupyterlab_pygments jupyter_server r-irkernel jupyterlab r=4.3.1

micromamba install -c conda-forge -c r -c bioconda -y \
r-png r-data.table r-systemfonts r-gdtools r-ggforce r-ggiraph bioconductor-xvector bioconductor-sparsearray \
bioconductor-biostrings bioconductor-delayedarray bioconductor-summarizedexperiment bioconductor-annotationdbi \
bioconductor-go.db bioconductor-keggrest bioconductor-fgsea bioconductor-deseq2 bioconductor-gosemsim pigz \
bioconductor-dose bioconductor-enrichplot bioconductor-clusterprofiler bioconductor-ggtree \
bioconductor-org.mm.eg.db bioconductor-org.hs.eg.db

mkdir path_to_HARMONICanalysis
wget -P path_to_HARMONICanalysis https://github.com/user-attachments/files/26292175/HARMONIC_FORMAT.tar.gz
### ex) wget -P /home/RNA/gitHARMONIC/ https://github.com/user-attachments/files/26292175/HARMONIC_FORMAT.tar.gz

cd path_to_HARMONICanalysis
pigz -dc -p 4 HARMONIC_FORMAT.tar.gz | tar -xf -

mv HARMONIC_FORMAT HARMONIC_project_name
### ex) mv HARMONIC_FORMAT HARMONIC_DIFF

cd HARMONIC_project_name/HARMONIC
pwd
### result of pwd is filepath we use in R

cd count
cp path_to_input_countMATRIX RAWcount.txt
### ex) cp RAWcount_DIFF.txt RAWcount.txt
```

### Input count matrix (`RAWcount.txt`)

- Tab-separated text file with a header row.
- First column: gene identifiers. Remaining columns: raw read counts for each sample.
- Sample columns must be **grouped by condition in the same order as `full_condition`**, with the number of columns per condition matching `number_of_rep`.

## Usage

### 1. Run HARMONIC

```r
# install.packages("remotes")
remotes::install_github("ERASMUSlab/HARMONIC")
library(haRmonic)

HARMONIC_ANALYSIS(
  filepath       = "/path/to/HARMONIC_project_name/HARMONIC",
  type           = "standard",
  full_condition = c("Control", "TreatmentA", "TreatmentB"),
  number_of_rep  = c(3, 3, 3),
  DEG_list_name  = "input_DEG_list.txt",
  mango_design   = c("TreatmentA_Control_UP", "TreatmentB_Control_UP"),
  core           = 4,
  ref_genome     = "mm",
  PASSED_NUM     = 3,
  PASSED_RATIO   = 25,
  similarity     = 60,
  FC             = 2,
  preprocessing  = "T"
)
```

#### Parameters

- `filepath` — absolute path to the `HARMONIC` folder of your project (the output of `pwd` in the setup step).
- `type` — DEG calling mode: `"standard"` (adjusted P < 0.05, |fold change| ≥ 2) or `"broad"` (adjusted P < 0.05, |fold change| ≥ 1.5).
- `full_condition` — names of all experimental conditions, in the same order as the sample columns of `RAWcount.txt`. Condition names must not contain underscores (`_`).
- `number_of_rep` — number of replicates for each condition, in the same order and length as `full_condition`.
- `DEG_list_name` — file name of the DEG list inside `filepath` (`input_DEG_list.txt` in the provided template).
- `mango_design` — comparisons to analyze, written as `<A>_<B>_UP` (genes up-regulated in A relative to B) or `<A>_<B>_DOWN` (genes down-regulated in A relative to B). `<A>` and `<B>` must be names listed in `full_condition`; `<B>` is used as the reference condition.
- `core` — number of CPU cores to use.
- `ref_genome` — reference organism: `"mm"` (mouse) or `"hs"` (human).
- `PASSED_NUM` — minimum number of enriched GO terms required for a tree to be retained as an active tree.
- `PASSED_RATIO` — minimum percentage of enriched GO terms within a tree required for it to be retained as an active tree.
- `similarity` — overlap threshold (%) above which active trees are merged into a single structure.
- `FC` — HWES ratio threshold between comparisons; trees whose HWES differs by at least this factor are classified as condition-specific, otherwise as common.
- `preprocessing` — `"T"` to run differential expression and GO enrichment before tree construction, or `"F"` to reuse existing preprocessing outputs.

**Recommended starting values:** `PASSED_NUM = 3`, `PASSED_RATIO = 25` and `similarity = 60`, as used in the HARMONIC manuscript.

### 2. Adjust parameters

Parameters can be tuned to the size of the enrichment input and the desired level of compression:

- **Increase `PASSED_NUM` or `PASSED_RATIO`** — stricter active-tree filtering; fewer trees are retained, each supported by more enriched terms. Decrease them for more permissive filtering.
- **Decrease `similarity`** — trees are merged at a lower overlap, producing fewer and broader structures. Increase it to merge less and keep more, narrower structures.
- **Decrease `FC`** — smaller HWES differences are sufficient to call a tree condition-specific, so more trees are classified as condition-specific and fewer as common. Increase it for the opposite effect.
- **Use `type = "broad"`** — if a comparison yields fewer than two active trees with `"standard"`, the more permissive DEG threshold can recover additional enrichment input.
- **Set `preprocessing = "F"`** — skip differential expression and GO enrichment when their outputs from a previous run are already present in `filepath`, for example when only the tree parameters are changed.

### 3. Dynamic analysis

Dynamic analysis identifies active trees that are selectively active in a chosen set of comparisons, for example the late time points of a differentiation time course, compared with the remaining comparisons. A tree is reported when it is active in every selected comparison and its HWES in the selected comparisons exceeds that in the remaining comparisons by at least `FC` (both minimum versus maximum and mean versus mean).

```r
HARMONIC_ANALYSIS(
  filepath          = "/path/to/HARMONIC_project_name/HARMONIC",
  type              = "standard",
  full_condition    = c("Day0", "Day4", "Day7", "Day14"),
  number_of_rep     = c(3, 3, 3, 3),
  DEG_list_name     = "input_DEG_list.txt",
  mango_design      = c("Day4_Day0_UP", "Day7_Day0_UP", "Day14_Day0_UP"),
  core              = 4,
  ref_genome        = "mm",
  PASSED_NUM        = 3,
  PASSED_RATIO      = 25,
  similarity        = 60,
  FC                = 2,
  condition         = 2:3,
  dynamic_analyisis = "T",
  preprocessing     = "F"
)
```

- `condition` — 1-based indices of the comparisons in `mango_design` to treat as the selected set. Use `2:3` for a consecutive range, or `c(1, 3)` for specific, non-consecutive comparisons.
- `dynamic_analyisis` — `"T"` to run dynamic analysis, `"F"` (default) to skip it. Note that the argument name is spelled `dynamic_analyisis` in the current release.

Results are written to `HARMONIC_SEPERATE_forMULTI_range.txt` in `filepath`.

## Documentation

Full documentation and tutorials:
https://erasmus-ngs.github.io/HARMONIC


## Citation

If you use HARMONIC in your research, please cite:
https://erasmuslab.github.io/HARMONIC/authors.html#citation
