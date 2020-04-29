This workflow was used to extract spike gene variants observed in GISAID SARS-CoV-2 sequences relative to the [reference sequence](https://www.ncbi.nlm.nih.gov/protein/1796318598) in [Pinto et al, __Structural and functional analysis of a potent sarbecovirus neutralizing antibody__, bioRxiv (2020)](https://www.biorxiv.org/content/10.1101/2020.04.07.023903v3). Tested on macOS 10.15.4.

## Environment setup

Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

Then, create a Conda environment:

```bash
conda env create -f environment.yml
```

Finally, activate the newly created environment:

```bash
conda activate gisaid-epicov
```

## Obtain SARS-CoV-2 isolate sequences from GISAID

Download CoV-2 genome sequences from [GISAID](https://www.gisaid.org/) after enabling the `complete (>29,000bp)` and `low coverage excl` filters. In the manuscript, we used a GISAID data snapshot obtained on April 27th 2020 at 10.00am Pacific time (n=11,839). Due to GISAID license terms, this step cannot be programmatically automated and users will need to obtain a GISAID data dump individually. We gratefully acknowledge the authors, originating and submitting laboratories of the sequences from GISAID’s EpiFlu Database on which this research is based.

Name the GISAID sequence dump `gisaid_cov2020_sequences.fa`.

## Cleanup GISAID sequences

```bash
cat gisaid_cov2020_sequences.fa | mac2unix | seqkit seq -u | seqkit replace -p " " -r "_" | seqkit replace -s -p "_" -r "-" | seqkit replace -s -p "   " -r "---" | seqkit replace -s -p " " -r "" | seqkit replace -s -p "U" -r "T" | seqkit replace -s -p "\?" -r "N" | seqkit grep -v -r -p "\/pangolin\/|\/bat\/" >gisaid_clean.fa
```

## Generate S protein alignments

Standard workflow:

```bash
cat gisaid_clean.fa | seqkit subseq -r 18000:28000 | seqkit replace -s -p "-" -r "N" >exonerate_input.fa

exonerate --cores 4 -m protein2dna --refine full --minintron 999999 --percent 30 --showalignment false --showvulgar false --ryo ">%ti\n%tcs" -q YP_009724390.1.fa -t exonerate_input.fa | grep -E -v "exonerate|Command line|Hostname" | seqkit translate -f 1 | grep -v "^1 1" | seqkit replace -p "_frame=1 " -r "" | seqkit sort -l -r | seqkit rmdup -n >s_aa.fa

mafft --thread 4 --amino --bl 80 --nomemsave --reorder --add s_aa.fa --keeplength YP_009724390.1.fa | seqkit grep -v -p "ref" >s_msa_default.fa
```

An alternative workflow based on genomic alignments:

```bash
mafft --thread 4 --nuc --anysymbol --nomemsave --adjustdirection --reorder --add gisaid_clean.fa --keeplength NC_045512.2.fa | seqkit grep -v -p "ref" >genome_msa.fa

cat genome_msa.fa | seqkit subseq -r 21563:25381 | seqkit translate -f 1 | grep -v "^1 1" | seqkit replace -p "_frame\=1 " -r "" >s_msa_alternative.fa
```

## Extract variants

Run the following R code:

```r
library(Biostrings)

ref <- readAAStringSet("YP_009724390.1.fa")
ref <- unname(unlist(strsplit(as.character(ref), "")))

msa <- readAAMultipleAlignment("s_msa_default.fa", format = "fasta")

if(ncol(msa) != length(ref)) stop("Mismatch between reference sequence and sequence alignment")

cm <- consensusMatrix(msa)
cm <- cm[intersect(rownames(cm), AA_STANDARD), ]
cm <- t(t(cm) / colSums(cm))

variants <- data.frame(
    coordinate = integer(),
    ref = character(),
    alt = character(),
    frequency = numeric()
)

for(i in 1:length(ref)) {
    for(j in rownames(cm)[cm[, i] > 0]) {
        if(j != ref[i]) {
            variants <- rbind(variants, data.frame(coordinate = i, ref = ref[i], alt = j, frequency = signif(cm[j, i], 3)))
        }
    }
}

write.csv(variants, file = "variants.csv", quote = F, row.names = F)
```

Variants can also be extracted from the alternative genomic alignment workflow by replacing `s_msa_default.fa` with `s_msa_alternative.fa` in the R code. Similarly, this workflow can be used with sequences arranged in a  multi-entry FASTA files other than GISAID's.
