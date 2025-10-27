Germline Variant Calling with GATK HaplotypeCaller
===================================================

# Overview

Germline variant calling identifies inherited variants present in all cells of an organism. Unlike somatic calling where we compare tumor to normal tissue, germline calling focuses on finding variants that differ from the reference genome but are consistent across all cells in our sample.

The GATK (Genome Analysis ToolKit) Best Practices workflow is the gold standard for germline variant calling. In this tutorial, we'll walk through the complete pipeline, from raw sequencing reads to filtered, annotated variants ready for downstream analysis.

**What you'll learn:**

1. How to align reads and prepare BAM files for variant calling
2. Why base quality score recalibration (BQSR) matters and how to apply it
3. The difference between single-sample and joint calling (and why joint calling is better)
4. How to filter variants using VQSR or hard filters
5. Quality control metrics to validate your callset
6. Common pitfalls and how to avoid them

**Expected variant counts for a human genome:**
- Whole genome sequencing (WGS): ~3-4 million SNVs, ~500k indels (generally, a roughly 10:1 ratio of SNVs to indels is expected)
- Whole exome sequencing (WES): ~25-35k SNVs in coding regions (again, roughly 10:1 SNVs to indels)

These numbers assume high-quality data from a sample of European ancestry. Samples from other populations may show different numbers due to reference bias (the reference genome is based primarily on European samples).

# Prerequisites and setup

## Reference genome and resources

You'll need the following files, all from the **same genome build** (mixing builds causes serious problems):

```bash
# Set up your reference paths
REF=~/references/Homo_sapiens_assembly38.fasta
DBSNP=~/references/dbsnp_138.hg38.vcf.gz
MILLS=~/references/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
KNOWN_INDELS=~/references/Homo_sapiens_assembly38.known_indels.vcf.gz
```

**What are these files?**
- `Homo_sapiens_assembly38.fasta`: The reference genome (GRCh38/hg38)
- `dbsnp_138.hg38.vcf.gz`: Database of known SNPs, used for annotation and BQSR
- `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`: High-confidence indels for BQSR
- `Homo_sapiens_assembly38.known_indels.vcf.gz`: Additional known indels for BQSR

All of these are available from the [Broad Resource Bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0).

## Reference indices

Before we can use the reference, we need to create several index files:

```bash
# BWA index (only needed once per reference)
bwa index ${REF}
# Runtime: ~1.5 hours on a modern machine

# Samtools faidx index
samtools faidx ${REF}
# Runtime: ~1 minute

# GATK sequence dictionary
gatk CreateSequenceDictionary -R ${REF}
# Runtime: ~2 minutes
```

**Note:** These commands only need to be run once per reference genome. The indices can be reused for all future samples.

## For whole exome sequencing

If you're working with exome data, you'll also need an interval list defining your capture regions:

```bash
INTERVALS=~/references/exome_targets.interval_list
```

Always use the interval list that matches your capture kit. Contact your sequencing provider if you don't have this file.

## Sample data

For this tutorial, we'll use paired-end Illumina reads:

```bash
SAMPLE=NA12878  # Sample identifier
R1=~/data/${SAMPLE}_R1.fastq.gz
R2=~/data/${SAMPLE}_R2.fastq.gz
THREADS=8  # Adjust based on your system
```

# The workflow

## Step 1: Read alignment with BWA-MEM

We'll use BWA-MEM (or the faster BWA-MEM2) to align our reads to the reference genome. This step is critical - poor alignment will lead to poor variant calls.

**Why BWA-MEM?** It's accurate, fast, and handles reads from 70bp to several megabases. It's also the aligner used to develop GATK Best Practices, so the downstream tools are optimized for BWA output.

```bash
bwa mem -t ${THREADS} \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
  -K 100000000 \
  -Y \
  ${REF} ${R1} ${R2} | \
  samtools sort -@ ${THREADS} -o ${SAMPLE}.sorted.bam -
```

Expected runtime: 2-48 hours for 30X WGS (depends on CPU cores)

**What do these parameters mean?**
- `-t ${THREADS}`: Use multiple CPU cores for speed
- `-R "@RG\t..."`: Read group information (CRITICAL - see box below)
- `-K 100000000`: Process 100M bases per batch (affects performance/memory tradeoff)
- `-Y`: Use soft clipping for supplementary alignments (recommended for HaplotypeCaller)
- `| samtools sort`: Pipe directly to sorting to save disk space

### Understanding read groups

Read groups (`@RG`) are **essential** for GATK. They tell the variant caller how reads are related:

- `ID`: Read group identifier (must be unique)
- `SM`: Sample name (this is what appears in your VCF - **must be identical** across all read groups for the same sample)
- `PL`: Platform (ILLUMINA, PACBIO, etc.)
- `LB`: Library identifier (used to detect duplicate artifacts)
- `PU`: Platform unit (e.g., flowcell lane)

**Common mistake:** Using different `SM` values for different lanes of the same sample. This will cause GATK to treat them as different samples!

## Step 2: Mark duplicates

PCR duplicates arise during library preparation and sequencing. Marking (not removing) them prevents them from biasing variant calls.

```bash
gatk MarkDuplicatesSpark \
  -I ${SAMPLE}.sorted.bam \
  -O ${SAMPLE}.markdup.bam \
  -M ${SAMPLE}.markdup.metrics.txt \
  --create-index true
```

Expected runtime: 15-30 minutes for 30X WGS

The `MarkDuplicatesSpark` version is faster than regular `MarkDuplicates` on multi-core machines. The output file `${SAMPLE}.markdup.metrics.txt` contains useful QC information:

```bash
# Check your duplicate rate
grep "^LIBRARY" -A 1 ${SAMPLE}.markdup.metrics.txt
```

**What's a good duplicate rate?**
- PCR-free libraries: 5-10% duplicates
- Standard libraries: 10-30% duplicates
- Amplicon sequencing: Can be much higher

High duplication (>40% for WGS) suggests library prep issues or low complexity.

## Step 3: Base Quality Score Recalibration (BQSR)

Sequencers report base quality scores, but these can be systematically wrong due to various technical factors. BQSR corrects these scores using known variants as "truth."

**Why does this matter?** Variant callers rely heavily on base quality scores. Miscalibrated scores lead to false positives and false negatives.

### Generate the recalibration table

```bash
gatk BaseRecalibrator \
  -R ${REF} \
  -I ${SAMPLE}.markdup.bam \
  --known-sites ${DBSNP} \
  --known-sites ${MILLS} \
  --known-sites ${KNOWN_INDELS} \
  -O ${SAMPLE}.recal.table
```

For WES, add: `-L ${INTERVALS}`

Expected runtime: 15-45 minutes for 30X WGS

**How does BQSR work?** It identifies positions that match known variants in dbSNP, assumes these are correct, then models the difference between reported and expected quality scores based on machine cycle, sequence context, and other covariates.

### Apply the recalibration

```bash
gatk ApplyBQSR \
  -R ${REF} \
  -I ${SAMPLE}.markdup.bam \
  --bqsr-recal-file ${SAMPLE}.recal.table \
  -O ${SAMPLE}.bqsr.bam
```

Expected runtime: 20-40 minutes for 30X WGS

### Index the final BAM

```bash
samtools index ${SAMPLE}.bqsr.bam
```

Expected runtime: 2-5 minutes

**Optional QC:** Compare before/after BQSR plots:

```bash
gatk AnalyzeCovariates \
  -bqsr ${SAMPLE}.recal.table \
  -plots ${SAMPLE}.recal.pdf
```

This generates plots showing the recalibration effect. Look for systematic patterns that get corrected.

## Step 4: Variant calling with HaplotypeCaller

Now we're ready to call variants. HaplotypeCaller uses a local de novo assembly approach that's particularly good at calling indels.

### Call variants in GVCF mode (recommended)

```bash
gatk HaplotypeCaller \
  -R ${REF} \
  -I ${SAMPLE}.bqsr.bam \
  -O ${SAMPLE}.g.vcf.gz \
  -ERC GVCF \
  --native-pair-hmm-threads ${THREADS}
```

For WES, add: `-L ${INTERVALS}`

Expected runtime: 4-12 hours for 30X WGS

**What is GVCF mode (`-ERC GVCF`)?**

A GVCF (genomic VCF) contains information about **every position** in the genome, not just variant sites. This includes:
- Sites that appear variant in this sample
- Reference confidence at every other position

**Why use GVCF?** It enables joint genotyping. Even if you only have one sample now, you might want to add more samples later. Joint genotyping dramatically improves accuracy, especially for rare variants.

**What if I'm only analyzing one sample?**

You still need to genotype the GVCF:

```bash
gatk GenotypeGVCFs \
  -R ${REF} \
  -V ${SAMPLE}.g.vcf.gz \
  -O ${SAMPLE}.vcf.gz
```

Expected runtime: 1-3 hours for one WGS sample

## Step 5: Joint genotyping (for cohorts)

If you have multiple samples, joint genotyping is much more accurate than merging per-sample VCFs.

### Create a GenomicsDB workspace

GenomicsDB efficiently stores GVCF data for many samples:

```bash
# Create a sample map file
ls *.g.vcf.gz | awk '{print $1 "\t" $1}' > cohort.sample_map

gatk GenomicsDBImport \
  --sample-name-map cohort.sample_map \
  --genomicsdb-workspace-path cohort_db \
  -L chr1
```

**Note:** GenomicsDBImport must be run per-chromosome or per-interval. For whole genomes, run this 24 times (once per chromosome) or use scatter-gather parallelization.

### Genotype the cohort

```bash
gatk GenotypeGVCFs \
  -R ${REF} \
  -V gendb://cohort_db \
  -O cohort.chr1.vcf.gz
```

Expected runtime: 2-6 hours per chromosome for 100 samples

Then merge the per-chromosome VCFs:

```bash
bcftools concat -Oz -o cohort.vcf.gz cohort.chr*.vcf.gz
bcftools index -t cohort.vcf.gz
```

## Step 6: Variant quality filtering

Raw variant calls contain many false positives. We need to filter them.

### Option A: VQSR (recommended for cohorts ≥30 samples)

Variant Quality Score Recalibration (VQSR) uses machine learning to separate true variants from artifacts. It requires enough variants to train the model.

**When to use VQSR:**
- ≥30 whole genomes, or
- ≥50-100 whole exomes

**VQSR for SNPs:**

```bash
gatk VariantRecalibrator \
  -R ${REF} \
  -V cohort.vcf.gz \
  -mode SNP \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
    ~/references/hapmap_3.3.hg38.vcf.gz \
  -resource:omni,known=false,training=true,truth=false,prior=12.0 \
    ~/references/1000G_omni2.5.hg38.vcf.gz \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 \
    ~/references/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
    ${DBSNP} \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  -O cohort.snp.recal \
  --tranches-file cohort.snp.tranches \
  --rscript-file cohort.snp.plots.R
```

**What do these annotations mean?**
- `QD` (QualByDepth): Variant quality normalized by depth (higher is better)
- `MQ` (MappingQuality): How well reads map (60 is perfect, <40 is concerning)
- `MQRankSum`: Difference in mapping quality between ref and alt reads
- `ReadPosRankSum`: Are variants near the ends of reads? (strand bias indicator)
- `FS` (FisherStrand): Another strand bias metric (higher = more bias)
- `SOR` (StrandOddsRatio): Strand bias estimated using a contingency table

**Apply the SNP recalibration:**

```bash
gatk ApplyVQSR \
  -R ${REF} \
  -V cohort.vcf.gz \
  -mode SNP \
  --recal-file cohort.snp.recal \
  --tranches-file cohort.snp.tranches \
  --truth-sensitivity-filter-level 99.0 \
  -O cohort.snp.recalibrated.vcf.gz
```

The `--truth-sensitivity-filter-level 99.0` means "keep variants until you've kept 99% of true positives (as determined by the training sets)."

**Repeat for indels:**

```bash
gatk VariantRecalibrator \
  -R ${REF} \
  -V cohort.snp.recalibrated.vcf.gz \
  -mode INDEL \
  -resource:mills,known=false,training=true,truth=true,prior=12.0 \
    ${MILLS} \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
    ${DBSNP} \
  -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
  -O cohort.indel.recal \
  --tranches-file cohort.indel.tranches \
  --rscript-file cohort.indel.plots.R

gatk ApplyVQSR \
  -R ${REF} \
  -V cohort.snp.recalibrated.vcf.gz \
  -mode INDEL \
  --recal-file cohort.indel.recal \
  --tranches-file cohort.indel.tranches \
  --truth-sensitivity-filter-level 99.0 \
  -O cohort.filtered.vcf.gz
```

### Option B: Hard filters (for small cohorts or single samples)

When you don't have enough samples for VQSR, use hard filters. These are simple cutoffs based on GATK recommendations:

```bash
# Filter SNPs
gatk VariantFiltration \
  -R ${REF} \
  -V ${SAMPLE}.vcf.gz \
  -O ${SAMPLE}.snp.filtered.vcf.gz \
  --filter-name "QD2" --filter-expression "QD < 2.0" \
  --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
  --filter-name "FS60" --filter-expression "FS > 60.0" \
  --filter-name "MQ40" --filter-expression "MQ < 40.0" \
  --filter-name "MQRankSum-12.5" --filter-expression "MQRankSum < -12.5" \
  --filter-name "ReadPosRankSum-8" --filter-expression "ReadPosRankSum < -8.0" \
  --filter-name "SOR3" --filter-expression "SOR > 3.0"

# Filter indels
gatk VariantFiltration \
  -R ${REF} \
  -V ${SAMPLE}.vcf.gz \
  -O ${SAMPLE}.indel.filtered.vcf.gz \
  --filter-name "QD2" --filter-expression "QD < 2.0" \
  --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
  --filter-name "FS200" --filter-expression "FS > 200.0" \
  --filter-name "ReadPosRankSum-20" --filter-expression "ReadPosRankSum < -20.0"
```

**Important:** These filters mark variants as filtered but don't remove them. To extract only PASS variants:

```bash
gatk SelectVariants \
  -R ${REF} \
  -V ${SAMPLE}.snp.filtered.vcf.gz \
  --exclude-filtered \
  -O ${SAMPLE}.pass.vcf.gz
```

## Step 7: Quality control

Before analyzing your variants, check that everything looks reasonable.

### Count variants

```bash
bcftools stats ${SAMPLE}.pass.vcf.gz | grep "^SN"
```

Look for lines like:
```
SN    0    number of SNPs:    3245678
SN    0    number of indels:  445632
```

**Expected values:**
- WGS: 3-4M SNPs, 400-600k indels
- WES: 25-35k SNPs, 2-4k indels

Much more or less suggests problems.

### Check Ti/Tv ratio

Transition/transversion ratio is a key quality metric:

```bash
bcftools stats ${SAMPLE}.pass.vcf.gz | grep "^TSTV"
```

**Expected values:**
- WGS: 2.0-2.1
- WES: 2.8-3.2 (higher because coding regions are more constrained)

Low Ti/Tv (<1.8) suggests many false positives. Very high Ti/Tv (>3.5 for WGS) might indicate over-filtering.

### Check het/hom ratio

```bash
bcftools stats ${SAMPLE}.pass.vcf.gz | grep "^PSC"
```

The het/hom ratio should be ~1.5-2.0 for most populations. Consanguineous samples will have more homozygous variants (lower ratio).

### Depth distribution

Check if your sequencing depth matches expectations:

```bash
bcftools query -f '%DP\n' ${SAMPLE}.pass.vcf.gz | \
  awk '{sum+=$1; n++} END {print "Mean DP:", sum/n}'
```

This should match your target coverage (e.g., ~30 for 30X WGS).

### Sample-level metrics with Picard

```bash
# For WGS
gatk CollectWgsMetrics \
  -I ${SAMPLE}.bqsr.bam \
  -O ${SAMPLE}.wgs.metrics.txt \
  -R ${REF}

# For WES
gatk CollectHsMetrics \
  -I ${SAMPLE}.bqsr.bam \
  -O ${SAMPLE}.hs.metrics.txt \
  -R ${REF} \
  -BAIT_INTERVALS ${INTERVALS} \
  -TARGET_INTERVALS ${INTERVALS}
```

Key metrics to check:
- `MEAN_COVERAGE`: Should match your target
- `PCT_10X` (WGS) / `PCT_TARGET_BASES_10X` (WES): Percentage covered at 10X (should be >90% for good data)
- `MEDIAN_COVERAGE` should be close to `MEAN_COVERAGE` (if not, coverage is uneven)
- `ZERO_CVG_TARGETS_PCT` (WES): Percentage of targets with no coverage (should be <5%)

# Common pitfalls and solutions

## Reference genome mismatches

**Problem:** Using dbSNP from hg19 with a hg38 reference, or mixing GRCh37 and GRCh38.

**Symptoms:**
- BQSR fails or produces bizarre results
- Very high or very low variant counts
- Variants at nonsensical positions

**Solution:** Download all resources from the same source (Broad Resource Bundle) for the same build. Double-check chromosome naming (chr1 vs 1).

## Incorrect read groups

**Problem:** Wrong or missing `SM` tags; different `SM` values for the same sample.

**Symptoms:**
- Joint genotyping fails
- Sample appears multiple times in VCF
- Duplicate marking produces strange results

**Solution:** Check your read groups with:

```bash
samtools view -H ${SAMPLE}.bam | grep "^@RG"
```

All read groups for one sample must have identical `SM` tags.

## Skipping GVCF mode

**Problem:** Calling variants directly without GVCF intermediate.

**Symptoms:**
- Can't do joint genotyping later
- Lower accuracy for rare variants
- Difficult to combine with additional samples

**Solution:** Always use `-ERC GVCF` with HaplotypeCaller, even for single samples.

## Underpowered BQSR

**Problem:** Running BQSR on very small panels or amplicons with few variants.

**Symptoms:**
- BQSR warnings about insufficient data
- No improvement or worsening after BQSR

**Solution:** For small targeted panels (<100 genes), consider skipping BQSR or using a larger known-sites resource.

## WES interval handling

**Problem:** Not using the correct capture intervals, or not including padding.

**Symptoms:**
- Many targets have zero coverage
- Poor indel calling near exon boundaries

**Solution:**
- Get the correct interval list from your sequencing provider
- Add 100bp padding to intervals:

```bash
gatk PreprocessIntervals \
  -R ${REF} \
  -L ${INTERVALS} \
  --padding 100 \
  -O ${INTERVALS}.padded
```

## Sex chromosome considerations

**Problem:** Calling X and Y chromosomes with the same ploidy as autosomes.

**Symptoms:**
- Males called as heterozygous for X and Y variants (should be hemizygous)
- Incorrect genotype likelihoods

**Solution:** For male samples, call non-PAR regions of X and Y separately:

```bash
gatk HaplotypeCaller \
  -R ${REF} \
  -I ${SAMPLE}.bqsr.bam \
  -O ${SAMPLE}.chrX.g.vcf.gz \
  -ERC GVCF \
  -L chrX -L chrY \
  --sample-ploidy 1
```

Then merge with autosomal calls.

## Over-filtering or under-filtering

**Problem:** Removing too many real variants, or keeping too many false positives.

**Symptoms:**
- Ti/Tv ratio way off expected values
- Missing known pathogenic variants
- Absurd variant counts

**Solution:**
- For VQSR, try different sensitivity levels (95.0-99.5)
- For hard filters, tune based on your data's characteristics
- Validate against a truth set (e.g., GIAB) if available

# Rules of thumb

**Coverage:**
- WGS: 30X minimum (clinical often uses 40X)
- WES: 50-100X for germline, 100-200X for clinical
- Aim for uniformity: 80% of bases at >20X is better than 50% at >50X

**Variant counts (human germline):**
- WGS SNVs: 3-4 million
- WGS indels: 400-600k
- WES SNVs: 25-35k
- WES indels: 2-4k

**Quality metrics:**
- Ti/Tv: WGS ~2.0-2.1, WES ~2.8-3.2
- Het/hom ratio: 1.5-2.0
- Indel:SNV ratio (WES): ~0.1-0.2

**Filtering:**
- VQSR needs ≥30 WGS or ≥50-100 WES samples
- Hard filters work for smaller cohorts but require more manual tuning
- Always check a known truth set or validated variants after filtering

# Tips for production pipelines

1. **Parallelize by chromosome:** Split your calling by chromosome and merge at the end. This cuts runtime dramatically.

2. **Use Spark versions where available:** `MarkDuplicatesSpark` is much faster than regular `MarkDuplicates`.

3. **Pin tool versions:** GATK output can change between versions. Pin your version and record it:

```bash
gatk --version > pipeline.version.txt
```

4. **Keep intermediate files:** Disk is cheap. Keep your BAMs and GVCFs in case you need to reprocess.

5. **Document sample metadata:** Keep a spreadsheet linking sample IDs to sequencing runs, coverage, and QC metrics.

6. **Validate with IGV:** Spot-check variants in IGV (Integrative Genomics Viewer). Load your BAM and VCF and browse to a few variants. Do they look real?

7. **Check for contamination:** Use `VerifyBamID2` or `ContEst` to check for sample swaps or cross-contamination.

8. **Use sex checks:** Infer sex from X/Y coverage and compare to reported sex:

```bash
samtools idxstats ${SAMPLE}.bam | awk '$1 ~ /chrX|chrY/ {print $1, $3}'
```

9. **Mendelian consistency:** If you have trios, check for Mendelian violations:

```bash
gatk PhaseByTransmission \
  -V family.vcf.gz \
  -ped family.ped \
  -O family.phased.vcf.gz
```

10. **Compare to known samples:** GIAB (Genome in a Bottle) provides truth sets for NA12878 and other samples. If you sequence these, compare your calls:

```bash
rtg vcfeval \
  -b giab_truth.vcf.gz \
  -c ${SAMPLE}.vcf.gz \
  -t hg38.sdf \
  -o vcfeval_results
```

# Where to go from here

You now have a filtered VCF file. Next steps typically include:

1. **Annotation:** Use VEP, ANNOVAR, or Funcotator to annotate variants with gene names, predicted effects, and population frequencies.

2. **Filtering for analysis:** Extract coding variants, rare variants, or variants in specific genes depending on your research question.

3. **Trio analysis:** If you have parent-child trios, look for de novo variants or perform linkage analysis.

4. **Phenotype association:** For case-control studies, test for associations between variants and phenotypes.

5. **Clinical interpretation:** For clinical sequencing, classify variants using ACMG guidelines.

# Additional resources

- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
- [GATK Forums](https://gatk.broadinstitute.org/hc/en-us/community/topics) (excellent for troubleshooting)
- [Broad Resource Bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)
- [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) (truth sets for validation)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) (clinical variant database)
