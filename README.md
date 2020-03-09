# CaReAl

![CaReAl](/imgs/CaReAl_logo.png)

 * CaReAl (Capturing Read Alignments) is a high-performance alignment capturing tool for visualizing the read-alignment status of nucleotide sequences and associated genome features.

 * Features
   * Visualizing full-depth of aligned reads.
   * Capturing read-level data rapidly and conveniently.
   * Optimized for the systematic exploration of regions of interest.
   * Useful for evaluating variant calls and detecting technical biases.

 * Overview  
 ![Overview](/imgs/CaReAl_overview.png)
 
 * How to interpret the plot?  
 ![CaReAl snapshot](/imgs/CaReAl_example.png)

 * Prerequisite
   * Python (version >= 2.7) with packages _tabix, multiprocessing, optparse_
   * R (version >= 3.3) with library _graphics_
   * SAMTools (version >= 0.1.19)
   
 * Usage
   - Synopsis  
   ```./careal -b <BAM or BAMs in TXT> -t <target or targets in TXT> -r <reference genome> [OPTIONS]```
   * Required fields  
        **-b, --bamfile \<FILE\>**: input a BAM or BAMs with directory paths in TXT, (i.e., sample_01.bam, samples.txt)  
        **-t, --target \<FILE\>**: chromosome and position or chromosome and range, or targets in TXT, (i.e., chr10:53933206, chr10:53933106-53933306, targets.txt)  
        **-r, --reference \<FILE\>**: reference genome file in FASTA, (i.e., genome.fasta)  
   * Optional fields  
        **-g, --gene \<FILE\>**: gene information file in BED indexed by Tabix  
        **-v, --with-vcf \<STRING\>**: query VCF \[TRUE or FALSE\], default is FALSE  
        **-o, --folder \<STRING\>**: output directory name, default is "OUTPUT"  
        **-n, --cpu \<NUMBER\>**: number of CPUs to use, default is 4  

 * Database
   * Provide gene symbols on a plot based on a given target.
   * It must be compressed by bgzip and indexed by tabix.
   * Gene information: [Homo_sapiens_hg19_75.txt.gz](/database/Homo_sapiens_hg19_75.txt.gz)
   * File structure (in BED)  
   ![BED format](/imgs/gene_db_structure.png)

 * FAQ
   * What are good targets for reads investigation and how do I get them?  
     * Variants with low depth coverage (i.e., DP < 5 for heterozygous SNV) are good targets and you can get targets using [VCFTools](https://vcftools.github.io/man_latest.html) (with '--minDP', '--maxDP' options).
     * [Minikel et al.](http://stm.sciencemag.org/content/8/322/322ra9) reviewed rare variants and provided IGV screenshots of the variants deemed to be genuine on [GitHub](https://github.com/ericminikel/prnp_penetrance/tree/master/supplement/igv).
   * What are the known error patterns of alignment?  
     * [Seo and Park et al.](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181304) reported four distinct error patterns from Ion Proton : (1) simplicity region, (2) SNV cluster, (3) peripheral sequence read, and (4) base inversion.
     * [Schirmer et al.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0976-y) evaluated variants from Illumina platforms and revealed a substantial bias related to motifs ending in "GG".
   * What is the maximum depth and width of a snapshot?
     * We tested a plot with a maximum depth of 2200x and width of 3,000 bp, but they depend on the graphics hardware.

 * Troubleshoot
   - BAM file should be indexed with SAMTools and please use the command below to build an index file:  
   ```$ samtools index Sample_001.bam```
   - To make sure that tview in samtools works properly, try the command below:  
   ```$ samtools tview -d T -p chr1:2409792-2409992 Sample_001.bam hg19.fa```
   - Make sure that the VCF is located in the same folder as the BAM file and that their filenames are identical.  
     For example, _Sample_001.bam, Sample_001.bai, Sample_001.vcf.gz,_ and _Sample_001.vcf.gz.tbi_ should be in the same folder.  
     ![File names](/imgs/bam_vcf_directory.png)
