# Verkko tutorial

Get a "Regular" instance: 4 cores, 8 G RAM, 30 GB Storage

## Learning goals

 - Understand verkko steps
 - Run verkko
 - Inspect logging output as it runs
 - Get familiarized with assembly outputs

 ## Where's the data?

 - Included in the GitPod when launched
 - Can be downloaded from https://obj.umiacs.umd.edu/training/bga2024_verkko_t2t_data.tar.gz 

## Generating the assembly

Let's try an assembly of some test data. First, let's get familiar with verkko syntax
```bash
verkko --help
```
now let's give it some data and see what it will do
```bash
 verkko -d test --hifi hifi.fasta.gz --nano ont.fasta.gz --hic1 hic.R1.fastq.gz --hic2 hic.R2.fastq.gz --snakeopts --dry-run|more
```

<details><summary><b>What if I have multiple input files?</b></summary>
Verkko will take arbitary lists of inputs for each parameter so wildcards are ok (<code>ont*.fastq.gz</code> for example). Only one caveat, the Hi-C pairs have to be sorted in the same order to maintain read pairing (that is if you give <code>file1_R1.fastq.gz file2_R2.fastq.gz</code> to --hic1 you cannot give <code>file2_R2.fastq.gz file1_R1.fastq.gz1</code> to --hic2).
</details>

Now, let's run it for real (it will take 10-15 min)
```bash
 verkko -d test --hifi hifi.fasta.gz --nano ont.fasta.gz --hic1 hic.R1.fastq.gz --hic2 hic.R2.fastq.gz --no-correction
```

<b>Note: we would normally not use `--no-correction` but we want it to run faster for this tutorial.</b> While waiting, we can go through presentation explaining how verkko works and motivation behind the steps.

### Let’s take a look at the outputs
```bash
ls -h test/
0-correction     8-hicPipeline                              assembly.paths.tsv
1-buildGraph     assembly.colors.csv                        assembly.scfmap
2-processGraph   assembly.disconnected.fasta                assembly.unassigned.fasta
3-align          assembly.fasta                             emptyfile
3-alignTips      assembly.haplotype1.fasta                  hifi-corrected.fasta.gz
4-processONT     assembly.haplotype2.fasta                  snakemake.sh
5-untip          assembly.homopolymer-compressed.gfa        verkko.yml
6-layoutContigs  assembly.homopolymer-compressed.layout
7-consensus      assembly.homopolymer-compressed.noseq.gfa
```
#### Some of the important outputs:
<dl>
<dt>assembly.fasta</dt>
<dd>The main output, all sequences as fasta</dd>
<dt>assembly.haplotype*.fasta</dt>
<dd>Output if you have phasing information like trio or Hi-C (like we did here)</dd>
<dt>assembly.homopolymer-compressed[.noseq].gfa</dt>
<dd><a href="https://github.com/GFA-spec/GFA-spec">GFA</a> formatted assembly results, in HPC space, useful for visualing the assembly</dd>
<dt>assembly.scfmap</dt>
<dd>Translation between final names and the initial paths</dd>
<dt>assembly.paths.tsv</dt>
<dd>Translation between the path name and the graph nodes making them up</dd>
</dl>

Let's look at these in detail:
```bash
seqtk comp test/assembly.fasta 
haplotype1-0000001      4108446 1137367 923158  912595  1135326 0       0       0       93992   0     00
haplotype2-0000002      4106285 1134608 912192  922701  1136784 0       0       0       93936   0     00
```

We have two assembled sequences, of approximately 4.1 Mbp for each haplotype. Let's translate them to the graph:
```bash
cat test/assembly.scfmap
path haplotype1-0000001 haplotype1_from_utig4-10
piece000001
end
path haplotype2-0000002 haplotype2_from_utig4-2
piece000002
end
path unassigned-0000003 na_unused_utig4-0
piece000003
end
```

We want to find out the path for `haplotype1-0000001` so we will search for `haplotype1_from_utig4-10` in the `assembly.paths.tsv`:
```bash
grep haplotype1_from_utig4-10 test/assembly.paths.tsv 

haplotype1_from_utig4-10        utig4-8-,utig4-4-,utig4-3+,utig4-7+,utig4-9+,utig4-10+  HAPLOTYPE1
```
and we can also check `haplotype2_from_utig4-2`
```bash
grep haplotype2_from_utig4-2 test/assembly.paths.tsv

haplotype2_from_utig4-2 utig4-2+,utig4-9-,utig4-6-,utig4-3-,utig4-5+,utig4-8+   HAPLOTYPE2
```
<details><summary><b>Verkko assembly graph</b></summary>
<img src="graph.jpg" alt="verkko bandage graph" /><br>
<figcaption><em>The two paths each use either the red (haplotype 1) or the blue (haplotype2) node. The other large gray nodes are homozygous (node the higher coverage relative to red/blue). The small bubbles (e.g. <code>utig4-[45]</code>) have no signal but are short so they are randomly assigned a haplotype.</em></figcaption>
</details>

#### A few helper scripts:
There are common things we do with verkko assemblies, such as alignment to a reference (if one exists) and looking for T2T contigs/scaffolds (telomeres on both end and a gap or no gaps). These scripts are available at the [MarBL training GitHub](https://github.com/marbl/training/tree/main/part2-assemble/docker/marbl_utils) and are also conveniently included in our GitPod. Let's run them:
```bash
cd test
# this just requires an assembly fasta file and generates assembly.t2t_ctgs, assembly.t2t_scfs, assembly.telomere.bed, assembly.gaps.bed
bash /workspace/marbl_utils/asm_validation/getT2T.sh assembly.fasta

# this take a reference which will be HPC-compressed if there isn't an HPC version already, an identity (default 99), and the assembly to align
# it outputs assembly.mashmap.out, translation_hap1, translation_hap2, and assembly.homopolymer-compressed.chr.csv
bash /workspace/marbl_utils/asm_validation/getChrNames.sh /workspace/chm13v2.0.fasta 99 assembly.fasta
```
