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
 verkko -d test --hifi hifi.fasta.gz --nano ont.fasta.gz --hic1 hic.R1.fastq.gz --hic2 hic.R2.fastq.gz --snakeopts --dry-run --screen human|more
```

<details><summary><b>What if I have multiple input files?</b></summary>
Verkko will take arbitary lists of inputs for each parameter so wildcards are ok (<code>ont*.fastq.gz</code> for example). Only one caveat, the Hi-C pairs have to be sorted in the same order to maintain read pairing (that is if you give <code>file1_R1.fastq.gz file2_R2.fastq.gz</code> to --hic1 you cannot give <code>file2_R2.fastq.gz file1_R1.fastq.gz1</code> to --hic2).
</details>

Now, let's run it for real (it will take 10-15 min)
```bash
 verkko -d test --hifi hifi.fasta.gz --nano ont.fasta.gz --hic1 hic.R1.fastq.gz --hic2 hic.R2.fastq.gz --no-correction --screen human
```

<b>Note: we would normally not use `--no-correction` but we want it to run faster for this tutorial.</b> While waiting, we can go through the [presentation](verkko.pptx) explaining how verkko works, motivation behind the steps, and some things that can go wrong.

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

<b>Let's look at these in detail:</b>
```bash
seqtk comp test/assembly.fasta 
haplotype1-0000001      4013576 1031011 993122  969324  1020119 0       0       0       163738  0       0       0
haplotype2-0000002      3991666 1027822 984853  964395  1014596 0       0       0       162904  0       0       0
```

We have two assembled sequences, of approximately 4 Mbp for each haplotype. Let's translate them to the graph:
```bash
cat test/assembly.scfmap
path haplotype1-0000001 haplotype1_from_utig4-0
piece000001
end
path haplotype2-0000002 haplotype2_from_utig4-3
piece000002
end
```

We want to find out the path for `haplotype1-0000001` so we will search for `haplotype1_from_utig4-0` in the `assembly.paths.tsv`:
```bash
grep -w haplotype1_from_utig4-0 test/assembly.paths.tsv 
haplotype1_from_utig4-0 utig4-1-,utig4-0+,utig4-4+      HAPLOTYPE1
```
and we can also check `haplotype2_from_utig4-3`
```bash
grep -w haplotype2_from_utig4-2 test/assembly.paths.tsv
haplotype2_from_utig4-3 utig4-2-,utig4-3+,utig4-5+      HAPLOTYPE2
```
<details><summary><b>Verkko assembly graph</b></summary>
<img src="graph.jpg" alt="verkko bandage graph" /><br>
<figcaption><em>The two paths each use either the red (haplotype 1) or the blue (haplotype2) node. The other large gray nodes are ambiguous and can be randomly assigned a haplotype. Homozygous nodes would also be gray but would have higher coverage, approximately 2x, relative to red/blue).</em></figcaption>
</details>

#### A few helper scripts:
There are common things we do with verkko assemblies, such as alignment to a reference (if one exists) and looking for T2T contigs/scaffolds (telomeres on both end and a gap or no gaps). These scripts are available at the [MarBL training GitHub](https://github.com/marbl/training/tree/main/part2-assemble/docker/marbl_utils) and are also conveniently included in our GitPod. Let's run them:
```bash
cd test
# this just requires an assembly fasta file and generates assembly.t2t_ctgs, assembly.t2t_scfs, assembly.telomere.bed, assembly.gaps.bed
bash /workspace/marbl_utils/asm_validation/getT2T.sh assembly.fasta
cat assembly.telomere.bed 
haplotype1-0000001      4002154 4013576 4013576
haplotype2-0000002      3985834 3987666 3987666
haplotype2-0000002      3987834 3991666 3991666
```
We have telomeres on the ends of both paths, this means utig4-[45] have telomeres but we don't have to parse that, we can add them to the graph automatically:
```bash
/workspace/??/remove_nodes_add_telomere.py --telo assembly.telomere.bed 
```

We can also assign the nodes to chromosomes if we have a previous reference available.
```bash
# this take a reference which will be HPC-compressed if there isn't an HPC version already, an identity (default 99), and the assembly to align
# it outputs assembly.mashmap.out, translation_hap1, translation_hap2, and assembly.homopolymer-compressed.chr.csv
bash /workspace/marbl_utils/asm_validation/getChrNames.sh /workspace/reference.fasta 99 assembly.fasta

cat translation_hap[12]
haplotype1-0000001      chr12   4013576 133324548
haplotype2-0000002      chr12   3991666 133324548

cat assembly.mashmap.out
haplotype1-0000001      4013576 0       2320000 +       chr12   133324548       129319114       131624898       18      2320000 26      id:f:0.99773    kc:f:1.19639
haplotype1-0000001      4013576 2330000 3710000 +       chr12   133324548       131621174       133001134       20      1380000 28      id:f:0.99857    kc:f:1.23998
haplotype1-0000001      4013576 3710000 4010000 +       chr12   133324548       133027651       133324459       19      300000  30      id:f:0.998912   kc:f:0.869261
haplotype2-0000002      3991666 0       3990000 +       chr12   133324548       129324978       133323934       20      3998956 28      id:f:0.99856    kc:f:1.16405
cd ..
```

<details><summary><b>Verkko assembly graph with telomeres and chr names</b></summary>
<img src="graph_tel.jpg" alt="verkko bandage graph" /><br>Same region as above but now we have added telomeric nodes to the graph (indicated in thick green). We also have labeled the nodes by their chromosome assignment based on thereference. This region is apparently from one end of Chr 12.</em></figcaption>
</details>

#### Editing an assembly
Lastly, let's say we decide the phasing is incorrect and we think utig4-4 should be in haplotype2 not 1 like it is now. We can edit the paths file:
```bash
cp test/8-hicPipeline/rukki.paths.gaf ./updated.gaf
vi updated.gaf
```
<details><summary><b>edited paths</b></summary>
<code>
name    path    assignment
haplotype1_from_utig4-0 <utig4-1>utig4-0        HAPLOTYPE1
haplotype2_from_utig4-3 <utig4-2>utig4-3>utig4-4        HAPLOTYPE2
na_unused_utig4-6       >utig4-6        NA
na_unused_utig4-7       >utig4-7        NA
na_unused_utig4-8       >utig4-8        NA
haplotype1_from_utig4-5 >utig4-5        HAPLOTYPE1
</code>
</details>

 Now that we have updated the paths, we can ask verkko to give us new consensus for these:
```bash
verkko -d cns --hifi chr12/hifi.fasta.gz --nano chr12/ont.fasta.gz --hic1 chr12/hic.R1.fastq.gz --hic2 chr12/hic.R2.fastq.gz --local --paths updated.gaf --assembly test
seqtk comp cns/assembly.fasta

haplotype1-0000001      3713519 946361  931851  903311  931996  0       0       0       156668  0       0       0
haplotype1-0000002      309933  87492   64610   67299   90532   0       0       0       7642    0       0       0
haplotype2-0000003      3995292 1028429 984872  966206  1015785 0       0       0       162910  0       0       0
```

There's lots more to learn about editing/finishing an assembly, including resolving remaining tangles, filling gaps, and patching with other assemblies. If you're interested in that, start by taking a look at the [MarBL training GitHub](https://github.com/marbl/training).
