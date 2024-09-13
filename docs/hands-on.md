# Verkko tutorial

Get a "Regular" instance: 4 cores, 8 G RAM, 30 GB Storage

## Learning objects

 - Understand verkko steps
 - Run verkko
 - Inspect logging output as it runs
 - Practice common troubleshooting tips
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
 verkko -d asm --hifi hifi.fasta.gz --nano ont.fasta.gz --hic1 hic.R1.fastq.gz --hic2 hic.R2.fastq.gz --snakeopts --dry-run|more
```
