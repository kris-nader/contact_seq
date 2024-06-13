## download th reference fasta and the GTF
wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz;
wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz;

## unzip
gzip -dk Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
gzip -dk Homo_sapiens.GRCh38.112.gtf.gz;

## copy a new file -- never modify the original
cp Homo_sapiens.GRCh38.dna.primary_assembly.fa hg38_fasta.fa
cp Homo_sapiens.GRCh38.112.gtf hg38_gtf.gtf

## append the fasta for the 2 genes
cat ./kristen_contact/GFP.txt >> hg38_fasta.fa
cat ./kristen_contact/mScarlet.txt  >> hg38_fasta.fa

## check if they have been appended to the fasta file-- there should be 2 new genes so 194+2
grep ">" hg38_fasta.fa
grep -c "^>" hg38_fasta.fa ## 196
grep -c "^>" Homo_sapiens.GRCh38.dna.primary_assembly.fa ## 194

## check length of the fasta genes
tail -n +2 ./kristen_contact/GFP.txt | wc | awk '{print $3-$1}' ## 717
tail -n +2 ./kristen_contact/mScarlet.txt | wc | awk '{print $3-$1}' ## 690

## check gtf format
echo -e "GFP\tcustom\tgene\t1\t717\t.\t+\t.\tgene_id \"id_GFP\"; gene_name \"GFP\"; gene_biotype \"protein_coding\";"
echo -e "GFP\tcustom\texon\t1\t717\t.\t+\t.\tgene_id \"id_GFP\"; gene_name \"GFP\"; gene_biotype \"protein_coding\";"
echo -e "mscarlet\tcustom\tgene\t1\t690\t.\t+\t.\tgene_id \"id_mscarlet\"; gene_name \"mscarlet\"; gene_biotype \"protein_coding\";"
echo -e "mscarlet\tcustom\texon\t1\t690\t.\t+\t.\tgene_id \"id_mscarlet\"; gene_name \"mscarlet\"; gene_biotype \"protein_coding\";"
## add in gtf format exon and gene
echo -e "GFP\tcustom\tgene\t1\t717\t.\t+\t.\tgene_id \"id_GFP\"; gene_name \"GFP\"; gene_biotype \"protein_coding\";"  >> hg38_gtf.gtf
echo -e "GFP\tcustom\texon\t1\t717\t.\t+\t.\tgene_id \"id_GFP\"; gene_name \"GFP\"; gene_biotype \"protein_coding\";"  >> hg38_gtf.gtf 
echo -e "mscarlet\tcustom\tgene\t1\t690\t.\t+\t.\tgene_id \"id_mscarlet\"; gene_name \"mscarlet\"; gene_biotype \"protein_coding\";" >> hg38_gtf.gtf
echo -e "mscarlet\tcustom\texon\t1\t690\t.\t+\t.\tgene_id \"id_mscarlet\"; gene_name \"mscarlet\"; gene_biotype \"protein_coding\";"  >> hg38_gtf.gtf 

## .fa and .gtf to .fa.gz and .gtf.gz
gzip -k hg38_gtf.gtf
gzip -k hg38_fasta.fa

## run mkref
split-pipe \
--mode mkref \
--genome_name hg38 \
--fasta hg38_fasta.fa.gz \
--genes hg38_gtf.gtf.gz \
--output_dir ./newgenome_no_genome_name/

## unzup and check if the 2 genes are there
gzip -dk ./newgenome_no_genome_name/process/genome.fas.gz
grep "GFP" ./newgenome_no_genome_name/process/genome.fas
grep "mScarlet" ./newgenome_no_genome_name/process/genome.fas



