# Some tips #

you may want to check these formats fa(fasta) vcf(Variant Call Format) gtf(The Gene transfer format)

### the sequence file .fa ###

* you can access any nucleotide by position, by dividing the position by the number of nucleotides per row, where the result is the row number, and the division remainder is the position in this row
* However, it's more efficient to use library that does it for you.. try seqIO

### the mutations file .vcf ###

* you will notice that it is a very big file.. to check it's contents you can use the command `head -100 mutations.vcf`	this will view the top 100 rows
* we don't care much now about the file header in this exercise
* it has mutations for all chromosomes, to reduce it's size you can only select mutations for chr22.. the command for this is `awk -v val='22' '$1 == val' mutations.vcf`
* you will also notice that some mutations are not point mutations.. meaning not a single nucleotide replaced by a single nucleotide. We need to filter these out, you can use this command `awk 'length($4) == 1 && length($5) == 1 && $1 == 22' mutations.vcf`
	* keep in mind you may still have a single nucleotide insertion or deletion after this step, you can handle this in you code.

### The annotations file .gtf ###

* it also has mutations for all chromosomes
* CDS means coding sequence, consider it the exon, the only difference is that an exon may have some non-coding sequence before or after
* one gene may have more than one exonic coding sequence, and it has a direction.. be careful the sequence in the fa file is in the forward strand (+ve strand), this is why regions coded in the reverse strand (-ve) will have a start position after the end position
* the exon(CDS) may end in the middle of a codon, the rest of the codon is in the following exon

### other tips ###

* To get the corresponding amino acid to the codon, you can use the attached json file. However, it's recommended to use a package, try seqIO
* To get the mass of an amino acid, you can create a similar json file.. amino acid masses should be easy to find.