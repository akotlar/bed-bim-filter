# bed-bim-filter
Takes a bed (or bim) file, and filters a .vcf file


Currently assumes that positions are 1-based, as in .bim file.

Can normalize the chromosome representation of the bim and vcf input files using ucscChr (will output as chrN).

Uses goroutines + channels to effectively use multiple cpu cores. Fast: can filter 22 million variants out of 40 million at approximately IO limits.
