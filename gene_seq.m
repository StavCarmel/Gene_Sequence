function [UTR5,ORF,UTR3,gene]= gene_seq(indices, strand, prev_gene, next_gene, chromseq)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% this function takes out the gene's sequence, concluding the 5'UTR, ORF, 3'UTR out of the genome sequence
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%indices   - vector of numbers indicates the gene's ORF indices from the genome sequence.
%%strand    - num 1 / -1, 1 indicates the gene is taken out from the primary strand, -1 indicates
% it's from the complementary.
%%prev_gene - the distance in nucleotides from the end of the stop codon of the previous gene to the start
% codon of the gene in the genome sequence in order to take the rightamount of nt for the 5'UTR.
%%next_gene - the distance in nucleotides from the end of the stop codon to the start codon of the next
% gene in the genome sequence in order to take the right amount of nt for the 5'UTR.
%%chromseq  - the genomes sequence.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ORF='';
if strand==1 %from primary strand
	cur_indeces=sort(indices);
    UTR5=chromseq(cur_indeces(1)-prev_gene+1:cur_indeces(1)-1);
    UTR3=chromseq(cur_indeces(end)+1:cur_indeces(end)+next_gene);
    for i=1:2:length(indices)-1
        ORF=[ORF,chromseq(cur_indeces(i):cur_indeces(i+1))];
    end
else % from complementary
    cur_indeces= fliplr(sort(indices));
    UTR5flip=chromseq(cur_indeces(1)+1:cur_indeces(1)+prev_gene);
    UTR5=seqrcomplement(UTR5flip);
    UTR3flip=chromseq(cur_indeces(end)-next_gene+1:cur_indeces(end)-1);
    UTR3=seqrcomplement(UTR3flip);
    for i=1:2:length(indices)-1
        cur_CDS=chromseq(cur_indeces(i+1):cur_indeces(i));
        ORF=[ORF,seqrcomplement(cur_CDS)];
    end
end
gene=[UTR5,ORF,UTR3];
end

