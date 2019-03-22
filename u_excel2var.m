function [ output ] = u_excel2var( folder )
%U_EXCEL2VAR Summary of this function goes here
%   Detailed explanation goes here
                        %1: passing sequence #
                        %2: sequence # (in original file)
                        %3: number of repeats
                        %4: average of individual alignments
                        %5: consensus alignment
                        %6: assigned barcode
                        %7: actual barcode
                        %8: barcode #
                        %9: quality score values
                        %10: average probability of base being correct
                        %11: ID (1-oligo, 2-GFP)
                        %12: rev_flag
                        %13: average base individual alignment
                        %14: average base consensus alignment
                        %15: avg swalignment
                        %16: sw_score post alignment
                        %17: SNPs
Files=dir(strcat(folder,'/*.*'));
file1 = strcat(folder,'/',Files(1).name)
file1 = strcat(folder,'/',Files(2).name)
file1 = strcat(folder,'/',Files(3).name)
output = table2cell(readtable(file1));
for k = 4:length(Files)
   %output = [output;readtable(strcat(folder,'/',Files(k).name),{'Passing #','Original #','NumRepeats','AvgNWScore','PostAlnNWScore','Consensus','AssBC','ActBC','BC#','QSV','AvgQ','ID','rev_flag','avgBase','ConsensusBase','AvgSWScore','PostAlnSWScore','SNPs'})];
   Files(k).name
   output = [output;table2cell(readtable(strcat(folder,'/',Files(k).name)))];
end

end

