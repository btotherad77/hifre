function [ output ] = f_countoutput( input )
%F_COUNTOUTPUT Summary of this function goes here
%   Detailed explanation goes here
k=1;
barcodes=cell2mat(input(:,7+k));
seqID=cell2mat(input(:,10+k));
SNPs=input(:,17);

output=zeros(14,6);
output(:,1)=0:13;
[x,y]=size(barcodes);
for i=1:x
    if(barcodes(i) == 99) %13 --> empty barcode in alignment
        tempbc = 13;
    else
        if(barcodes(i) == -1) %-1 --> unassignable barcode in alignment
            tempbc = 1;
        else
            tempbc = barcodes(i) + 1;
        end
    end
    
    output(tempbc,seqID(i)+1)=output(tempbc,seqID(i)+1)+1;
    
    if(strcmp(char(SNPs(i)),'GTC') == 1)
        output(tempbc,5)=output(tempbc,5)+1;
    end
    if(strcmp(char(SNPs(i)),'CAG') == 1)
        output(tempbc,4)=output(tempbc,4)+1;
    end
end

%Oligo/(GFP+GFPSNP)
output(:,6)=output(:,2)./(output(:,2)+output(:,3));
%GFPSNP/(GFP+GFPSNP)
output(:,7)=output(:,4)./(output(:,4)+output(:,5));
end

