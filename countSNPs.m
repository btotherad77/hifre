function [ SNPcount ] = countSNPs( output )
%COUNTSNPS Summary of this function goes here
%   Detailed explanation goes here
SNPs = output(:,17);
IDs = num2str(cell2mat(output(:,8)));
[m,~]=size(output);

SNPcount=containers.Map;

for i = 1:m
    
    x=char(strcat(IDs(i,:),SNPs(i)));
    x = x(1:5);
    if(SNPcount.isKey(x))
        SNPcount(x)=SNPcount(x)+1;
    else
        SNPcount(x)=1;
    end
    
end

end

