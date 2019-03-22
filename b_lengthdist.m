function [output]=lengthdist(seqarr)
[a b]=size(seqarr);
output=zeros(1,b);
for i=1:b
    output(1,i)=length(char(seqarr(1,i)));
end
end