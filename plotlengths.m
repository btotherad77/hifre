function [output]=plotlengths(seqarr,rows,columns)
[a b]=size(seqarr)
for q=1:a
    output=zeros(1,b);
    for i=1:b
        output(1,i)=length(char(seqarr(1,i)));
    end
    subplot(rows,columns,q)
    histogram(output)
    title(strcat('plot #',num2str(q)))
end
end