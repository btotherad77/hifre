function [  ] = y_writetoexcel( filename,data,start )
%Y_WRITETOEXCEL Summary of this function goes here
%   Detailed explanation goes here
[m,n]=size(data)
count=0;
s=50000;
NUM=floor(m/s)
for i=0:NUM
    [(i*s+1), ((i+1)*s)]
    if(((i+1)*s)>m)
        xlswrite(strcat('2018_05_15-Analysis(',num2str(i),')'),data((i*s+1):m,:))
    else
        xlswrite(strcat('2018_05_15-Analysis(',num2str(i),')'),data((i*s+1):((i+1)*s),:))
    end
end
end

