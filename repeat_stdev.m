function [ output, output2 ] = repeat_stdev( input )
%REPEAT_STDEV Summary of this function goes here
%   Detailed explanation goes here
output = zeros(max(input(:,1)),6);
output2 = zeros(max(input(:,1)),18);
last=[];
last_i = 1;
count = 0;
for i = min(input(:,1)):5:max(input(:,1))
    %temp = input(input(:,1)==i,:);
    temp = [input(input(:,1)==i,:);input(input(:,1)==i+1,:);input(input(:,1)==i+2,:);input(input(:,1)==i+3,:);input(input(:,1)==i+4,:)];
    [r,~]=size(temp);
    output(i,1)=i;
    output(i,2)=r;
    output(i,3)=mean(temp(:,2));
    output(i,4)=std(temp(:,2));
    output(i,5)=mean(temp(:,3));
    output(i,6)=std(temp(:,3));
     if(count>0 && length(temp) > 3)
             [length(last),length(temp)]
             
%          if(length(temp)>length(last))
%              [h1, p1, ci1, stats1] = ttest2(temp(1:length(last),2),last(:,2));
%              [h2, p2, ci2, stats2] = ttest2(temp(1:length(last),3),last(:,3));
%              output2(i,:) = [i, last_i, h1, p1, ci1', h2, p2, ci2'];
%          else
%              [h1, p1, ci1, stats1] = ttest2(temp(:,2),last(1:length(temp),2));
%              [h2, p2, ci2, stats2] = ttest2(temp(:,3),last(1:length(temp),3));
%              output2(i,:) = [i, last_i, h1, p1, ci1', h2, p2, ci2'];
%          end
             [h1, p1, ci1, stats1] = ttest2(temp(:,2),last(:,2),'Vartype','unequal');
             [h2, p2, ci2, stats2] = ttest2(temp(:,3),last(:,3),'Vartype','unequal');
             output2(i,:) = [i, last_i, output(i,3:6), output(i-5,3:6), h1, p1, ci1', h2, p2, ci2'];
     end
     
     last = temp;
     last_i = i;
    count = count + 1;
end

end

