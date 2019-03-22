function [ output ] = getRandom( length, num )
%GETRANDOM Summary of this function goes here
%   Detailed explanation goes here
output=[];
for i=1:num
    output=[output;getrandomseq(length)];
end

end

function [s] = getrandomseq(n)
s='';
for i=1:n
    s=strcat(s,getbase(randi([0,3])));
end

end
function [base] = getbase(i)
base='X';
    if(i==0)
        base='A';
    else
        if(i==1)
            base='T';
        else
            if(i==2)
                base='C';
            else
                base='G';
            end
        end
    end
end