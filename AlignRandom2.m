function [ average max min h stdev MAX] = AlignRandom2( N, FP, normalize_y_n )
%ALIGNRANDOM Summary of this function goes here
%   N = number of times to run
%   FP = sequence to randomize over
%   normalize_y_n = normalize, yes/no?

t=FP;
[a b]=size(t);
stdev=0;
max=0;
min=1000;
total=0;
h=zeros(N,1);
MAX=swalign(t,t,'Alphabet','NT');
for i=1:N
rs=getrandomseq(b);
[score align position]=swalign(rs,t,'Alphabet','NT'); %max = 97
if(score>max)
    max=score;
end
if(score<min)
    min=score;
end
total=total+score;
h(i,1)=score;
end
if(normalize_y_n==0)
    average=total/N;
    max=max;
    min=min;
else
    average=total/N/MAX;
    max=max/MAX;
    min=min/MAX;
    h=h/MAX;    
end

%count=0;
%for i=1:N
%    count=count+(average-h(i,1))^2;
%end
%stdev=sqrt(count/N);
stdev=std(h)

[average, max, min, stdev, MAX]
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