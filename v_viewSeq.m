function [ output_args ] = v_viewSeq( seqs )
%V_VIEWSEQ Summary of this function goes here
%   Detailed explanation goes here
[~,n]=size(seqs);
%x = multialign(seqs,'terminalGapAdjust',false);
%seqalignviewer(x)
x = multialign(seqs,'terminalGapAdjust',true);
seqalignviewer(x)


[s1,a1,n1] = swalign(x(n,:),'AGAGACAGNNNNNNNNNNGTCTACATCC');
[s2,a2,n2] = swalign(x(n,:),'TGTAGATANNNNNNNNNNAAGGTTATGC');

[~,len1]=size(a1);
len1=len1-18;
start1=n1(1)+8;
[~,len2]=size(a2);
len2=len2-18;
start2=n2(1)+8;

seqs2=[];

for i=1:(n)
    s = x(i,:);
    seqs2 = [seqs2; strcat(s(start1:(start1+len1)),s(start2:(start2+len2)))];
end 
x2 = multialign(seqs2,'terminalGapAdjust',false,'GAPOPEN',10);
%seqalignviewer(x2)

seqs3 = [];
for i=1:(n)
    s = x(i,:);
    seqs3 = [seqs3; strcat(s(25:50),s(70:95))];
end 
x3 = multialign(seqs3,'terminalGapAdjust',false);
%seqalignviewer(x3)
%generateTree(seqs3)

output_args = [start1, len1,start2, len2];

end

