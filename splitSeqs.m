function [ positions, lengths ] = splitSeqs( seqs )
%SPLITSEQS Summary of this function goes here
%   Detailed explanation goes here
positions = [];
lengths = [];

x = multialign(seqs,'terminalGapAdjust',false);

[m,n] = size(seqs)

seqs3 = [];
for i=1:(n)
    s = x(i,:);
    seqs3 = [seqs3; strcat(s(25:50),s(70:95))];
end 
x3 = multialign(seqs3,'terminalGapAdjust',false);

seqs = seqs3;

threshold = 0.4;

[~,cell_count] = size(seqs);
distances = seqpdist(seqs(1:(cell_count)),'method','jukes-cantor','indels','pair');

AA = zeros(cell_count);
for i=1:cell_count
    for j=1:cell_count
        if(i>j)
            AA(i,j)=distances((j-1)*(cell_count-j/2)+i-j);
        else
            if(i~=j)
                AA(i,j)=distances((i-1)*(cell_count-i/2)+j-i);
            end
        end
    end
end

AA
positions = [positions;1];
lengths = [lengths;1];
count = 1;
cell_count
for i = 1 : (cell_count - 1)
    j = i + 1;
   AA(i,j)
   
    if(AA(i,j) >= threshold)
        positions = [positions;j];
        lengths = [lengths; 1];
        count = count + 1;
    else
        lengths(count) = lengths(count) + 1;
    end
end


end

