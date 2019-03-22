function [ report ] = examineOutput( output , figure_start)
%EXAMINEOUTPUT Summary of this function goes here
%   Detailed explanation goes here

    %1: passing sequence #
    %2: sequence # (in original file)
    %3: number of repeats
    %4: average of individual alignments
    %5: consensus alignment
    %6: assigned barcode
    %7: actual barcode
    %8: barcode #
    %9: quality score values
    %10: average probability of base being correct
    %11: ID (oligo 1 or oligo 2)
    %12: rev_flag
    %13: average base individual alignment
    %14: average base consensus alignment
    %15: avg swalignment
    %16: sw_score post alignment
    %17: SNPs

k=1;
c1 = [0.5,0.5,0.5];
c2 = [0.9,0,0];
subplot2 = @(m,n,p) subtightplot (m, n, p, [0.01 0.02], [0.06 0.06], [0.06 0.03]);
%% SW alignment improvement vs length
figure(8+figure_start);
subplot(2,1,1)
p1 = scatter(cell2mat(output(:,3)),cell2mat(output(:,15+k))/max(cell2mat(output(:,15+k))),'o','MarkerEdgeColor','k','MarkerFaceColor',c2);
%hold on
subplot(2,1,2)
p2 = scatter(cell2mat(output(:,3)),cell2mat(output(:,14+k))/max(cell2mat(output(:,15+k))),'o','MarkerEdgeColor','k','MarkerFaceColor',c1);
%p2.MarkerFaceAlpha = .5;
%p2.MarkerEdgeAlpha = .5;

ylim([-0.1,1.01])
hold off
title('swalignment improvement vs. length')

figure(9+figure_start);
p1 = scatter(cell2mat(output(:,3)),cell2mat(output(:,15+k))/max(cell2mat(output(:,15+k))),'o','MarkerEdgeColor','k','MarkerFaceColor',c2);
hold on
p2 = scatter(cell2mat(output(:,3)),cell2mat(output(:,14+k))/max(cell2mat(output(:,15+k))),'o','MarkerEdgeColor','k','MarkerFaceColor',c1);
title('swalignment improvement vs. length')




%% ID-assigned alignment improvement
figure(1+figure_start);
h1 = histogram(cell2mat(output(:,3+k)));
hold on
h2 = histogram(cell2mat(output(:,4+k)));
h1.Normalization = 'probability';
h1.BinWidth = 5;
h2.Normalization = 'probability';
h2.BinWidth = 5;
xlim([0 160]);
hold off

title('alignment improvement')

%% base alignment improvement
figure(2+figure_start)
h3 = histogram(cell2mat(output(:,12+k)));
hold on
h4 = histogram(cell2mat(output(:,13+k)));
h3.Normalization = 'probability';
h3.BinWidth = 5;
h4.Normalization = 'probability';
h4.BinWidth = 5;
xlim([0 100]);
hold off
title('(base) alignment improvement')

%% demultiplexing
figure(3+figure_start)


subplot2(2,1,1)
h99=histogram(cell2mat(output(:,7+k)),'NumBins',14);
h99.Normalization = 'probability';
xlim([-1 12]);
set(gca,'xtick',[])
set(gca,'ytick',[])
ylim([.15 .4]);
subplot2(2,1,2)
h99.BinWidth=1;
h100=histogram(cell2mat(output(:,7+k)),'NumBins',14);
h100.Normalization = 'probability';
xlim([-1 12]);
ylim([0 .1]);
set(gca,'xtick',[])
set(gca,'ytick',[])
h100.BinWidth=1;
title('demultiplexing -- reads per sample')


%% counting IDs
A=f_countoutput(output);

%% SW alignment improvement
figure(5+figure_start);
h5 = histogram(cell2mat(output(:,14+k))/max(cell2mat(output(:,15+k))));
hold on
h6 = histogram(cell2mat(output(:,15+k))/max(cell2mat(output(:,15+k))));
h5.Normalization = 'probability';
h5.BinWidth = .02;
h5.FaceColor = c1;
h6.Normalization = 'probability';
h6.BinWidth = .02;
h6.FaceColor = c2;
xlim([0 1]);
hold off
title('swalignment improvement')




temp = output;
report = zeros(120,8);
for i=3:140
    X = cell2mat(output(:,3)); 
    A = cell2mat(output(:,15));
    B = cell2mat(output(:,16));
    if(~isempty(A(X==i)))
        report(i-2,:)=[i, median(A(X==i)), median(B(X==i)), mean(A(X==i)), mean(B(X==i)), min(A(X==i)), min(B(X==i)), sum(X==i)];
    end
end


end

