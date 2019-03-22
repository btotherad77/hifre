function [output,alignments,calignments,seqs]=c_locateRepeats(seqarrs,FP,rows,columns,figurenum,ploty_n,check_rev,debug,write_to_file)

output = [];
[ a b ] = size(seqarrs);
%a is the number of fastq files read in
%b is the number of sequences in the files!

%FP='GGCAGCGTCAGATGTGTATAAGAGACAG';

tic
start_time = tic;

% set up counter
file_counter = 0;
passing = zeros(a,b);
d_matrix = 'NUC44';
maxFP=swalign(FP,FP,'Alphabet','NT');
thresh = 0.423;%% determined threshold for alignment 

up_flag_fp = 1;
fp_count = 0;
alignments = [];
calignments = [];
alignments_sw = [];
calignments_sw = [];
base_alignments = [];
base_calignments = [];
seqs = [];
pass_count = 0;

%repeat length:
Length = 156;
CIP =           'GGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGGTTATGCAACTAATCTCCGAGCTAATATGATCACTAATGTCTGCANNNNNNNNNNNNACTTACTCTAGGTATG';

CIP_OLG1 =       'GGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNTGTAGAGATGACTATCGTACCCTGCACGGTACNNNNNNNNNNAAGGTTATGCAACTAATCTCCGAGCTAATATGATCACTAATGTCTGCANNNNNNNNNNNNACTTACTCTAGGTATG';
CIP_OLG2 =       'GGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNGTCTACATCCCTGAGTCTATATGATGTAGATANNNNNNNNNNAAGGTTATGCAACTAATCTCCGAGCTAATATGATCACTAATGTCTGCANNNNNNNNNNNNACTTACTCTAGGTATG';
CIP_OLG2_SNP =   'GGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNGTCTACATCGCTGAGTCTTTATGATGTACATANNNNNNNNNNAAGGTTATGCAACTAATCTCCGAGCTAATATGATCACTAATGTCTGCANNNNNNNNNNNNACTTACTCTAGGTATG';


max_base_alignment = nwalign(CIP,CIP,'Alphabet','NT')
max_ID_alignment_gfp = nwalign(CIP_OLG2,CIP_OLG2,'Alphabet','NT')
miss_ID_alignment = nwalign(CIP_OLG1,CIP_OLG2,'Alphabet','NT')
snp_ID_alignment = nwalign(CIP_OLG2_SNP,CIP_OLG2,'Alphabet','NT')

subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.02], [0.03 0.03], [0.02 0.01]);
fig_count=1;
plot_count=1;
errorcount_before=containers.Map;
errorcount_after=containers.Map;
SNPfinder=containers.Map;

fp_len=length(FP);
for q=1:a %iterate through files
    for i=1:b %iterate through sequences
        failed=0;
        FPpoints_x = [];
        FPpoints_y = [];  %y is the alignment value
        FPpoints_i = [];  %i is the additional length to take subsequence from
        
        rev_flag=1;
        tmp=char(seqarrs(q,i)); %get full read #i from file a
        z=length(tmp); % z = length of read
        
        if(mod(i,500)==1)%print status update
            [i, pass_count, toc ]
        end 
        
        if(z > 156*3 ) %length cutoff --> skip read if too short...
            for REV = 1:2 %go through again and test revcomp
                if(rev_flag == 1)
                    FPpoints_x = [];
                    FPpoints_y = [];
                    FPpoints_i = []; 
                    
                    increment = 0;
                    
                    f = zeros(1,z-fp_len); %alignment score based on position in read
                    for p=1:(z-fp_len) % iterate basewise through sequence
                            f(1,p)=swalign(tmp(p:p+fp_len-1),FP,'Alphabet','NT')/maxFP; 
                            if(f(1,p)>=thresh) %alignment at position p passes defined threshold
                                if(up_flag_fp==1) %starts at 1
                                    %if current alignment is greater than the last
                                    %alignment, then replace it with this one
                                    fp_count=fp_count+1;
                                    if(fp_count==1)
                                       FPpoints_x=[FPpoints_x;p];       %x is the position
                                       FPpoints_y=[FPpoints_y;f(1,p)];  %y is the alignment value
                                       FPpoints_i=[FPpoints_i;0];  %y is the alignment value
                                    else
                                        num = length(FPpoints_y);
                                        if(num>0) %prevents it from messing up if sequence starts in middle of FP
                                           if(f(1,p)>=FPpoints_y(num))
                                            if(f(1,p)==FPpoints_y(num))
                                                %FPpoints_i(num) = FPpoints_i(num) + 1;  %y is the alignment value
                                            else
                                                %FPpoints_i(num) = 0;  %y is the alignment value
                                                FPpoints_x(length(FPpoints_x))=p;       %x is the position
                                                %^^ only increment position
                                                %marker if the alignment is
                                                %strictly greater
                                            end
                                            %FPpoints_x(length(FPpoints_x))=p;       %x is the position
                                            FPpoints_y(length(FPpoints_y))=f(1,p);  %y is the alignment value
                                           else
                                           %if alignment is decreasing overall, set flag to
                                           %false and set fp counter to zero
                                               up_flag_fp=0;
                                               fp_count=0;
                                           end
                                        end %end if num>0
                                    end%end if fp_count == 1
                                end %end if up_flag == 1
                            else
                                %reset flag once dropping below threshold
                                up_flag_fp=1;
                                fp_count=0;
                            end               
                    end
                    if( length(FPpoints_x) < 3 && check_rev == 1 )
                        tmp=seqrcomplement(tmp);
                    else
                        rev_flag = 0; %prevent from rechecking for rev comp.
                    end %end of check for rev comp
                end
            end
        %%%%%%%%%%%%%% PLOTTING FP ALIGNMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if( plot_count <= 2100 && ploty_n == 1 )
                        if( failed == 0 )
                            x=[1;z-fp_len];

                                ymax=[0.5,0.5];
                                ymaxmax=[1;1];
                                ymean=[0.2;0.2];
                                y3s=[0.443;0.443];
                            figure(figurenum+ceil(plot_count/(rows*columns)))
                            subplot_index=mod(plot_count,rows*columns);
                            if(subplot_index==0)
                                subplot_index=rows*columns;
                                if(1==1)
                                    set(gcf,'Position',[000,50,1900,940])
                                    saveas(gcf,strcat('Images/2018_04_25/FindRepeats_OmitLessThan6Hits_1-200_FP+FPc_',int2str(fig_count),'.png'));
                                end
                                fig_count = fig_count+1;
                            end
                            subplot(rows,columns,subplot_index)
                            if(rev_flag==0)
                                plot(1:(z-fp_len),f,'b',x,ymax,x,ymean,x,y3s,x,ymaxmax,FPpoints_x,FPpoints_y,'k*')
                            else
                                plot(1:(z-fp_len),f,'r',x,ymax,x,ymean,x,y3s,x,ymaxmax,FPpoints_x,FPpoints_y,'k*')
                            end
                            title(strcat('sequence',num2str(i)))
                            plot_count = plot_count+1;
                        end
                    end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%                          sequence alignment
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 SEQS = {};
                 temp_alignments = [];
                 temp_base_alignments = [];
                 temp_alignments_sw = [];
                 strcell={};
                 cell_count=1;
                 if(length(FPpoints_x) > 3) %if more than 3 repeats found in sequence
                     passing(q,i) = 1;
                     for j=1:(length(FPpoints_x))
                         [~,nn]=size(tmp);
                         if( FPpoints_x(j) + Length + 0*FPpoints_i(j) - 1 <= nn )
                            temp_seq=tmp(FPpoints_x(j):(FPpoints_x(j) + Length + 0*FPpoints_i(j) -1));
                            
                            %save sequences for seqviewer
                            strcell{cell_count}=temp_seq;
                            SEQS{cell_count}=temp_seq;
                            cell_count=cell_count+1;
                                                
                            %%assign sequence to oligo or GFP        
                            [ind,realign_test_seq]=assignID(temp_seq,CIP_OLG1,CIP_OLG2);
                            
                            if(debug)
                                [j,length(temp_seq),length(realign_test_seq)]
                            end
                            if(debug)
                                temp_seq
                            end
                            [sc,aln]=nwalign(temp_seq,realign_test_seq,'Alphabet','NT');
                            [sc_sw,aln_sw]=swalign(temp_seq,realign_test_seq,'Alphabet','NT');
                            if(debug)
                               % sc
                               % aln
                            end
                            alignments=[alignments;sc];
                            temp_alignments=[temp_alignments;sc];
                            temp_alignments_sw=[temp_alignments_sw;sc_sw];
                            %alignments_sw=[alignments_sw;sc_sw];
                            [testing1,testing2]=nwalign(temp_seq,CIP,'Alphabet','NT');%COUNT BASE ALIGNMENTS BEFORE ID ASSIGNMENT
                            if(debug)
                              %  testing1
                               % testing2
                            end
                            base_alignments=[base_alignments;testing1];
                            temp_base_alignments=[temp_base_alignments;testing1];
                         end
                     end

                     oo = cell_count-1;
                     if(oo>=3)%IF THERE ARE GREATER THAN 3 FOUND SEQUENCES
                        alignment=multialign(SEQS,'Weights','THG','ScoringMatrix',d_matrix);    %look into this

                        [x,qq,qs]=compilestringalign(alignment,0,length(calignments));         %look into this
                        % ^^ CONSENSUS,quality,Qscore
                         
                        if(length(x)>Length) %limit repeat to having length = Length (156)
                            x=x(1:Length);
                        end
                        
                        
                        %%assign sequence to oligo or GFP
                        [ind,realign_test_seq]=assignID(x,CIP_OLG1,CIP_OLG2);
                        
                        strcell{cell_count}=x;%consensus
                        if(ind == 1)%original
                            strcell{cell_count+1}=CIP_OLG1;
                        else
                            strcell{cell_count+1}=CIP_OLG2;
                        end
                        
                        if(debug)
                            %x
                             v_viewSeq(strcell);
                        end
                        
                        [c,aln]=nwalign(x,realign_test_seq,'Alphabet','NT');
                        [c_sw,aln_sw]=swalign(x,realign_test_seq,'Alphabet','NT');
                        if(debug)
                            c
                            aln
                        end
                        
                        
                        bc = nwalign(x,CIP,'Alphabet','NT');
                        for ggg = 1:oo %copy consensus alignment times the number of repeats so that it has the same # of entries as individual alignments
                            calignments = [calignments;c];
                            calignments_sw = [calignments_sw;c_sw];
                            base_calignments = [base_calignments;bc];
                        end
                        %% assign SNP (attempt #1)
                        SNPs = 'N/A';
                        if(ind == 2) %if the alignment was to GFP and not to OLIGO
                            SNPs = findSNPs(x,CIP_OLG2);
                        end
                        if(SNPfinder.isKey(strcat(SNPs)))
                            SNPfinder(strcat(SNPs))=SNPfinder(strcat(SNPs))+1;
                        else
                            SNPfinder(strcat(SNPs))=1;
                        end
                        
                        
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
                        %11: ID (1-oligo, 2-GFP)
                        %12: rev_flag
                        %13: average base individual alignment
                        %14: average base consensus alignment
                        %15: avg swalignment
                        %16: sw_score post alignment
                        %17: SNPs
                        pass_count = pass_count + 1;
                        output = [output;pass_count,i,length(FPpoints_x),mean(temp_alignments),c,x_assignbarcode2(x,debug),qq,qs,ind,rev_flag,mean(temp_base_alignments),bc,mean(temp_alignments_sw),c_sw,cellstr(SNPs)];
                        cell_count
                       % AA = zeros(cell_count,cell_count);
                        if(debug && pass_count == 5)
                            %generateTree(strcell)
                            v_viewSeq(strcell);
                        disp('Press a key !')  % Press a key here.You can see the message 'Paused: Press any key' in        % the lower left corner of MATLAB window.
                        pause;
                        end
                        %original__consensus
                        seqs = [seqs;cellstr(tmp),cellstr(x)];
                        
                        if(write_to_file == 1)
                            s = 1000;
                            mod(pass_count,s);
                            if(mod(pass_count,s) == 0)
                                [pass_count s]
                                strcat('2018_05_15-Analysis(',num2str(file_counter),')',num2str((file_counter*s+1)),'-',num2str(pass_count))
                                xlswrite(strcat('2018_05_15-Analysis(',num2str(file_counter),')',num2str((file_counter*s+1)),'-',num2str(pass_count)),output((file_counter*s+1):pass_count,:));
                                file_counter = file_counter + 1;
                            end
                            
                        end
                     end
                 end
        end
    end
end

end
function [ID,sequence] = assignID(x,s1,s2)
    [~,ID]=max([swalign(x,s1,'Alphabet','NT'),swalign(x,s2,'Alphabet','NT')]); %assign ID to consensus sequence
    if(ID==1)
           sequence=s1;
    else
           sequence=s2;
    end
end
function [output,quality,Qscore] = compilestringalign(STRS,ploty_n,q)
%[x y]=size(STRS(:,1:10));
[x, y] = size(char(STRS));
%counter=zeros(1,y);
for j=1:y
    counter('A',j)=0;
    counter('T',j)=0;
    counter('G',j)=0;
    counter('C',j)=0;
    counter('-',j)=0;
end
for i=1:x
    for j=1:y
        if(STRS(i,j)=='A'||STRS(i,j)=='T'||STRS(i,j)=='G'||STRS(i,j)=='C'||STRS(i,j)=='-')
            counter(STRS(i,j),j)=counter(STRS(i,j),j)+1/x;
        end
    end
end

output='';
letters=['A' 'T' 'C' 'G' '-'];
MAX=zeros(y,2);
I=1:y;
quality='';
for j=1:y
    [m, i]=max([counter('A',j) counter('T',j) counter('C',j) counter('G',j) counter('-',j)]);
    MAX(j,:)=[j,m];
    I(j)=i;
    if(i~=5)
        quality=strcat(quality,char(65-20*log10(m)));
        output=strcat(output,letters(i));
    end
end
Qscore=mean(MAX(:,2));
if(ploty_n~=0)
    N=(ploty_n-2)/1000;
    subplot(1,1,1)
    figure(ploty_n)
    plot(1:y,counter('A',:),'r',1:y,counter('T',:),'c',1:y,counter('G',:),'g',1:y,counter('C',:),'y',1:y,MAX,'m--');
    set(gca,'Color',[0 0 0]);
    ylim([0 1]);
    figure(ploty_n+1)
    subplot(2,1,2)
    subplot(1,2,2)
    b=bar(MAX');
    colors=['r' 'c' 'y' 'g' 'w'];
    for j=1:y
        set(b(j),'FaceColor',colors(I(j)));
        set(b(j),'EdgeColor','k');
    end
    xlim([1.6 2.4])
    set(gca,'Color',[0 0 0]);
    ylim([0 1]);
end
end