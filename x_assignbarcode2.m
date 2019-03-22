function [ output ] = x_assignbarcode2( s , debug)
%D_DEMULIPLEX Summary of this function goes here
%   Detailed explanation goes here
b1 = 'TAGTCATCTCTA';
b2 = 'TACCAGGTCCTA';
b3 = 'TATCCATCCTTA';
b4 = 'TAAGCTCGCATA';
%%%%%%%%%%% NEED TO BE ASSIGNED %%%%%%%%%%%%%%%%%%
b5 = 'ATCCTCTCCTCA';
b6 = 'CGTCTACGATGC';
b7 = 'CACGAAGTGGAA';
b8 = 'GTTCCCTGTCCC';
b9 = 'ACGTTAAGGCCA';
b10 = 'TCCTCGTGAGGT';
b11 = 'CATCTAACCTAG';
b12 = 'TCGGAATTGGCT';%QBGJ866BK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bX = '------------';
b1c = seqrcomplement(b1);
b2c = seqrcomplement(b2);
b3c = seqrcomplement(b3);
b4c = seqrcomplement(b4);
b5c = seqrcomplement(b5);
b6c = seqrcomplement(b6);
b7c = seqrcomplement(b7);
b8c = seqrcomplement(b8);
b9c = seqrcomplement(b9);
b10c = seqrcomplement(b10);
b11c = seqrcomplement(b11);
b12c = seqrcomplement(b12);
%CIP='GGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNGTCTACATCCCTGAGTCTATATGATGTAGATANNNNNNNNNNAAGGTTATGCAACTAATCTCCGAGCTAATATGATCACTAATGTCTGCANNNNNNNNNNNNACTTACTCTAGGTATG';
CIP='GGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGGTTATGCAACTAATCTCCGAGCTAATATGATCACTAATGTCTGCANNNNNNNNNNNNACTTACTCTAGGTATG';
%CIPr='CATACCTAGAGTAAGTNNNNNNNNNNNNTGCAGACATTAGTGATCATATTAGCTCGGAGATTAGTTGCATAACCTTNNNNNNNNNNTATCTACATCATATAGACTCAGGGATGTAGACNNNNNNNNNNCTGTCTCTTATACACATCTGACGCTGCC';
CIPr='CTGTCTCTTATACACATCTGACGCTGCCCATACCTAGAGTAAGTNNNNNNNNNNNNTGCAGACATTAGTGATCATATTAGCTCGGAGATTAGTTGCATAACCTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN';
failed=0;


BC='*';
l=12;
    [score, alignment] = nwalign( char(s) , CIP , 'Alphabet' , 'NT' );
    if(debug)
        alignment
    end
%     if( score < 40 )%%REVERSE%%
%         revflag = 1;
%         [score, alignment] = nwalign( char(s) , CIPr , 'Alphabet' , 'NT' );
%         if(debug)
%            revflag
%            alignment
%         end
%         k = strfind(alignment(3,:),'TNNNNNNNNNNNNT');
%         if( score < 100 )
%             failed = 1;
%         end
%     else%%FORWARD%%
        k = strfind(alignment(3,:),'ANNNNNNNNNNNNA');
        if( isempty(k) )
            k = strfind(alignment(3,:),'A-NNNNNNNNNNNNA');
            l=13;
            bX = '-------------';
            if( isempty(k) )
                k = strfind(alignment(3,:),'ANNNNNNNNNNNN-A');
                l=13;
                bX = '-------------';
                if( isempty(k) )
                    k = strfind(alignment(3,:),'A-NNNNNNNNNNNN-A');
                    l=14;
                       bX = '--------------';
                end
            end
        end            
        revflag = 0;
   % end
    
    if(debug)
       k
    end
    %alignment
    barcode=alignment(1,k+1:k+l);
  
       
    temp=0;
        
    cutoff=6.775;
    
    
    ind = -1;                        
    if(isempty(barcode))
        BC='N/A';
    else
        [a,ind]=max([swalign(barcode,b1c,'Alphabet','NT'),swalign(barcode,b2c,'Alphabet','NT'),swalign(barcode,b3c,'Alphabet','NT'),swalign(barcode,b4c,'Alphabet','NT'),swalign(barcode,b5c,'Alphabet','NT'),swalign(barcode,b6c,'Alphabet','NT'),swalign(barcode,b7c,'Alphabet','NT'),swalign(barcode,b8c,'Alphabet','NT'),swalign(barcode,b9c,'Alphabet','NT'),swalign(barcode,b10c,'Alphabet','NT'),swalign(barcode,b11c,'Alphabet','NT'),swalign(barcode,b12c,'Alphabet','NT'),swalign(barcode,bX,'Alphabet','NT')]);
        if(a>cutoff)
            switch ind
                case 1  
                    BC=b1c;
                case 2  
                    BC=b2c;
                case 3  
                    BC=b3c;
                case 4  
                    BC=b4c;
                case 5  
                    BC=b5c;
                case 6  
                    BC=b6c;
                case 7  
                    BC=b7c;
                case 8  
                    BC=b8c;
                case 9  
                    BC=b9c;
                case 10  
                    BC=b10c;
                case 11  
                    BC=b11c;
                case 12  
                    BC=b12c;
                case 13  %% empty alignment
                    BC=bX;
                otherwise BC='N/A'; ind = 99; %% weird, should never happen 
            end
        else %% unassignable barcode
            BC='N/A';
            ind = 0;        
        end
    end
       
    BC=cellstr(BC);
    output=[BC,cellstr(barcode),ind];
end

