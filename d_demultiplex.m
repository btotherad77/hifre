function [ BCs, b1arr,b2arr,b3arr,b4arr,bXarr,NAarr ] = d_demultiplex( s )
%D_DEMULIPLEX Summary of this function goes here
%   Detailed explanation goes here
b1 = 'GTCATCTC';
b2 = 'CCAGGTCC';
b3 = 'TCCATCCT';
b4 = 'AGCTCGCA';
bX = '--------';
b1c = seqrcomplement(b1);
b2c = seqrcomplement(b2);
b3c = seqrcomplement(b3);
b4c = seqrcomplement(b4);
CIP='GGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNGTCTACATCCCTGAGTCTATATGATGTAGATANNNNNNNNNNAAGGTTATGCAACTAATCTCCGAGCTAATATGATCACTAATGTCTGCATANNNNNNNNTAACTTACTCTAGGTATGGGCAGCGTCAGATGTGTATAAGAGACAG';
CIPr='CTGTCTCTTATACACATCTGACGCTGCCCATACCTAGAGTAAGTTANNNNNNNNTATGCAGACATTAGTGATCATATTAGCTCGGAGATTAGTTGCATAACCTTNNNNNNNNNNTATCTACATCATATAGACTCAGGGATGTAGACNNNNNNNNNNCTGTCTCTTATACACATCTGACGCTGCC';
failed=zeros(length(s),1);
% MBCs('b1arr',1:2)=[cellstr('start'),cellstr('start')]';
% MBCs('b2',1:2)=[cellstr('start'),cellstr('start')];
% MBCs('b3',1:2)=[cellstr('start'),cellstr('start')];
% MBCs('b4',1:2)=[cellstr('start'),cellstr('start')];
% MBCs('bX',1:2)=[cellstr('start'),cellstr('start')];
% MBCs('N/A',1:2)=[cellstr('start'),cellstr('start')];
b1arr=[];
b2arr=[];
b3arr=[];
b4arr=[];
bXarr=[];
NAarr=[];
BCs=[];
% BCtest=[];
%   for i=1:4
%         if(i==1)
%             barcode=b1;
%         elseif(i==2)
%             barcode=b2;
%         elseif(i==3)
%             barcode=b3;
%         elseif(i==4)
%             barcode=b4;
%         end
%         BCtest=[BCtest;cellstr(barcode),swalign(barcode,b1c, 'Alphabet' , 'NT'),swalign(barcode,b2c, 'Alphabet' , 'NT'),swalign(barcode,b3c, 'Alphabet' , 'NT'),swalign(barcode,b4c, 'Alphabet' , 'NT'),swalign(barcode,bX, 'Alphabet' , 'NT')];
%        BCtest=[BCtest;cellstr(barcode),swalign(barcode,b1, 'Alphabet' , 'NT'),swalign(barcode,b2, 'Alphabet' , 'NT'),swalign(barcode,b3, 'Alphabet' , 'NT'),swalign(barcode,b4, 'Alphabet' , 'NT'),swalign(barcode,bX, 'Alphabet' , 'NT')];
%   end
%     BCtest
for i=1:length(s)
    [score, alignment] = swalign( char(s(i)) , CIP , 'Alphabet' , 'NT' );
    k = strfind(alignment(3,:),'TANNNNNNNNTA');
    k1 = strfind(alignment(3,:),'AGNNNNNNNNNNGT');
    k2 = strfind(alignment(3,:),'TANNNNNNNNNNAA');
    revflag = 0;
    if( score < 100 )
        revflag = 1;
        [score, alignment] = swalign( char(s(i)) , CIPr , 'Alphabet' , 'NT' );
        k = strfind(alignment(3,:),'TANNNNNNNNTA');
        k1 = strfind(alignment(3,:),'AGNNNNNNNNNNGT');
        k2 = strfind(alignment(3,:),'TANNNNNNNNNNAA');
        
        if( score < 100 )
        
            failed(i) = 1;
        end
    else
    end
    %alignment
    %k
    barcode=alignment(1,k+2:k+9);
    MBC1=alignment(1,k1+2:k1+11);
    MBC2=alignment(1,k2+2:k2+11);
  
       
        
    cutoff=6.775;
    if(isempty(barcode))
        BCs=[BCs;cellstr(barcode),-1,-1,-1,-1,-1,-1,cellstr('N/A'),-1];
    else
        if(revflag == 0)
            if(swalign(barcode,b1c, 'Alphabet' , 'NT') >= cutoff)
                result = cellstr('b1');
                result2 = 1;
                b1arr=[b1arr;cellstr(MBC1),cellstr(MBC2)];
            else
                if(swalign(barcode,b2c, 'Alphabet' , 'NT') >= cutoff)
                    result = cellstr('b2');
                    result2 = 2;
                    b2arr=[b2arr;cellstr(MBC1),cellstr(MBC2)];
                else
                    if(swalign(barcode,b3c, 'Alphabet' , 'NT') >= cutoff)
                        result = cellstr('b3');
                        result2 = 3;
                        b3arr=[b3arr;cellstr(MBC1),cellstr(MBC2)];
                    else
                        if(swalign(barcode,b4c, 'Alphabet' , 'NT') >= cutoff)
                            result = cellstr('b4');
                            result2 = 4;
                            b4arr=[b4arr;cellstr(MBC1),cellstr(MBC2)];
                        else
                            if(barcode == '--------')
                                result = cellstr('bX');
                                result2 = 0;
                                bXarr=[bXarr;cellstr(MBC1),cellstr(MBC2)];
                            else
                                result = cellstr('N/A');
                                NAarr=[NAarr;cellstr(MBC1),cellstr(MBC2)];
                                result2 = -1;
                            end
                        end
                    end
                end
            end
            BCs=[BCs;cellstr(barcode),revflag,swalign(barcode,b1c, 'Alphabet' , 'NT'),swalign(barcode,b2c, 'Alphabet' , 'NT'),swalign(barcode,b3c, 'Alphabet' , 'NT'),swalign(barcode,b4c, 'Alphabet' , 'NT'),swalign(barcode,bX, 'Alphabet' , 'NT'),result,result2];
            
        else
          if(swalign(barcode,b1, 'Alphabet' , 'NT') >= cutoff)
                result = cellstr('b1');
                result2 = 1;
                b1arr=[b1arr;cellstr(MBC1),cellstr(MBC2)];
            else
                if(swalign(barcode,b2, 'Alphabet' , 'NT') >= cutoff)
                    result = cellstr('b2');
                    result2 = 2;
                    b2arr=[b2arr;cellstr(MBC1),cellstr(MBC2)];
                else
                    if(swalign(barcode,b3, 'Alphabet' , 'NT') >= cutoff)
                        result = cellstr('b3');
                        result2 = 3;
                        b3arr=[b3arr;cellstr(MBC1),cellstr(MBC2)];
                    else
                        if(swalign(barcode,b4, 'Alphabet' , 'NT') >= cutoff)
                            result = cellstr('b4');
                            result2 = 4;
                            b4arr=[b4arr;cellstr(MBC1),cellstr(MBC2)];
                        else
                            if(barcode == '--------')
                                result = cellstr('bX');
                                result2 = 0;
                                bXarr=[bXarr;cellstr(MBC1),cellstr(MBC2)];
                            else
                                result = cellstr('N/A');
                                NAarr=[NAarr;cellstr(MBC1),cellstr(MBC2)];
                                result2 = -1;
                            end
                        end
                    end
                end
            end
            BCs=[BCs;cellstr(barcode),revflag,swalign(barcode,b1, 'Alphabet' , 'NT'),swalign(barcode,b2, 'Alphabet' , 'NT'),swalign(barcode,b3, 'Alphabet' , 'NT'),swalign(barcode,b4, 'Alphabet' , 'NT'),swalign(barcode,bX, 'Alphabet' , 'NT'),result,result2];
        end
    end
end
% x=b1arr(:,1)
% y=b1arr(:,2)
% b1arr=b1arr(~(cellfun('isempty',x)&&cellfun('isempty',y)));
% x=b2arr(:,1)
% y=b2arr(:,2)
% b2arr=b2arr(~cellfun('isempty',x)&&~cellfun('isempty',y));
% x=b3arr(:,1)
% y=b3arr(:,2)
% b3arr=b3arr(~cellfun('isempty',x)&&~cellfun('isempty',y));
% x=b4arr(:,1)
% y=b4arr(:,2)
% b4arr=b4arr(~cellfun('isempty',x)&&~cellfun('isempty',y));
% x=bXarr(:,1)
% y=bXarr(:,2)
% bXarr=bXarr(~cellfun('isempty',x)&&~cellfun('isempty',y));

end

