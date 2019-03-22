function [output]=findFP(seqarrs,FP,RP,rows,columns,figurenum,normalize_y_n)

[a b]=size(seqarrs);
%a is the number of fastq files read in
%b is the number of sequences in the files!
passing=zeros(a,b);
figure(figurenum);
maxFP=swalign(FP,FP,'Alphabet','NT');
maxRP=swalign(RP,RP,'Alphabet','NT');
for q=1:a
    for i=1:b
        tmp=char(seqarrs(q,i))
        z=length(tmp);
        l1=length(FP);
        l2=length(RP);
        o=zeros(1,z-l1);
        r=zeros(1,z-l2);
        for p=1:(z-l1)
            if(normalize_y_n==0)
                o(1,p)=swalign(tmp(p:p+l1),FP,'Alphabet','NT');
                r(1,p)=swalign(tmp(p:p+l2),RP,'Alphabet','NT');
            else
                o(1,p)=swalign(tmp(p:p+l1),FP,'Alphabet','NT')/maxFP;
                r(1,p)=swalign(tmp(p:p+l2),RP,'Alphabet','NT')/maxRP;                
            end
        end
        x=[1;z-l1];
        x2=[1;z-l2];
        if(normalize_y_n==0)
            ymax=[17;17];
            ymaxmax=[38.82];
            ymean=[7.75;7.75];
            y3s=[7.75+3*1.8;7.75+3*1.8];
        else
            ymax=[0.5,0.5];
            ymaxmax=[1;1];
            ymean=[0.2;0.2];
            y3s=[0.443;0.443];            
        end
        subplot(rows,columns,i)
        plot(1:(z-l1),o,'b',x,ymax,x,ymean,x,y3s,x,ymaxmax,1:(z-l2),r,'r')
        title(strcat('sequence',num2str(i)))
    end
end

end