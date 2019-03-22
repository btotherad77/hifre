function [ SNPs ] = findSNPs( x,matching_strand )
    [~,snp_alignment,~] = nwalign(x,matching_strand,'Alphabet','NT');
    %i1 = 48 C-->G
    %i2 = 57 A-->T
    %i3 = 67 G-->C
    
    %expected results:
    %GFP ==> 'CAG'
    %SNP ==> 'GTC'
    
    i_star = 0;
    raw = snp_alignment(1,:);
    template = snp_alignment(3,:);
    SNPs='';
    for index = 1:length(snp_alignment)
        if( template(index) ~= '-')
            i_star = i_star + 1;
        end
        if(i_star == 48)
            SNPs = strcat(SNPs,raw(index));
        end
        if(i_star == 57)
            SNPs = strcat(SNPs,raw(index));
        end
        if(i_star == 67)
            SNPs = strcat(SNPs,raw(index));
        end
    end
end

