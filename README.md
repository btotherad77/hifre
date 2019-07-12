# hifre
Code to accompany the HiFRe manuscript (link)

Published in Analytical Chemistry: https://pubs.acs.org/doi/pdf/10.1021/acs.analchem.9b00856

*****************************************************************************************************
***************   note: this code is shared in interest of scientific transparency.   ***************

***************    it is not user friendly or readily adaptable to other scenarios    ***************

*****************************************************************************************************


The code is designed to take fastq files directly from the output of the MinKNOW software with no preprocessing on other software. Also, please note that modications will be required to adapt the code to different sequences.

1. Put the .fastq files into a folder in the MATLAB workspace. Call InitializeVars(); to set up constants.
2. Call [header, sequences] = a_readAll('your folder name'); 
3. dist = b_lengthdist(sequences); will give you the distribution of lengths in your unprocessed sequencing data.
4. [output,alignments,calignments,seqs,errorcount_before,errorcount_after]=c_locateRepeats(sequences,FP,1,1,1,0,1,0,0)
  the output variable is the bread and butter and contains the majority of the analytical information--the other variables are mostly for error checking.
  output is a 17-column cell matrix, the rows correspond to different raw reads and the columns contain the following information: 
      1. pass count
      2. original sequence # (i.e. the location of this raw read in the original sequences array)
      3. number of repeats
      4. average alignment of repeats
      5. alignment of consensus sequence
      6. barcode read
      7. assigned demultiplexed barcode
      8. index of assigned barcode
      9. quality score (using a pseudo-phred score that I made up for internal analysis)
      10. quality score (again, kind of made up, but it's kind of like an average accuracy of consensus)
      11. ind -- ID of the sequence -- either oligo 1 or oligo 2
      12. forward/reverse read?
      13. average of alignments over just the MIP region (NeedlemanWunsch) 
      14. alignment of consensus over just the MIP region (NeedlemanWunsch)
      15. average alignment of repeats (Smith Waterman)
      16. alignment of consensus sequence (Smith Waterman)
      17. identity of SNV locations]  
5. [ report] = examineOutput( output , 1)
      This will plot a number of useful figures


  
  
