function SaveAnalyData(Seq1,Seq2)
    SeqPath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\';
    DataSavePath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';
    AnalyName=[Seq1,'_',Seq2];
    AlignScoreMatrix=LocalSeqAlign(Seq1,Seq2,SeqPath,DataSavePath);
    SortedScoreMatrix=ChrNumOpt(AlignScoreMatrix,DataSavePath,AnalyName);
    GroupChr(SortedScoreMatrix,DataSavePath,AnalyName);
end