function SaveData(Seq1,Seq2)
    SeqPath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\';
    SavePath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';DataPath=SavePath;
    AnalyName=[Seq1,'_',Seq2];
%     disp('Executing LocalSeqAlign');AlignScoreMatrix=LocalSeqAlign(Seq1,Seq2,SeqPath,DataSavePath);
    AlignIdentMatrix=load([DataPath,'\AlignIdentMatrix_',AnalyName,'.mat']); % ASM: Alignment Score Matrix
    AlignIdentMatrix=AlignIdentMatrix.AlignIdentMatrix;
    disp('Executing ChrNumOpt');SortedScoreMatrix=ChrNumOpt(AlignIdentMatrix,SavePath,AnalyName);
    disp('Executing GroupChr');GroupChr(SortedScoreMatrix,SavePath,AnalyName);
end