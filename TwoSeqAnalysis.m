function TwoSeqAnalysis(Seq1,Seq2,SaveFig,DataSave)
    if (strcmp(DataSave,'on'))
        SaveData(Seq1,Seq2);
    end
    DataPath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';
    % LET_LET Analysis
    AnalyName=[Seq1,'_',Seq2];AnalyFigName=[Seq1,'-',Seq2];
    % Local Sequence Alignment
    AlignScoreMatrix=load([DataPath,'\AlignScoreMatrix_',AnalyName,'.mat']); % ASM: Alignment Score Matrix
    AlignScoreMatrix=AlignScoreMatrix.AlignScoreMatrix;
    
    SortedInfoMatrix=load([DataPath,'\SortedInfoMatrix_',AnalyName,'.mat']);
    SortedInfoMatrix=SortedInfoMatrix.SortedInfoMatrix;
    
    GroupedChr=load([DataPath,'\ChrGrouping_',AnalyName,'.mat']);
    GroupedChr=GroupedChr.GroupedChr;
    
    % Alignment Score
    imagesc(AlignScoreMatrix);
    colormap(hot);colorbar;title('LET-LET Sequence Comparison');
    xlabel('Chromosome Number');ylabel('Chromosome Number');
    % Pairwise Chromosome Grouping
    figure('Name',['Pairwise Grouping:',AnalyFigName]);
    Marker='o';LineStyle='--';Color='k';
    GroupingPlot_I(GroupedChr,Marker,LineStyle,Color,AnalyFigName,SaveFig);
    % Priority Analysis
    Priorities=[5,9,15];PriorPlots(Priorities,SortedInfoMatrix,AnalyFigName,SaveFig)
    % Chromosome Number Optimization
    ChrNumOptPlot_I(SortedInfoMatrix,AnalyFigName,SaveFig);
end