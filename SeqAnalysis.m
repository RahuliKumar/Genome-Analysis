function SeqAnalysis(Seq1,Seq2,Priorities,SaveFig,DataSave)
    if (strcmp(DataSave,'on'))
        SaveData(Seq1,Seq2);
    end
    DataPath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';
    % LET_LET Analysis
    AnalyDataName=[Seq1,'_',Seq2];AnalyFigName=[Seq1,'-',Seq2];
    % Local Sequence Alignment
%     AlignScoreMatrix=load([DataPath,'\AlignScoreMatrix_',AnalyDataName,'.mat']); % ASM: Alignment Score Matrix
%     AlignScoreMatrix=AlignScoreMatrix.AlignScoreMatrix;
    
    AlignIdentMatrix=load([DataPath,'\AlignIdentMatrix_',AnalyDataName,'.mat']); % ASM: Alignment Score Matrix
    AlignIdentMatrix=AlignIdentMatrix.AlignIdentMatrix;
    
    SortedInfoMatrix=load([DataPath,'\SortedInfoMatrix_',AnalyDataName,'.mat']);
    SortedInfoMatrix=SortedInfoMatrix.SortedInfoMatrix;
    
    GroupedChr=load([DataPath,'\ChrGrouping_',AnalyDataName,'.mat']);
    GroupedChr=GroupedChr.GroupedChr;
    
    % Alignment Score
    figure;imagesc(AlignIdentMatrix);
    colormap(hot);colorbar;title([AnalyFigName,' Sequence Comparison']);
    xlabel('Chromosome Number');ylabel('Chromosome Number');
    if (strcmp(SaveFig,'on'))
        SaveImage(['Alignment Score',AnalyFigName]);
    end
    % Pairwise Chromosome Grouping
    figure('Name',['Pairwise Grouping:',AnalyFigName]);
    Marker='o';LineStyle='--';Color='k';
    GroupingPlot_I(GroupedChr,Marker,LineStyle,Color,AnalyFigName,SaveFig);
    % Priority Analysis
    PriorPlots(Priorities,SortedInfoMatrix,AnalyFigName,SaveFig)
    % Chromosome Number Optimization
    ChrNumOptPlot_I(SortedInfoMatrix,AnalyFigName,SaveFig);
end