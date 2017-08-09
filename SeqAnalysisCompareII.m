function SeqAnalysisCompareII(Analy1,Analy2,Priorities,SaveFig)
    
    DataPath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';
    
    AnalyFigNameI=[Analy1(1:3),'-',Analy1(5:7)];
    AnalyFigNameII=[Analy2(1:3),'-',Analy2(5:7)];
    AnalyFigName=[AnalyFigNameI,' & ',AnalyFigNameII];
    
    % Pairwise Chromosome Grouping Compare
    
    GroupedChrI=load([DataPath,'\ChrGrouping_',Analy1,'.mat']);
    GroupedChrI=GroupedChrI.GroupedChr;
    GroupedChrII=load([DataPath,'\ChrGrouping_',Analy2,'.mat']);
    GroupedChrII=GroupedChrII.GroupedChr;
    
    figure('Name',['Pairwise Grouping Compare:',Analy1,'&',Analy2]);
    MarkerI='o';LineStyle='--';ColorI='k'; MarkerII='p';ColorII='b';
    GroupingPlot_II(GroupedChrI,GroupedChrII,MarkerI,MarkerII,LineStyle,ColorI,ColorII,AnalyFigNameI,AnalyFigNameII,SaveFig);
     
    % Priority Analysis Comaprision
    SortedInfoMatrixI=load([DataPath,'\SortedInfoMatrix_',Analy1,'.mat']);
    SortedInfoMatrixI=SortedInfoMatrixI.SortedInfoMatrix;
    SortedInfoMatrixII=load([DataPath,'\SortedInfoMatrix_',Analy2,'.mat']);
    SortedInfoMatrixII=SortedInfoMatrixII.SortedInfoMatrix;
    PriorComparePlots(Priorities,SortedInfoMatrixI,SortedInfoMatrixII,AnalyFigName,SaveFig);
   
end
