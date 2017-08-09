function SeqAnalysisCompare(Analy1,Analy2,Analy3,SaveFig)
    
    DataPath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';
    
    AnalyFigNameI=[Analy1(1:3),'-',Analy1(5:7)];
    AnalyFigNameII=[Analy2(1:3),'-',Analy2(5:7)];
    AnalyFigNameIII=[Analy3(1:3),'-',Analy3(5:7)];
    LegendName={AnalyFigNameI,AnalyFigNameII,AnalyFigNameIII};
    
    % Chromosome Number Optimization
    SortedInfoMatrixI=load([DataPath,'\SortedInfoMatrix_',Analy1,'.mat']);
    SortedInfoMatrixI=SortedInfoMatrixI.SortedInfoMatrix;
    SortedInfoMatrixII=load([DataPath,'\SortedInfoMatrix_',Analy2,'.mat']);
    SortedInfoMatrixII=SortedInfoMatrixII.SortedInfoMatrix;
    SortedInfoMatrixIII=load([DataPath,'\SortedInfoMatrix_',Analy3,'.mat']);
    SortedInfoMatrixIII=SortedInfoMatrixIII.SortedInfoMatrix;
    ChrNumOptPlot_II(SortedInfoMatrixI,SortedInfoMatrixII,SortedInfoMatrixIII,LegendName,SaveFig)
    
end
