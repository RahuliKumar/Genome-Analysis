function PriorComparePlots(Priorities,SortedInfoMatrixI,SortedInfoMatrixII,AnalyFigName,SaveFig)
    S=size(Priorities);
    for P=1:S(2)
        NumOfPrior=Priorities(P);
        % 2D AND 3D PLOT
        [NumOfPriorMatrix,PriorMatrix]=PriorAnalyCompare(SortedInfoMatrixI,SortedInfoMatrixII,NumOfPrior);
        Prior3DPlot(NumOfPriorMatrix,PriorMatrix,NumOfPrior,AnalyFigName,SaveFig);
        Prior2DPlot(NumOfPriorMatrix,PriorMatrix,NumOfPrior,AnalyFigName,SaveFig);
    end
end