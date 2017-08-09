function PriorPlots(Priorities,SortedInfoMatrix,AnalyFigName,SaveFig)
    S=size(Priorities);
    for P=1:S(2)
        NumOfPrior=Priorities(P);
        % 2D AND 3D PLOT
        [NumOfPriorMatrix,PriorMatrix]=PriorAnaly(SortedInfoMatrix,NumOfPrior);
        Prior3DPlot(NumOfPriorMatrix,PriorMatrix,NumOfPrior,AnalyFigName,SaveFig);
        Prior2DPlot(NumOfPriorMatrix,PriorMatrix,NumOfPrior,AnalyFigName,SaveFig);
        % Priority Score
        ScoreList=PriorScore(NumOfPriorMatrix,NumOfPrior);
        scatter(1:16,ScoreList,'o','filled','LineWidth',2,'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',[0 0 0]);xlabel('Chromosome Number');ylabel('Priority Score [Weighted Sum]');
        xlim([0,17]);hold on;plot([0,17],[mean(ScoreList),mean(ScoreList)],'--','LineWidth',1.5);
        title([AnalyFigName,' Score:',num2str(NumOfPrior),' Priorities']);
        if (strcmp(SaveFig,'on'))
            SaveImage(['Score',AnalyFigName,num2str(NumOfPrior),' Priorities']);
        end
        % Cumulative Priority
        figure;CumPriority(NumOfPriorMatrix,NumOfPrior,AnalyFigName,SaveFig);
    end
end