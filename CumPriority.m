function CumPriority(NumOfPriorMatrix,NumOfPrior,AnalyFigName,Save)
    ScoreList=zeros(16,15);TempMatrix=zeros(16,2);hold all;
    for Chr=1:16
        TempScore=0;
        for Prior=1:NumOfPrior
            Weights=(1/NumOfPrior)*(NumOfPrior+1-Prior);%  Weights uniformly distributed among each priority.
            TempScore=TempScore+(Weights)*NumOfPriorMatrix(Chr,Prior);
            ScoreList(Chr,Prior)=TempScore;
        end
        plot(ScoreList(Chr,1:Prior),'LineWidth',1.5);
        TempMatrix(Chr,:)=[Chr,ScoreList(Chr,Prior)];
    end
    TempMatrix=sortrows(TempMatrix,2);
    for Chr=1:16
        if (rem(Chr,2))
            text(NumOfPrior+0.05,TempMatrix(Chr,2),num2str(TempMatrix(Chr,1)));
        else
            text(NumOfPrior+0.1,TempMatrix(Chr,2),num2str(['\rightarrow ',num2str(TempMatrix(Chr,1))]));
        end
    end
    hold off;
    title([AnalyFigName,':',num2str(NumOfPrior),' Priorities']);xlim([0,NumOfPrior+1]);xlabel('Priority Number');ylabel('Cumulative Score [Weighted Sum]');
    if (strcmp(Save,'on'))
        SaveImage(['CumPriority_',AnalyFigName,num2str(NumOfPrior),'Priorities']);
    end
    
end