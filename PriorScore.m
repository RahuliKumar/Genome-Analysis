function ScoreList=PriorScore(NumOfPriorMatrix,NumOfPrior)
    ScoreList=zeros(16,1);
    for Chr=1:16
        TempScore=0;
        for Prior=1:NumOfPrior
            Weights=(1/NumOfPrior)*(NumOfPrior+1-Prior); %  Weights uniformly distributed among each priority.
            TempScore=TempScore+(Weights)*NumOfPriorMatrix(Chr,Prior);
        end
        ScoreList(Chr,1)=TempScore;
    end
end