function [SortedScoreMatrix,SortedInfoMatrix]=ChrNumOpt(ScoreMatrix,SavePath,AnalyName)
    % For Chromosome Number Optimization;
    
    %**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
    % SortedInfoMatrix:Stores sorted information like alignment score, cummulative sum and corresponding chr num.
    % SortedScoreMatrix: Stores sorted alignment score with chr num
    %**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
    SortedInfoMatrix=cell(1);SortedScoreMatrix=cell(16,1);
    for ChrNum=1:16
        ColNum=1+(ChrNum-1)*4;
        TempMatrix=zeros(16,2);TempMatrix(:,1)=1:16;
        TempMatrix(:,2)=transpose(ScoreMatrix(ChrNum,:));
        TempMatrix(ChrNum,:)=[];TempMatrix=sortrows(TempMatrix,-2);
        CumSum=cumsum(TempMatrix(:,2));SortedScoreMatrix{ChrNum,1}=TempMatrix;
        SortedInfoMatrix{1}(:,ColNum:ColNum+3)=[repmat(ChrNum,15,1),TempMatrix,CumSum];
    end
    save([SavePath,'\SortedInfoMatrix_',AnalyName,'.mat'],'SortedInfoMatrix');
    save([SavePath,'\SortedScoreMatrix_',AnalyName,'.mat'],'SortedScoreMatrix');
end
