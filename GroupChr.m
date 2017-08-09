function GroupChr(SortedScoreMatrix,SavePath,AnalyName)
    
    %**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
    % GroupedChr: Stores Grouping result
    % UngroupedChr: Chromosomes to be grouped
    % GroupedChrNum: Chromosomes already grouped
    %**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
    
    %**==**==**==**==**==**==**[Algorithm Description]**==**==**==**==**==**==**
    % While loop searches for global maximum score from the all the local alignment scores and terminates if all the
    % chromosomes have been grouped;for loop searches for global maxima score from SortedScoreMatrix.If: Global maxima
    % is found and the grouping of chromosome is not clashing with any of already grouped chromosome, then grouping
    % occours and then these chromosomes are deleted from the ungrouped list. Else: If there is any clash with already
    % grouped chromosome, then the loop continues by deleting the clashing row from SortedScoreMatrix
    %**==**==**==**==**==**==**[Algorithm Description]**==**==**==**==**==**==**
    GroupedChr=zeros(8,2);UngroupedChr=1:16;ChrNum=16;GroupedChrNum=0;
    
    while (ChrNum~=0)
        TempMax=[-inf,0,0];
        for Chr=1:ChrNum
            TempChr=UngroupedChr(Chr);
            if (TempMax(1)<SortedScoreMatrix{TempChr,1}(1,2))
                TempMax=[SortedScoreMatrix{TempChr,1}(1,2),TempChr,SortedScoreMatrix{TempChr,1}(1,1)];
            end
        end
        
        if (any(UngroupedChr==TempMax(3)))
            UngroupedChr=UngroupedChr(UngroupedChr~=TempMax(2) & UngroupedChr~=TempMax(3)); % deleting the grouped chromosome
            GroupedChr(GroupedChrNum+1,1:2)=[TempMax(2),TempMax(3)];
            GroupedChrNum=GroupedChrNum+1;
            ChrNum=ChrNum-2;
         else
            SortedScoreMatrix{TempMax(2),1}(1,:)=[];
        end
    end
    save([SavePath,'\ChrGrouping_',AnalyName,'.mat'],'GroupedChr');
end