figure('Name','Cumulative Sum Plot','NumberTitle','on');
p=plot(CumSum);p.Marker='s';

% for naming of chromosome on each point
for TempName=1:15
    txt =['  ',num2str(TempMatrix(TempName,1))];
    TempY=CumSum(TempName);
    if (TempName>8)
        TempY=TempY-(TempY)/35;
        text(TempName-.3,TempY,txt)
    else
        text(TempName,TempY,txt)
    end
    if (ChrNum==1)
        if (TempName==10)
            txt={'No drastic change','in slope observed'};
            text(3,1.2*TempY,txt);
        end
        if (TempName==5)
            txt={'\leftarrow This point corresponds to','     cumulative sum of LETLAS','     for Chromosome 1 with 3,4,','     14,13,15 & 12'};
            text(6.8,TempY-TempY/40,txt);
        end
    end
    % information for chromosome 2
    if (ChrNum==2 && TempName==9)
        txt='Drastic Change of Slope \rightarrow';
        text(2.9,TempY+(TempY)/30,txt);
    end
end
title(['Chromosome:',num2str(ChrNum)]);
xlabel('Number of Chromosomes');ylabel('LETLAS Cumulative Sum');
xlim([0,16]);ylim([0,CumSum(15)+CumSum(15)/10]);
