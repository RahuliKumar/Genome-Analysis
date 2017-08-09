function [NumOfPriorMatrix,PriorMatrix]=PriorAnaly(SortedInfoMatrix,NumOfPrior)
    %**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
    % NumOfPrior: Number of priorities to be analyzed
    % NumOfChr: Number of chromosomes
    % PriorMatrix: Store the priority details
    % NumOfPriorMatrix: Stores the number of chromosomes at ith priority
    %**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
    NumOfChr=16;
    PriorMatrix=cell(NumOfChr,NumOfPrior);
    NumOfPriorMatrix=zeros(NumOfChr,NumOfPrior);
    % Extracting the chromosome at ith priority for a chromosome
    for Chr=1:NumOfChr
        TempCol=2+(Chr-1)*4;
        for Prior=1:NumOfPrior
            TempChrNum=SortedInfoMatrix{1,1}(Prior,TempCol);
            TempNum=NumOfPriorMatrix(TempChrNum,Prior)+1;
            PriorMatrix{TempChrNum,Prior}(TempNum,1)=Chr;
            NumOfPriorMatrix(TempChrNum,Prior)=TempNum;
        end
    end
    %Total number of chromosomes upto ith (i=NumOfPrior) priority
    for Chr=1:NumOfChr
        NumOfPriorMatrix(Chr,NumOfPrior+2)=sum(NumOfPriorMatrix(Chr,1:NumOfPrior));
    end
end
