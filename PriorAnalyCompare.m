function [NumOfPriorMatrix,PriorMatrix]=PriorAnalyCompare(SortedInfoMatrixI,SortedInfoMatrixII,NumOfPrior)
    
    [~,PriorMatrixI]=PriorAnaly(SortedInfoMatrixI,NumOfPrior);
    [~,PriorMatrixII]=PriorAnaly(SortedInfoMatrixII,NumOfPrior);
    
    NumOfChr=16;PriorMatrix=cell(NumOfChr,NumOfPrior);
    NumOfPriorMatrix=zeros(NumOfChr,NumOfPrior);
    % Extracting the chromosome at ith priority for a chromosome
    for Chr=1:NumOfChr
        for Prior=1:NumOfPrior
            PriorMatrix{Chr,Prior}=intersect(PriorMatrixI{Chr,Prior},PriorMatrixII{Chr,Prior});
            if (~isempty(PriorMatrix{Chr,Prior}))
                Size=size(PriorMatrix{Chr,Prior});
                NumOfPriorMatrix(Chr,Prior)=Size(1);
            end
        end
    end
    %Total number of chromosomes upto ith (i=NumOfPrior) priority
    for Chr=1:NumOfChr
        NumOfPriorMatrix(Chr,NumOfPrior+2)=sum(NumOfPriorMatrix(Chr,1:NumOfPrior));
    end
end