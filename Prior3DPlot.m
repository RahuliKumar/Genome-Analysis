function Prior3DPlot(NumOfPriorMatrix,PriorMatrix,NumOfPrior,AnalyFigName,Save)
    figure('Name',['3DPriorPlot:',num2str(NumOfPrior),' Priorities']);
    
    for Chr=1:16
        for Priority=1:NumOfPrior
            TempNumb=NumOfPriorMatrix(Chr,Priority);
            for TempChrNum=1:TempNumb
                TempChrName=PriorMatrix{Chr,Priority}(TempChrNum,1);
                scatter3(Chr,TempChrName,Priority,100,'O','filled'); hold on;
            end
        end
    end
    title([AnalyFigName,':',num2str(NumOfPrior), ' Priorities']);
    hold off;xlabel('Chromosome Number');ylabel('Priority Chromosome');zlabel('Priority Number');
    if(strcmp(Save,'on'))
        SaveImage(['Prior3DPlot_',num2str(NumOfPrior),'_',AnalyFigName]);
    end
end