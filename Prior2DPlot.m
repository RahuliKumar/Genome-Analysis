function Prior2DPlot(NumOfPriorMatrix,PriorMatrix,NumOfPrior,AnalyFigName,Save)
    Plots=zeros(NumOfPrior,1);
    PriorityName={'Priority:1','Priority:2','Priority:3','Priority:4','Priority:5','Priority:6','Priority:7','Priority:8','Priority:9'...
        ,'Priority:10','Priority:11','Priority:12','Priority:13','Priority:14','Priority:15'};
    MarkerType=['o','+','*','x','s','d','^','v','>','<','p','h'];
    figure('Name',['2DPriorPlot:',num2str(NumOfPrior),' Priorities']);
    % */*/*/*/*/*/*/*/*/*/*/[Trick for plotting legends]*/*/*/*/*/*/*/*/*/*/*/*
    for Prior=1:NumOfPrior
        if (Prior==13)
            Plots(Prior)=scatter(20,15,70,'o','filled','LineWidth',1.5);
        elseif (Prior==14)
            Plots(Prior)=scatter(20,15,70,'s','filled','LineWidth',1.5);
        elseif (Prior==15)
            Plots(Prior)=scatter(20,15,70,'d','filled','LineWidth',1.5);
        else
            Plots(Prior)=scatter(20,15,70,MarkerType(Prior),'LineWidth',1.5); hold on;
        end
    end
    legend(Plots(1:Prior),PriorityName(1:Prior),'Location','northeast');
    % */*/*/*/*/*/*/*/*/*/*/[Trick for plotting legends]*/*/*/*/*/*/*/*/*/*/*/*
    for Chr=1:16
        for Prior=1:NumOfPrior
            TempNumb=NumOfPriorMatrix(Chr,Prior);
            if (TempNumb>0)
                X_Val=repmat(Chr,1,TempNumb);
                Y_Val=PriorMatrix{Chr,Prior}(:,1);
                if (Prior==13)
                    scatter(X_Val,Y_Val,70,'o','filled','LineWidth',1.5);
                elseif (Prior==14)
                    scatter(X_Val,Y_Val,70,'s','filled','LineWidth',1.5);
                elseif (Prior==15)
                    scatter(X_Val,Y_Val,70,'d','filled','LineWidth',1.5);
                else
                    scatter(X_Val,Y_Val,70,MarkerType(Prior),'LineWidth',1.5); hold on;
                end
            end
        end
    end
    axis([0 24 0 17]);hold off;grid on;grid minor;
    title([AnalyFigName,':',num2str(NumOfPrior), ' Priorities']);
    xlabel('Chromosome Number');ylabel('Priority Chromosome');
    if(strcmp(Save,'on'))
        SaveImage(['Prior2DPlot_',num2str(NumOfPrior),'_',AnalyFigName]);    
    end
end