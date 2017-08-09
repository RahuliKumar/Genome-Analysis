function ChrNumOptPlot_II(SortedInfoMatrixI,SortedInfoMatrixII,SortedInfoMatrixIII,LegendName,Save)
    for Chr=1:16
        figure;
        Norm_Matrix=normc([SortedInfoMatrixI{1,1}(:,4+(Chr-1)*4),SortedInfoMatrixII{1,1}(:,4+(Chr-1)*4),SortedInfoMatrixIII{1,1}(:,4+(Chr-1)*4)]);
        plot(Norm_Matrix(:,1),'Marker','o','LineWidth',1.5);hold on;
        plot(Norm_Matrix(:,2),'Marker','s','LineWidth',1.5);
        plot(Norm_Matrix(:,3),'Marker','p','LineWidth',1.5);
        xlim([0,17]);title(['Chromosome:',num2str(Chr)]);legend(LegendName,'Location','northwest');
        xlabel('Number Of Chromosomes');ylabel('Cumulative Score [Normalized]');
        if (strcmp(Save,'on'))
            SaveImage(['ChrNumOpt_',LegendName{1},'_',LegendName{2},'_',LegendName{3},'_',num2str(Chr)]);    
        end
    end
end
