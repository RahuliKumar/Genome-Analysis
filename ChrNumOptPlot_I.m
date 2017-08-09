function ChrNumOptPlot_I(SortedInfoMatrix,LegendName,Save)
    for Chr=1:16
        figure;
        Norm_Matrix=normc(SortedInfoMatrix{1,1}(:,4+(Chr-1)*4));
        plot(Norm_Matrix,'Marker','o','LineWidth',1.5);hold on;
        for ChrName=1:15
            TempNum=Norm_Matrix(ChrName);
            text(ChrName+0.2,TempNum-TempNum/50,num2str(SortedInfoMatrix{1,1}(ChrName,2+(Chr-1)*4)));
        end
        xlim([0,17]);title(['Chromosome:',num2str(Chr)]);legend(LegendName,'Location','northwest');
        xlabel('Number Of Chromosomes');ylabel('Cumulative Score[Normalized]');
        if (strcmp(Save,'on'))
            SaveImage(['ChrNumOpt_',LegendName,'_',num2str(Chr)]);    
        end
    end
end
