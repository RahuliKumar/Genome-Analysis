function GroupingPlot_II(GroupedChrI,GroupedChrII,MarkerI,MarkerII,LineStyle,ColorI,ColorII,AnalyFigNameI,AnalyFigNameII,SaveFig)
    hold all;
    for Grp=1:8
        plot(GroupedChrI(Grp,:),[Grp,Grp],LineStyle,'Marker',MarkerI,'LineWidth',1.3,'Color',ColorI,'MarkerSize',12);
        plot(GroupedChrII(Grp,:),[Grp,Grp],LineStyle,'Marker',MarkerII,'LineWidth',1.3,'Color',ColorII,'MarkerSize',12);
    end
    grid on;legend({AnalyFigNameI,AnalyFigNameII});xlim([0,18]);ylim([0,9]);xlabel('Chromosome');ylabel('Pairwise Group');title('Pairwise Grouping');hold off;
    if (strcmp(SaveFig,'on'))
        SaveImage(['GroupingPlot_II_',AnalyFigNameI,'_',AnalyFigNameII]);
    end
end