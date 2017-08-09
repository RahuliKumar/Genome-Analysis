function GroupingPlot_I(GroupedChr,Marker,LineStyle,Color,AnalyFigName,Save)
    for Grp=1:8
        plot(GroupedChr(Grp,:),[Grp,Grp],LineStyle,'Marker',Marker,'LineWidth',1,'Color',Color,'MarkerFaceColor',Color);hold on;
    end
    title([AnalyFigName,' Pairwise Group']);xlim([0,17]);ylim([0,9]);xlabel('Chromosome');ylabel('Pairwise Group');
    if (strcmp(Save,'on'))
        SaveImage(['GroupingPlot_I_',AnalyFigName]);
    end
end