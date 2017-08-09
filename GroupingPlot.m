function GroupingPlotI(GroupedChr,Marker,LineStyle,Color)
    for Grp=1:8
        plot(GroupedChr(Grp,:),[Grp,Grp],LineStyle,'Marker',Marker,'LineWidth',1.3,'Color',Color,'MarkerFaceColor',Color);hold on;
    end
    xlim([0,17]);ylim([0,9]);xlabel('Chromosome');ylabel('Pairwise Group');
end