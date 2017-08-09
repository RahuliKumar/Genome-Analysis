function SaveImage(FigName)
    SavePath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\Analysis Figure Ident';
    set(gca,'box','off','tickdir','out')
    set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
    set(gcf,'PaperSize',[pos(3), pos(4)]);
    print(gcf,[SavePath,'\',FigName],'-dpng','-r600');
end
