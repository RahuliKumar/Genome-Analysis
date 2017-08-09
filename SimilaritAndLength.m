clc;
LETEnd=[801,6608,1098,904,6473,5530,781,5505,7784,7767,807,12085,6344,7428,847,7223];
% Right End Telomere (RET)
DataPath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';
SortedInfoMatrix=load([DataPath,'\SortedInfoMatrix_LET_LET']);
SortedInfoMatrix=SortedInfoMatrix.SortedInfoMatrix;
SimilarityLength=cell(1);k=1;
for i=1:16
    for j=1:15
        Sim=SortedInfoMatrix{1,1}(j,3+(i-1)*4);
%         if (Sim>=90)
            Chr2=SortedInfoMatrix{1,1}(j,2+(i-1)*4);
            Diff=abs(LETEnd(i)-LETEnd(Chr2))*100/(max(LETEnd(i),LETEnd(Chr2)));
            SimilarityLength{1,1}(k,1:6)=[i,Chr2,Sim,LETEnd(i),LETEnd(Chr2),Diff];
            k=k+1;
%         end
    end
end
SimilarityLength=sortrows(SimilarityLength{1,1},3);
A=SimilarityLength(:,3);
B=SimilarityLength(:,6);
close all;
plot(A,B,'o');
title('LET-LET Comparision');
xlabel('Sequence Similarity (%)');ylabel('Length Difference (%)');

set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
print(gcf,[DataPath,'\','LET_Similrity_Length'],'-dpng','-r600');

%%
DataPath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';
RETStart=[229411,812379,315783,1524625,569599,269731,1083635,556105,439068,744902,665904,1064281,923541,783278,1083922,942396];
RETEnd=[230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948010];

SortedInfoMatrix=load([DataPath,'\SortedInfoMatrix_RET_RET']);
SortedInfoMatrix=SortedInfoMatrix.SortedInfoMatrix;
SimilarityLength2=cell(1);k=1;
for i=1:16
    for j=1:15
        Sim=SortedInfoMatrix{1,1}(j,3+(i-1)*4);
%         if (Sim>=90)
            Chr1Length=RETEnd(i)-RETStart(i)+1;
            Chr2=SortedInfoMatrix{1,1}(j,2+(i-1)*4);
            Chr2Length=RETEnd(Chr2)-RETStart(Chr2)+1;
            Diff=abs(Chr1Length-Chr2Length)*100/(max(Chr2Length,Chr1Length));
            SimilarityLength2{1,1}(k,1:6)=[i,Chr2,Sim,Chr1Length,Chr2Length,Diff];
            k=k+1;
%         end
    end
end
SimilarityLength2=sortrows(SimilarityLength2{1,1},3);
A=SimilarityLength2(:,3);
B=SimilarityLength2(:,6);
figure;
plot(A,B,'o');
title('RET-RET Comparision');
xlabel('Sequence Similarity (%)');ylabel('Length Difference (%)');
set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
print(gcf,[DataPath,'\','RET_Similrity_Length'],'-dpng','-r600');

