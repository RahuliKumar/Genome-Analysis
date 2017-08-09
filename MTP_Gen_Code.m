
%% To execute the whole of the following code

CodeInfo=' 1.LET_LET; 2.RET_RET; 3.RET_LET; 4.Centromere; 12.LET_LET+RET_RET; 123. LET_LET+RET_RET+Centromere and so on';
disp(CodeInfo);
Analysis=input('Enter an option: ');


%% 0.Code Description; 

% General information of code sections: 
%--------------------------------------------------------------------------
% 0.Code Description; 
% 1.Sequence Information; 
% 2.Sequence Storing; 
% 3.Sequence Alignment; 
% 4.Chromosome Number Optimization; 
% 5.Chromosome Grouping based on local alignment; 
% 6.Priority Chromosome Analysis (PCA) based on local alignment; 
%--------------------------------------------------------------------------

% 0.1**********************************************************************
% Comment containing information about the code in this code section.
% 0.1**********************************************************************

%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
% Description of important variables being created in this section and might be used in furthur code section; Variable Class
%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**

% *=*=*=*=*=*=*=*=*=Resusable variable description [RVD]*=*=*=*=*=*=*=*=*=*=
% RVD: Variable being used from above sub-section(s).Number denotes origin code section.
% *=*=*=*=*=*=*=*=*=Resusable variable description [RVD]*=*=*=*=*=*=*=*=*=


%% 1.Sequence Information; Status:D

%**************************************************************************
% Organism: Saccharomyces cerevisiae [S288c strain]; Chromosomes: 16 chromosomes; Data Source: GenBank [NCBI].
% GenBank Accesion Number: {'BK006935',....,'BK006949'};['BK006934' is for Chr VIII (except this, all in proper sequence)]
%**************************************************************************

%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
% 1.1 LETEnd: Stores end position of Left End Telomere (LET) for each chromosome; 
%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**

%Chromosome accession number
AccessionNum={'BK006935','BK006936','BK006937','BK006938','BK006939','BK006940','BK006941','BK006934','BK006942','BK006943','BK006944','BK006945','BK006946','BK006947','BK006948','BK006949'};
% Left End Telomere (LET) [Starting position for each Left End Telomere (LET) is 1 so its not required to be stored]
LETEnd=[801,6608,1098,904,6473,5530,781,5505,7784,7767,807,12085,6344,7428,847,7223]; 
% Right End Telomere (RET)
RETStart=[229411,812379,315783,1524625,569599,269731,1083635,556105,439068,744902,665904,1064281,923541,783278,1083922,942396]; 
RETEnd=[230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948010]; 
%Centromere 
CENStart=[151465,238207,114385,449711,151987,148510,496920,105586,355629,436307,440129,150828,268031,628758,326584,555957];
CENEnd=[151582,238323,114501,449821,152104,148627,497038,105703,355745,436425,440246,150947,268149,628875,326702,556073];

Path='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\'; % Save Folder
for ChrNum=1:16
    disp(['Saving Sequence for Chr:',num2str(ChrNum)]);
    ChrSeq=getgenbank(AccessionNum{ChrNum},'SequenceOnly',true,'ToFile',[Path,'Genome\Chromosome_',num2str(ChrNum),'.fsa']);
    LETSeq=getgenbank(AccessionNum{ChrNum},'SequenceOnly',true,'PARTIALSEQ',[1,LETEnd(ChrNum)],'ToFile',[Path,'LeftTelomere\Chromosome_',num2str(ChrNum),'.fsa']);
    RETSeq=getgenbank(AccessionNum{ChromosomeNumber},'SequenceOnly',true,'PARTIALSEQ',[RETStart(ChrNum),RETEnd(ChrNum)],'ToFile',[Path,'RightTelomere\Chromosome_',num2str(ChrNum),'.fsa']);
    CENSeq=getgenbank(AccessionNum{ChrNum},'SequenceOnly',true,'PARTIALSEQ',[CENStart(ChrNum),CENEnd(ChrNum)],'ToFile',[Path,'Centromere\Chromosome_',num2str(ChrNum),'.fsa']);
end

%% 3(a).Sequence Alignment; Status:D

%**************************************************************************
% Following for loop is for LET Local Alignment (LETLA)
%**************************************************************************

%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
% 3(a).1 AlignScoreMatrix: Matrix to store score from LETLA using swalign inbuilt MATLAB function; Data type:Array
%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
AlignScoreMatrix=zeros(16);

% */*/*/*/*/*/*/*/*/*/*/*/* [To store LETLA score]*/*/*/*/*/*/*/*/*/*/*/*/*

for ChromosomeNumberT1=1:16
    TempSeq1=fastaread(['D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\LeftTelomere\Chromosome_',num2str(ChromosomeNumberT1),'.fsa']);
    for ChromosomeNumberT2=1:16
        TempSeq2=fastaread(['D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\LeftTelomere\Chromosome_',num2str(ChromosomeNumberT2),'.fsa']);
        TempScore=swalign(TempSeq1,TempSeq2);
        AlignScoreMatrix(ChromosomeNumberT1,ChromosomeNumberT2)=TempScore;
    end
end
save('AlignScoreMatrix.mat','AlignScoreMatrix');
% */*/*/*/*/*/*/*/*/*/*/*/* [To store LETLA score]*/*/*/*/*/*/*/*/*/*/*/*/*

% **--**--**--**--**--** [RESULT FROM THIS SECTION] **--**--**--**--**--**
% RESULT: AlignScoreMatrix stored as AlignScoreMatrix.mat file
% **--**--**--**--**--** [RESULT FROM THIS SECTION] **--**--**--**--**--**


%% 3(b).
% 3.1 *'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*
% Data from above for loop has been stored and being loaded as variable A to save time.
A=load('AlignScoreMatrix.mat');
AlignScoreMatrix=A.AlignScoreMatrix;
% 3.1 *'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*'*

% */*/*/*/*/*/*/[To visualize the LETLA score in matrix form]*/*/*/*/*/*/*/
imagesc(AlignScoreMatrix);
colormap(hot);colorbar;title('LET Sequence Comparison');
xlabel('Chromosome Number');ylabel('Chromosome Number');
% */*/*/*/*/*/*/[To visualize the LETLA score in matrix form]*/*/*/*/*/*/*/

% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
% print(gcf,'comparision','-dpng','-r600');
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*

% **--**--**--**--**--** [RESULT FROM THIS SECTION] **--**--**--**--**--**
% RESULT: Heatmap generated for LETLA matrix for visualization
% **--**--**--**--**--** [RESULT FROM THIS SECTION] **--**--**--**--**--**

%% 4.Chromosome Number Optimization; Status:D

%**************************************************************************
% The following aproach is being used to select number of sequences for grouping analysis.
%**************************************************************************

% *=*=*=*=*=*=*=*=*=Resusable variable description [RVD]*=*=*=*=*=*=*=*=*=
% 3(b). AlignScoreMatrix: Matrix to store score from LETLA using swalign inbuilt MATLAB function; Data type:Array
% *=*=*=*=*=*=*=*=*=Resusable variable description [RVD]*=*=*=*=*=*=*=*=*=
A=load('AlignScoreMatrix.mat');
AlignScoreMatrix=A.AlignScoreMatrix;
close all;
%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
% SelectionMatrixI:  Stores LETLA sorted score in one way (all stored in different cells)
% SelectionMatrixII: Stores LETLA sorted score in another way (all stored in one cell)
%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
SelectionMatrixI=cell(16,1);
SelectionMatrixII=cell(1);
SubPlotNumber=1;ColumnNumber=1;

% Loop for each 16 chromosomes
for ChromosomeNumberS1=1:16
    TempMatrix=zeros(16,2); TempMatrix(:,1)=1:16;
    TempMatrix(:,2)=transpose(AlignScoreMatrix(ChromosomeNumberS1,:));
    TempMatrix(ChromosomeNumberS1,:)=[];
    TempMatrix=sortrows(TempMatrix,-2);
    SelectionMatrixI{ChromosomeNumberS1,1}=TempMatrix;
    SelectionMatrixII{1}(:,ColumnNumber:ColumnNumber+1)=TempMatrix;
    if (ChromosomeNumberS1<16)
        SelectionMatrixII{1}(:,ColumnNumber+2)=ChromosomeNumberS1+1;
    end
    ColumnNumber=ColumnNumber+3;
    TempCumSum=cumsum(TempMatrix(:,2));
    
    % If condition to create new figure after each 4 chromosome analysis as each figure have 4 chromosome analysis
%     if (rem(ChromosomeNumberS1,4)==1)
%         SubPlotNumber=1;
        figure('Name','Cumulative Sum Plot','NumberTitle','on');
%         suptitle(['Chromosome:',num2str(ChromosomeNumberS1),' - ',num2str(ChromosomeNumberS1+3)]);
%     end
    
%     subplot(2,2,SubPlotNumber);
    p=plot(TempCumSum);p.Marker='s';
     
    % for naming of chromosome on each point
    for TempName=1:15
        txt =['  ',num2str(TempMatrix(TempName,1))];
        TempY=TempCumSum(TempName);
        if (TempName>8)
            TempY=TempY-(TempY)/35;
            text(TempName-.3,TempY,txt)
        else
            text(TempName,TempY,txt)
        end
        if (ChromosomeNumberS1==1)
            if (TempName==10)
                txt={'No drastic change','in slope observed'};
                text(3,1.2*TempY,txt);
            end
            if (TempName==5)
                txt={'\leftarrow This point corresponds to','     cumulative sum of LETLAS','     for Chromosome 1 with 3,4,','     14,13,15 & 12'};
                text(6.8,TempY-TempY/40,txt);
            end
        end
        % information for chromosome 2
        if (ChromosomeNumberS1==2 && TempName==9)
            txt='Drastic Change of Slope \rightarrow';
            text(2.9,TempY+(TempY)/30,txt);
        end
    end
    title(['Chromosome:',num2str(ChromosomeNumberS1)]);
    xlabel('Number of Chromosomes');ylabel('LETLAS Cumulative Sum');
    xlim([0,16]);ylim([0,TempCumSum(15)+TempCumSum(15)/10]);
%   SubPlotNumber=SubPlotNumber+1;
     
    % Save each figure after for a set of 4 chromosomes after the subplotting
    % */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
%     if (rem(ChromosomeNumberS1,4)==0)
        %pause;
        set(gca,'box','off','tickdir','out')
        set(gcf,'PaperPositionMode','auto','Units','inches');
        pos=get(gcf,'pos');
        set(gcf,'PaperSize',[pos(3), pos(4)]);
%         print(gcf,['LET_NIce_Modified_CumSum',num2str(ChromosomeNumberS1)],'-dpng','-r600');
%     end
    % */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
end

% save('SelectionMatrixI.mat','SelectionMatrixI');
% save('SelectionMatrixII.mat','SelectionMatrixII');

% **--**--**--** [RESULT OF THE ANALYSIS FROM THIS SECTION] **--**--**--**
% RESULT: 9 is the optimum number of chromosomes to be compared for second approach
% **--**--**--** [RESULT OF THE ANALYSIS FROM THIS SECTION] **--**--**--**



%% 5.Chromosome Grouping based on LET Local Alignment (LETLA); Status:D

%**************************************************************************
% Chromosome grouping based on single best score.
%**************************************************************************

%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
% GroupedChromosome): Result storing variable from Chromosome Grouping matrix
%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**

% *=*=*=*=*=*=*=*=*=Resusable variable description [RVD]*=*=*=*=*=*=*=*=*=
% 4. SelectionMatrixI: Stores LETLA sorted score in one way
% *=*=*=*=*=*=*=*=*=Resusable variable description [RVD]*=*=*=*=*=*=*=*=*=

GroupedChromosome=zeros(8,2);

UngroupedChromosome=1:16; % Chromosomes to be grouped
TempChromosomeNumber=16;  % Number of chromosomes to be grouped
GroupedChromosomeNumber=0; % Chromosomes already grouped

% While loop searches for global maximum score from the all the local alignment
% scores and terminates if all the chromosomes have been grouped
while (TempChromosomeNumber~=0)
    TempMax=[-inf,0,0];
    
    % for loop searches for global maxima score from SelectionMatrix
    for Chromosome=1:TempChromosomeNumber
        TempChromosome=UngroupedChromosome(Chromosome);
        if (TempMax(1)<SelectionMatrixI{TempChromosome,1}(1,2))
            TempMax=[SelectionMatrixI{TempChromosome,1}(1,2),TempChromosome,SelectionMatrixI{TempChromosome,1}(1,1)];
        end
    end
    
    %if global maxima is found and the grouping of chromosome is not
    %clashing with any of already grouped chromosome, then grouping occours
    if (any(UngroupedChromosome==TempMax(3)))
        UngroupedChromosome=UngroupedChromosome(UngroupedChromosome~=TempMax(2) & UngroupedChromosome~=TempMax(3));
        GroupedChromosome(GroupedChromosomeNumber+1,1:2)=[TempMax(2),TempMax(3)];
        GroupedChromosomeNumber=GroupedChromosomeNumber+1;
        TempChromosomeNumber=TempChromosomeNumber-2;
%         disp('Ungrouped Chromosome:');disp(UngroupedChromosome);
        
        % else if there is any clash with already grouped chromosome, then the
        % loop continues by deleting the clashing row from SelectionMatrix
    else
        SelectionMatrixI{TempMax(2),1}(1,:)=[];
    end
end



%% 6(a).Priority Chromosome Analysis (PCA) for each Chromosome based on LETLA; Status:ND

%**************************************************************************
% Priority Chromosoem Plotting based of LETLA from above analysis
%**************************************************************************

clear;clc;close all;

%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
% NumberOfPriorities: Upto how many priority u want to analyse
% ChromosomeNumber: Numeral name of the chromosomes
% NumberOfChromosomes: Number of chromosomes for which the analysis is being done
% PrioirtyChart: Vriable to store the priority chart
% NumberOfPrioitiesChart: Stores the number of chromosomes at ith priority
% PriorityScoreListI: List to score priority score (Ist approach)
% Weights: Weights for respective priority
% MarkerType: Shapes of marker for different priorities
%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
close all;
NumberOfPriorities=5;NumberOfChromosomes=16;
PriorityChart=cell(NumberOfChromosomes,NumberOfPriorities);
NumberOfPrioritiesChart=zeros(NumberOfChromosomes,NumberOfPriorities);
PriorityScoreListI=cell(16,1);
PriorityScoreListII=zeros(16,1);
Weights=[0.6,0.5,0.4,0.3,0.2];
MarkerType=['o','+','*','x','s','d','^','v','>','<','p','h'];
PriorityName={'Priority:1','Priority:2','Priority:3','Priority:4','Priority:5','Priority:6','Priority:7','Priority:8','Priority:9'...
    ,'Priority:10','Priority:11','Priority:12','Priority:13','Priority:14','Priority:15'};
MeanOfPriorityScore=zeros(1);


% *=*=*=*=*=*=*=*=*=Resusable variable description [RVD]*=*=*=*=*=*=*=*=*=
% 4. SelectionMatrixII: Stores LETLA sorted score in one way
% *=*=*=*=*=*=*=*=*=Resusable variable description [RVD]*=*=*=*=*=*=*=*=*=
S=load('SelectionMatrixII.mat');
SelectionMatrixII=S.SelectionMatrixII{1,1};

% Extracting the chromosome at ith priority for a chromosome
for Chromosome=1:NumberOfChromosomes
    TempCol=1+(Chromosome-1)*3;
    for Priorities=1:NumberOfPriorities
        TempChromosomeNumber=SelectionMatrixII(Priorities,TempCol);
        TempNumber=NumberOfPrioritiesChart(TempChromosomeNumber,Priorities)+1;
        PriorityChart{TempChromosomeNumber,Priorities}(TempNumber,1)=Chromosome;
        NumberOfPrioritiesChart(TempChromosomeNumber,Priorities)=TempNumber;
    end
end

%Total number of chromosome at upto ith (i=NumberOfPriorities) priority
for Chromosome=1:NumberOfChromosomes
    NumberOfPrioritiesChart(Chromosome,NumberOfPriorities+2)=sum(NumberOfPrioritiesChart(Chromosome,1:NumberOfPriorities)); 
end

% For Plotting the 3D graph
figure('Name',['3D Priority Dot Plot:',num2str(NumberOfPriorities),' Priorities']);
for Chromosome=1:NumberOfChromosomes
    disp(Chromosome);
    for Priority=1:NumberOfPriorities
        TempNumber=NumberOfPrioritiesChart(Chromosome,Priority);
        for TempChromosomeNumber=1:TempNumber
            TempChromosomeName=PriorityChart{Chromosome,Priority}(TempChromosomeNumber,1);
            scatter3(Chromosome,TempChromosomeName,Priority,100,'O','filled'); hold on;
        end
    end
end
hold off;xlabel('Chromosome Number');ylabel('Priority Chromosome');zlabel('Priority Number');

% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
% print(gcf,['LET_3D_',num2str(NumberOfPriorities),'_Priorities'],'-dpng','-r600');
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*

% For Plotting the 2D graph
fig=figure('Name',['2D Priority Dot Plot:',num2str(NumberOfPriorities),' Priorities']);

% */*/*/*/*/*/*/*/*/*/*/[Trick for plotting legends]*/*/*/*/*/*/*/*/*/*/*/*
for Priority=1:NumberOfPriorities
    if (Priority==13)
        Plots(Priority)=scatter(20,15,70,'o','filled','LineWidth',1.5);
    elseif (Priority==14)
        Plots(Priority)=scatter(20,15,70,'s','filled','LineWidth',1.5);
    elseif (Priority==15)
        Plots(Priority)=scatter(20,15,70,'d','filled','LineWidth',1.5);
    else
        Plots(Priority)=scatter(20,15,70,MarkerType(Priority),'LineWidth',1.5); hold on;
    end
end
legend(Plots(1:Priority),PriorityName(1:Priority),'Location','northeast');
% */*/*/*/*/*/*/*/*/*/*/[Trick for plotting legends]*/*/*/*/*/*/*/*/*/*/*/*
for Chromosome=1:NumberOfChromosomes
    disp(Chromosome);
    for Priority=1:NumberOfPriorities
        TempNumber=NumberOfPrioritiesChart(Chromosome,Priority);
        if (TempNumber>0)
            X_Val=repmat(Chromosome,1,TempNumber);
            Y_Val=PriorityChart{Chromosome,Priority}(:,1);
            if (Priority==13)
                scatter(X_Val,Y_Val,70,'o','filled','LineWidth',1.5);
            elseif (Priority==14)
                scatter(X_Val,Y_Val,70,'s','filled','LineWidth',1.5);
            elseif (Priority==15)
                scatter(X_Val,Y_Val,70,'d','filled','LineWidth',1.5);
            else
                scatter(X_Val,Y_Val,70,MarkerType(Priority),'LineWidth',1.5); hold on;
            end
        end
    end
end
axis([0 24 0 17]);
% set(h,'FontSize',8);
hold off;grid on;grid minor;
xlabel('Chromosome Number');ylabel('Priority Chromosome');

% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
% print(gcf,['LET_2D_',num2str(NumberOfPriorities),'_Priorities'],'-dpng','-r600');
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*

%% For quantification of the priority [Weighted Average]
figure('Name',['Priority Score Plot: ',num2str(NumberOfPriorities),' Priorities']);
for Chromosome=1:NumberOfChromosomes
    TempScore=0;
    for Priorities=1:NumberOfPriorities
        TempScore=TempScore+(Weights(Priorities))*NumberOfPrioritiesChart(Chromosome,Priorities);
    end
    if (TempScore~=0)
        TempScore=TempScore/NumberOfPrioritiesChart(Chromosome,Priorities+2);
        PriorityScoreListI{Chromosome,1}=TempScore;
        PriorityScoreListII(Chromosome,1)=TempScore;
    else
        PriorityScoreListI{Chromosome,1}=0;
        PriorityScoreListII(Chromosome,1)=0;
    end
end
MeanOfPriorityScore(1)= mean(PriorityScoreListII);
plot(PriorityScoreListII,'o','color','r','MarkerSize',8,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0]);hold on;
xlabel('Chromosome Number'); ylabel('Priority Score');xlim([0,17])
line([0,16],[MeanOfPriorityScore(1),MeanOfPriorityScore(1)],'color','b','LineWidth',2);hold off;
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
% print(gcf,['LET_',num2str(NumberOfPriorities),'_priorities Score [Weighted Average]'],'-dpng','-r600');
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*

%% For quantification of the priority [Weighted Sum]
figure('Name',['Priority Score Plot: ',num2str(NumberOfPriorities),' Priorities']);
for Chromosome=1:NumberOfChromosomes
    TempScore=0;
    for Priorities=1:NumberOfPriorities
        TempScore=TempScore+(Weights(Priorities))*NumberOfPrioritiesChart(Chromosome,Priorities);
    end
    PriorityScoreListII(Chromosome,1)=TempScore;
end
MeanOfPriorityScore(1)=mean(PriorityScoreListII);
plot(PriorityScoreListII,'o','color','r','MarkerSize',8,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0]);hold on;
xlabel('Chromosome Number'); ylabel('Priority Score');
line([0,16],[MeanOfPriorityScore(1),MeanOfPriorityScore(1)],'color','b','LineWidth',2);hold off;

% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
% print(gcf,['LET_',num2str(NumberOfPriorities),'_priorities Score [Weighted Sum]'],'-dpng','-r600');
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*



%% 6(b).Priority Chromosome Analysis (PCA) for each Chromosome based on LETLA; Status:ND

%**************************************************************************
% Priority Chromosoem Plotting based of LETLA from above analysis
%**************************************************************************

clear;clc;close all;

%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**
% NumberOfChromosomes: Number of chromosomes for which the analysis is being done
% PrioirtyChart: Vriable to store the priority chart
% NumberOfPrioitiesChart: Stores the number of chromosomes at ith priority
% PriorityScoreListI: List to score priority score (Ist approach)
% Weights: Weights for respective priority
% MarkerType: Shapes of marker for different priorities
%**==**==**==**==**==**==**[Variable Description]**==**==**==**==**==**==**

NumberOfPriorities=15;NumberOfChromosomes=16;
PriorityChart=cell(NumberOfChromosomes,NumberOfPriorities);
NumberOfPrioritiesChart=zeros(NumberOfChromosomes,NumberOfPriorities);
PriorityScoreListI=cell(16,1);
PriorityScoreListII=zeros(16,1);
% Weights=[1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0];
% MarkerType=['o','+','*','x','s','d','^','>','<','p','h'];
MarkerType=['o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o'];
PriorityName={'Priority:1','Priority:2','Priority:3','Priority:4','Priority:5','Priority:6','Priority:7','Priority:8','Priority:9'...
    ,'Priority:10','Priority:11','Priority:12','Priority:13','Priority:14','Priority:15'};
MeanOfPriorityScore=zeros(1);


% *=*=*=*=*=*=*=*=*=Resusable variable description [RVD]*=*=*=*=*=*=*=*=*=
% 4. SelectionMatrixII: Stores LETLA sorted score in one way
% *=*=*=*=*=*=*=*=*=Resusable variable description [RVD]*=*=*=*=*=*=*=*=*=
S=load('SelectionMatrixII.mat');
SelectionMatrixII=S.SelectionMatrixII{1,1};

% Extracting the chromosome at ith priority for a chromosome
for Chromosome=1:NumberOfChromosomes
    TempCol=1+(Chromosome-1)*3;
    for Priorities=1:NumberOfPriorities
        TempChromosomeNumber=SelectionMatrixII(Priorities,TempCol);
        TempNumber=NumberOfPrioritiesChart(TempChromosomeNumber,Priorities)+1;
        PriorityChart{TempChromosomeNumber,Priorities}(TempNumber,1)=Chromosome;
        NumberOfPrioritiesChart(TempChromosomeNumber,Priorities)=TempNumber;
    end
end

%Total number of chromosome at upto ith (i=NumberOfPriorities) priority
for Chromosome=1:NumberOfChromosomes
    NumberOfPrioritiesChart(Chromosome,NumberOfPriorities+2)=sum(NumberOfPrioritiesChart(Chromosome,1:NumberOfPriorities)); 
end

% For Plotting the 3D graph
figure('Name',['3D Priority Dot Plot:',num2str(NumberOfPriorities),' Priorities']);
for Chromosome=1:NumberOfChromosomes
    disp(Chromosome);
    for Priority=1:NumberOfPriorities
        TempNumber=NumberOfPrioritiesChart(Chromosome,Priority);
        for TempChromosomeNumber=1:TempNumber
            TempChromosomeName=PriorityChart{Chromosome,Priority}(TempChromosomeNumber,1);
            scatter3(Chromosome,TempChromosomeName,Priority,50,MarkerType(Priority),'filled'); hold on;
            
        end
    end
end
hold off;xlabel('Chromosome Number');ylabel('Priority Chromosome');zlabel('Priority Number');

% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
% print(gcf,['LET_3D_',num2str(NumberOfPriorities),'_Priorities'],'-dpng','-r600');
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*

% For Plotting the 2D graph
figure('Name',['2D Priority Dot Plot:',num2str(NumberOfPriorities),' Priorities']);
for Chromosome=1:NumberOfChromosomes
    disp(Chromosome);
    for Priority=1:NumberOfPriorities
        TempNumber=NumberOfPrioritiesChart(Chromosome,Priority);
        for TempChromosomeNumber=1:TempNumber
            TempChromosomeName=PriorityChart{Chromosome,Priority}(TempChromosomeNumber,1);
            scatter(Chromosome,TempChromosomeName,40,MarkerType(Priority),'filled','MarkerFaceColor','b'); hold on;
        end
    end
end
xlabel('Chromosome Number');ylabel('Priority Chromosome');axis([0,17,0,17]);hold off;

% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
% print(gcf,['LET_2D_',num2str(NumberOfPriorities),'_Priorities'],'-dpng','-r600');
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*



%% For quantification of the priority [Weighted Average]
figure('Name',['Priority Score Plot: ',num2str(NumberOfPriorities),' Priorities']);
for Chromosome=1:NumberOfChromosomes
    TempScore=0;
    for Priorities=1:NumberOfPriorities
        Weights=(1/NumberOfPriorities)*(NumberOfPriorities+1-Priorities);
        TempScore=TempScore+(Weights)*NumberOfPrioritiesChart(Chromosome,Priorities);
    end
    if (TempScore~=0)
        TempScore=TempScore/NumberOfPrioritiesChart(Chromosome,Priorities+2);
        PriorityScoreListI{Chromosome,1}=TempScore;
        PriorityScoreListII(Chromosome,1)=TempScore;
    else
        PriorityScoreListI{Chromosome,1}=0;
        PriorityScoreListII(Chromosome,1)=0;
    end
end
MeanOfPriorityScore(1)=mean(PriorityScoreListII);
plot(PriorityScoreListII,'o','color','r','MarkerSize',8,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0]);hold on;
xlabel('Chromosome Number');ylabel('Priority Score');
line([0,16],[MeanOfPriorityScore(1),MeanOfPriorityScore(1)],'color','b','LineWidth',2);hold off;

% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
% print(gcf,['LET_',num2str(NumberOfPriorities),'_priorities Score [Weighted Average]'],'-dpng','-r600');
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*


%% For quantification of the priority [Weighted Sum]
figure('Name','Summation Priority Score Plot');
for Chromosome=1:NumberOfChromosomes
    TempScore=0;
    for Priorities=1:NumberOfPriorities
        Weights=(1/NumberOfPriorities)*(NumberOfPriorities+1-Priorities);
        TempScore=TempScore+(Weights)*NumberOfPrioritiesChart(Chromosome,Priorities);
    end
    PriorityScoreListII(Chromosome,1)=TempScore;
end
MeanOfPriorityScore(1)=mean(PriorityScoreListII);
plot(PriorityScoreListII,'o','color','r','MarkerSize',8,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0]);hold on;
xlabel('Chromosome Number');
ylabel('Priority Score');
line([0,16],[MeanOfPriorityScore(1),MeanOfPriorityScore(1)],'color','b','LineWidth',2);hold off;
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
% print(gcf,['LET_',num2str(NumberOfPriorities),'_priorities Score [Weighted Sum]'],'-dpng','-r600');
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*


% Chromosome Importance Based on Priority score  
%%

PriorityScoreListIII=zeros(16,15);
Plots=zeros(16,15);
Name={'Chromosome:1','Chromosome:2','Chromosome:3','Chromosome:4','Chromosome:5','Chromosome:6',...
    'Chromosome:7','Chromosome:8','Chromosome:9','Chromosome:10','Chromosome:11','Chromosome:12',...
    'Chromosome:13','Chromosome:14','Chromosome:15','Chromosome:16'};
hold all;xlabel('Number of Priorities');ylabel('Priority Score [Weighted Sum]');xlim([1,16]);ylim([0,14]);
WithoutBracket=[13,6,9,10,5];

for Chromosome=1:NumberOfChromosomes
    TempScore=0;
    for Priorities=1:NumberOfPriorities
        Weights=(1/NumberOfPriorities)*(NumberOfPriorities+1-Priorities);
        TempScore=TempScore+(Weights)*NumberOfPrioritiesChart(Chromosome,Priorities);
        PriorityScoreListIII(Chromosome,Priorities)=TempScore;
    end
    Plots(Chromosome)= plot(PriorityScoreListIII(Chromosome,1:Priorities),'LineWidth',2); 
    if (sum(find(WithoutBracket==Chromosome)))
        Txt = ['\rightarrow ',num2str(Chromosome)];
        text(15.15,PriorityScoreListIII(Chromosome,Priorities),Txt);
    else
        text(15.1,PriorityScoreListIII(Chromosome,Priorities),num2str(Chromosome));
    end
   
end
plot([5,5],[0,15],'--','color','b','LineWidth',1.2);hold off;
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*
set(gca,'box','off','tickdir','out')
set(gcf,'PaperPositionMode','auto','Units','inches'); pos=get(gcf,'pos');
set(gcf,'PaperSize',[pos(3), pos(4)]);
% print(gcf,'Cumulative Priorities Score ','-dpng','-r600');
% */*/*/*/*/*/*/*/*/*/*/[To save high quality image]*/*/*/*/*/*/*/*/*/*/*/*













