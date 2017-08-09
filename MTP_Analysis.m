
% TELOMERE ANALYSIS

DataPath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';SavePath=DataPath;

%% LET_LET Analysis

AnalyDataName='LET_LET';AnalyFigName='LET-LET';

ASM_LET_LET=load([DataPath,'\AlignScoreMatrix_LET_LET.mat']); % ASM: Alignment Score Matrix
ASM_LET_LET=ASM_LET_LET.AlignScoreMatrixLETLET;

imagesc(ASM_LET_LET);
colormap(hot);colorbar;title('LET-LET Sequence Comparison');
xlabel('Chromosome Number');ylabel('Chromosome Number');

% Pairwise Chromosome Grouping
GroupedChr=load([DataPath,'\ChrGrouping_LET_LET.mat']);
figure('Name',['Pairwise Grouping:',AnalyFigName]);
Marker='o';LineStyle='--';
GroupingPlot(GroupedChr,Marker,LineStyle)


ChrNumOpt(ASM_LET_LET,SavePath,AnalyDataName)




%% RET_RET Analysis
A=load([DataPath,'\AlignScoreMatrixRETRET.mat']);
AlignScoreMatrixRETRET=A.AlignScoreMatrixRETRET;

imagesc(AlignScoreMatrixRETRET);
colormap(hot);colorbar;title('RET_RET Sequence Comparison');
xlabel('Chromosome Number');ylabel('Chromosome Number');


%% LET_RET Analysis

A=load([DataPath,'\AlignScoreMatrixLETRET.mat']);
AlignScoreMatrixLETRET=A.AlignScoreMatrixLETRET;

imagesc(AlignScoreMatrixLETRET);
colormap(hot);colorbar;title('LET_RET Sequence Comparison');
xlabel('Chromosome Number');ylabel('Chromosome Number');

%% LET_LET & RET_RET Comparision 


%%  LET_LET & LET_RET Comparision 

%%  RET_LET & RET_RET Comparision 