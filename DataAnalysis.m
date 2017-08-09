tic;clc;
Priorities=[5,9,15]; SaveFig='on'; SaveData='off';

% SeqAnalysis(Seq1,Seq2,Priorities,SaveFig,DataSave);
% close all;SeqAnalysis('LET','LET',Priorities,SaveFig,SaveData);
% close all;SeqAnalysis('RET','RET',Priorities,SaveFig,SaveData);
% close all;SeqAnalysis('CEN','CEN',Priorities,SaveFig,SaveData);
% close all;SeqAnalysis('RET','LET',Priorities,SaveFig,SaveData);

% SeqAnalysisCompareII(Analy1,Analy2,Priorities,SaveFig)
close all;SeqAnalysisCompareII('LET_LET','RET_RET',Priorities,SaveFig);
close all;SeqAnalysisCompareII('LET_LET','CEN_CEN',Priorities,SaveFig);
close all;SeqAnalysisCompareII('RET_RET','CEN_CEN',Priorities,SaveFig);

% SeqAnalysisCompare(Analy1,Analy2,Analy3,SaveFig)
close all;SeqAnalysisCompare('LET_LET','RET_RET','CEN_CEN',SaveFig);
% close all;SeqAnalysisCompare('LET_LET','RET_LET','CEN_CEN',Priorities,SaveFig);
% close all;SeqAnalysisCompare('RET_RET','RET_LET','CEN_CEN',Priorities,SaveFig);
close all;
toc;

%All of above analysis takes about 717.423856 seconds [with SaveFig='on' and SaveData='off']