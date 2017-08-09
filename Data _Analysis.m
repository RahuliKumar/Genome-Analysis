tic;
Priorities=[5,9,15]; SaveFig='on'; SaveData='off';

% TELOMERE ANALYSIS

%SeqAnalysis(Seq1,Seq2,SaveFig,SaveData);
close all;clear;
SeqAnalysis('RET','LET',Priorities,SaveFig,SaveData);
SeqAnalysis('RET','RET',Priorities,SaveFig,SaveData);
SeqAnalysis('RET','LET',Priorities,SaveFig,SaveData);

SeqAnalysisCompareII(Analy1,Analy2,SaveFig)
clc;close all;SeqAnalysisCompare('LET_LET','RET_RET',Priorities,SaveFig);
SeqAnalysisCompare('LET_LET','RET_LET',Priorities,SaveFig);
SeqAnalysisCompare('RET_RET','RET_LET',Priorities,SaveFig);

%CENTROMERE ANALYSIS
SeqAnalysis('CEN','CEN','off',Priorities,SaveFig,SaveData);

toc;

