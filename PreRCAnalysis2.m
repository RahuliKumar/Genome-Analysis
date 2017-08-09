tic;clc;

SavePath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';

PreRC_Chr_ScoreAndIdent_Matrix=cell(16,35);
PreRC_Chr_Ident_Matrix=zeros(16,9);
ORC_ORC_ScoreAndIdent_Matrix=zeros(6,17);
ORC_ORC_Ident_Matrix=zeros(6);

% PreRCs=[ORC1,...,ORC6,CDC6,CDT1,MCM2]
PreRCStart=[142210,360652,141073,866715,155100,344321,69338,522048,174920];
PreRCEnd=[144954,362514,142923,868304,156539,345628,70879,523862,177526];
PreRCChrNum=[13,2,12,16,14,8,10,10,2];

PreRCs={'ORC1','ORC2','ORC3','ORC4','ORC5','ORC6','CDC6','CDT1','MCM2'};
SeqPath1='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\Genome\';
SeqPath2='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\PRE-RC\';

%checking for similarity in data
% for PreRC=1:9
%     TempSeq1=fastaread([SeqPath1,'\Chromosome_',num2str(PreRCChrNum(PreRC)),'.fsa']);
%     TempSeq1=TempSeq1.Sequence;
%     TempSeq1=TempSeq1(PreRCStart(PreRC):PreRCEnd(PreRC));
%     TempSeq2=fastaread([SeqPath2,'\',PreRCs{PreRC},'.fsa']);
%     [AlignScore,Alignment]=swalign(TempSeq1,TempSeq2);
%     Size=size(Alignment);
%     Identical=0;
%     for i=1:Size(2)
%         if (strcmp(Alignment(2,i),'|'))
%             Identical=Identical+1;
%         end
%     end
%     Identical=(Identical*100)/Size(2);
%     disp('****************');
%     disp([PreRC,Identical]);
%     disp(size(TempSeq1));
%     disp(size(TempSeq2.Sequence));
% end



for PreRC=1:2
    disp('************');
    disp(PreRC);
    TempSeq1=fastaread([SeqPath1,'\Chromosome_',num2str(PreRCChrNum(PreRC)),'.fsa']);
    TempSeq1=TempSeq1.Sequence;
    TempSeq1=TempSeq1(PreRCStart(PreRC):PreRCEnd(PreRC));
    %PreRC-Chr Analysis Chromosome
    for ChrNum=1:2
        disp(ChrNum);
        TempSeq2=fastaread([SeqPath1,'\Chromosome_',num2str(ChrNum),'.fsa']);
        [AlignScore,Alignment]=swalign(TempSeq1,TempSeq2);
        TempSeq2=TempSeq2.Sequence;
        Repeats=strfind(TempSeq2,TempSeq1);
        PreRC_Chr_ScoreAndIdent_Matrix{ChrNum,3+(PreRC-1)*4}=Repeats;
        PreRC_Chr_ScoreAndIdent_Matrix{ChrNum,1+(PreRC-1)*4}=AlignScore;
        
        Size=size(Alignment);
        Identical=0;
        for i=1:Size(2)
            if (strcmp(Alignment(2,i),'|'))
                Identical=Identical+1;
            end
        end
        PreRC_Chr_ScoreAndIdent_Matrix{ChrNum,2+(PreRC-1)*4}=(Identical*100)/Size(2);
        PreRC_Chr_Ident_Matrix(ChrNum,PreRC)=(Identical*100)/Size(2);
    end
    % ORC-ORC Analysis
%     if (PreRC<7)
        for ORC=1:0
            TempSeq2=fastaread([SeqPath1,'\Chromosome_',num2str(PreRCChrNum(ORC)),'.fsa']);
            TempSeq2=TempSeq2.Sequence;
            TempSeq2=TempSeq2(PreRCStart(ORC):PreRCEnd(ORC));
            [AlignScore,Alignment]=swalign(TempSeq1,TempSeq2);
            ORC_ORC_ScoreAndIdent_Matrix(PreRC,1+(ORC-1)*3)=AlignScore;
            Size=size(Alignment);
            Identical=0;
            for i=1:Size(2)
                if (strcmp(Alignment(2,i),'|'))
                    Identical=Identical+1;
                end
            end
            ORC_ORC_ScoreAndIdent_Matrix(PreRC,2+(ORC-1)*3)=(Identical*100)/Size(2);
            ORC_ORC_Ident_Matrix(PreRC,ORC)=(Identical*100)/Size(2);
        end
%     end
end
% save([SavePath,'\PreRC_Chr_ScoreAndIdent_Matrix.mat'],'PreRC_Chr_ScoreAndIdent_Matrix');
% save([SavePath,'\PreRC_Chr_Ident_Matrix.mat'],'PreRC_Chr_Ident_Matrix');
% save([SavePath,'\ORC_ORC_ScoreAndIdent_Matrix.mat'],'ORC_ORC_ScoreAndIdent_Matrix');
% save([SavePath,'\ORC_ORC_Ident_Matrix.mat'],'ORC_ORC_Ident_Matrix');
toc

% The above analysis took 2930.285335 seconds.