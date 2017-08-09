tic
SavePath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';

PreRC_Chr_ScoreAndIdent_Matrix=zeros(16,26);
PreRC_Chr_Ident_Matrix=zeros(16,9);
ORC_ORC_ScoreAndIdent_Matrix=zeros(6,17);
ORC_ORC_Ident_Matrix=zeros(6);
PreRCs={'ORC1','ORC2','ORC3','ORC4','ORC5','ORC6','CDC6','CDT1','MCM2'};
SeqPath1='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\Genome\';
SeqPath2='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\PRE-RC\';
for PreRC=1:9
    TempSeq1=fastaread([SeqPath2,'\',PreRCs{PreRC},'.fsa']);
    %PreRC-Chr Analysis Chromosome
    for ChrNum=1:16
        TempSeq2=fastaread([SeqPath1,'\Chromosome_',num2str(ChrNum),'.fsa']);
        [AlignScore,Alignment]=swalign(TempSeq1,TempSeq2);
        PreRC_Chr_ScoreAndIdent_Matrix(ChrNum,1+(PreRC-1)*3)=AlignScore;
        Size=size(Alignment);
        Identical=0;
        for i=1:Size(2)
            if (strcmp(Alignment(2,i),'|'))
                Identical=Identical+1;
            end
        end
        PreRC_Chr_ScoreAndIdent_Matrix(ChrNum,2+(PreRC-1)*3)=(Identical*100)/Size(2);
        PreRC_Chr_Ident_Matrix(ChrNum,PreRC)=(Identical*100)/Size(2);
    end
    % ORC-ORC Analysis
    if (PreRC<7)
        for ORC=1:6
            TempSeq2=fastaread([SeqPath2,'\',PreRCs{ORC},'.fsa']);
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
    end
end
save([SavePath,'\PreRC_Chr_ScoreAndIdent_Matrix.mat'],'PreRC_Chr_ScoreAndIdent_Matrix');
save([SavePath,'\PreRC_Chr_Ident_Matrix.mat'],'PreRC_Chr_Ident_Matrix');
save([SavePath,'\ORC_ORC_ScoreAndIdent_Matrix.mat'],'ORC_ORC_ScoreAndIdent_Matrix');
save([SavePath,'\ORC_ORC_Ident_Matrix.mat'],'ORC_ORC_Ident_Matrix');
toc

% The above analysis took 2885.398993 seconds.