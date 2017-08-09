AlignScoreMatrix=zeros(16);AlignIdentMatrix=zeros(16);
SeqPath1='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\RET\';
SeqPath2=SeqPath1;
for ChrNum1=10:10
    TempSeq1=fastaread([SeqPath1,'\Chromosome_',num2str(ChrNum1),'.fsa']);
    for ChrNum2=16:16
        TempSeq2=fastaread([SeqPath2,'\Chromosome_',num2str(ChrNum2),'.fsa']);
        [AlignScore,Alignment]=swalign(TempSeq1,TempSeq2);
        AlignScoreMatrix(ChrNum1,ChrNum2)=AlignScore;
        Size=size(Alignment);
        Identical=0;
        for i=1:Size(2)
            if (strcmp(Alignment(2,i),'|'))
                Identical=Identical+1;
            end
        end
        AlignIdentMatrix(ChrNum1,ChrNum2)=(Identical*100)/Size(2);
    end
end
% save(['D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data','\AlignIdentMatrix_RET_RET.mat'],'AlignIdentMatrix');
