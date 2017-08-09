function AlignMatrix=LocalSeqAlign(Seq1,Seq2,SeqPath,SavePath)
    AlignMatrix=zeros(16);AlignIdentMatrix=zeros(16);
    SeqPath1=[SeqPath,'\',Seq1];SeqPath2=[SeqPath,'\',Seq2];
    for ChrNum1=1:16
        TempSeq1=fastaread([SeqPath1,'\Chromosome_',num2str(ChrNum1),'.fsa']);
        for ChrNum2=1:16
            TempSeq2=fastaread([SeqPath2,'\Chromosome_',num2str(ChrNum2),'.fsa']);
            [Score,Alignment]=swalign(TempSeq1,TempSeq2);
            AlignMatrix(ChrNum1,ChrNum2)=Score;
        end
        Size=size(Alignment);
        Identical=0;
        for i=1:Size(2)
            if (strcmp(Alignment(2,i),'|'))
                Identical=Identical+1;
            end
        end
        AlignIdentMatrix(ChrNum1,ChrNum2)=(Identical*100)/Size(2);
    end
    save([SavePath,'\AlignScoreMatrix_',Seq1,'_',Seq2,'.mat'],'AlignScoreMatrix');
    save([SavePath,'\AlignIdentMatrix',Seq1,'_',Seq2,'.mat'],'AlignIdentMatrix');
end