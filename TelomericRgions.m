X_element_L_Start=[337,6139,640,435,6279,5068,318,5052,7322,7305,347,11626,5875,6964,389,6778];
X_element_L_End=[801,6608,1098,904,6473,5530,781,5505,7784,7767,807,12085,6344,7428,847,7223];

X_element_Comb_Repeat_L_Start=[63,5849,361,155,0,4820,35,4952,6949,6932,68,11359,5586,6692,118,6525];
X_element_Comb_Repeat_L_End=[336,6138,639,434,0,5067,317,5051,7321,7304,346,11625,5874,6963,388,6777];

Y_Prime_element_L_Start=[0,1,0,0,1,1,0,35,78,61,0,76,52,1,0,1];
Y_Prime_element_L_End=[0,5848,0,0,6278,4684,0,4951,6948,6931,0,5661,5537,6617,0,6524];

Y_Prime_element_XII_L_Start_End=[5724,11196];
*****************************************
X_element_R_Start=[229411,812379,315783,1524625,569599,269731,1083635,556105,439068,744902,665904,1064281,923541,783278,1083922,942396];
X_element_R_End=[229871,812848,316239,1525089,570060,270112,1084092,556575,439533,745361,666364,1064747,924005,783740,1084367,942771];

X_element_Comb_Repeat_R_Start=[0,812849,316240,1525090,570061,0,1084093,556576,439534,745362,666365,1064748,924006,783741,1084368,0];
X_element_Comb_Repeat_R_End=[0,813137,316521,1525370,570329,0,1084346,556831,439815,745662,666608,1065030,924306,784037,1084620,0];

Y_Prime_element_R_Start=[0,0,0,1525467,570330,0,1084347,556986,0,0,0,1065156,0,0,1084621,942810];
Y_Prime_element_R_End=[0,0,0,1531933,576874,0,1090940,562456,0,0,0,1071645,0,0,1091273,948010];

Y_Prime_element_XII_R_Start_End=[1071654,1078177];

AccessionNum={'BK006935','BK006936','BK006937','BK006938','BK006939','BK006940','BK006941','BK006934','BK006942','BK006943','BK006944','BK006945','BK006946','BK006947','BK006948','BK006949'};

Path='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\';
for Chr=1:2
    disp(Chr);disp('*********');
    getgenbank(AccessionNum{Chr},'SequenceOnly',true,'PARTIALSEQ',[X_element_L_Start(Chr),X_element_L_End(Chr)],'ToFile',[Path,'LET_X\Chromosome_',num2str(Chr),'.fsa']);
    getgenbank(AccessionNum{Chr},'SequenceOnly',true,'PARTIALSEQ',[X_element_R_Start(Chr),X_element_R_End(Chr)],'ToFile',[Path,'RET_X\Chromosome_',num2str(Chr),'.fsa']);
    if (X_element_Comb_Repeat_L_Start(Chr)~=0)
        getgenbank(AccessionNum{Chr},'SequenceOnly',true,'PARTIALSEQ',[X_element_Comb_Repeat_L_Start(Chr),X_element_Comb_Repeat_L_End(Chr)],'ToFile',[Path,'LET_XCR\Chromosome_',num2str(Chr),'.fsa']);
    end
    if (X_element_Comb_Repeat_R_Start(Chr)~=0)
        getgenbank(AccessionNum{Chr},'SequenceOnly',true,'PARTIALSEQ',[X_element_Comb_Repeat_R_Start(Chr),X_element_Comb_Repeat_R_End(Chr)],'ToFile',[Path,'RET_XCR\Chromosome_',num2str(Chr),'.fsa']);
    end
    if (Y_Prime_element_L_Start(Chr)~=0)
        getgenbank(AccessionNum{Chr},'SequenceOnly',true,'PARTIALSEQ',[Y_Prime_element_L_Start(Chr),Y_Prime_element_L_End(Chr)],'ToFile',[Path,'LET_Y\Chromosome_',num2str(Chr),'.fsa']);
    end
    if (Y_Prime_element_R_Start(Chr)~=0)
        getgenbank(AccessionNum{Chr},'SequenceOnly',true,'PARTIALSEQ',[Y_Prime_element_R_Start(Chr),Y_Prime_element_R_End(Chr)],'ToFile',[Path,'RET_Y\Chromosome_',num2str(Chr),'.fsa']);
    end
end
%
getgenbank(AccessionNum{12},'SequenceOnly',true,'PARTIALSEQ',[1071654,1078177],'ToFile',[Path,'RET_Y\Chromosome_12_2.fsa']);
getgenbank(AccessionNum{12},'SequenceOnly',true,'PARTIALSEQ',[5724,11196],'ToFile',[Path,'LET_Y\Chromosome_12_2.fsa']);
% Local Sequence Alignment LET_X
SeqPath1='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\LET_X';
SeqPath2=SeqPath1;
AlignMatrix_LET_X=zeros(16);AlignIdentMatrix_LET_X=zeros(16);NormAlignIdentMatrix_LET_X=zeros(16);
for ChrNum1=1:16
    TempSeq1=fastaread([SeqPath1,'\Chromosome_',num2str(ChrNum1),'.fsa']);
    for ChrNum2=1:16
        TempSeq2=fastaread([SeqPath2,'\Chromosome_',num2str(ChrNum2),'.fsa']);
        [Score,Alignment]=swalign(TempSeq1,TempSeq2);
        AlignMatrix_LET_X(ChrNum1,ChrNum2)=Score;
        
        Size=size(Alignment);
        Identical=0;
        for i=1:Size(2)
            if (strcmp(Alignment(2,i),'|'))
                Identical=Identical+1;
            end
        end
        AlignIdentMatrix_LET_X(ChrNum1,ChrNum2)=(Identical*100)/Size(2);
        TempSize=min(X_element_L_End(ChrNum1)-X_element_L_Start(ChrNum1),X_element_L_End(ChrNum2)-X_element_L_Start(ChrNum2));
        NormAlignIdentMatrix_LET_X(ChrNum1,ChrNum2)=AlignIdentMatrix_LET_X(ChrNum1,ChrNum2)*(TempSize/469);
    end
end
% LET_X Sorting
AlignIdentMatrix_LET_X_Sorted=zeros(1,3);
Row=1;
format shortg;
for Align=1:15
    for Align2=Align+1:16
        AlignIdentMatrix_LET_X_Sorted(Row,1:3)=[Align,Align2,round(AlignIdentMatrix_LET_X(Align,Align2),2)];
        Row=Row+1;
    end
end
AlignIdentMatrix_LET_X_Sorted=sortrows(AlignIdentMatrix_LET_X_Sorted,-3);
% 90% cutoff LET_X
Row=1;AlignIdent_LET_X_Sorted_90=zeros(1,3);
for i=1:120
    if (AlignIdentMatrix_LET_X_Sorted(i,3)>90)
        AlignIdent_LET_X_Sorted_90(Row,:)=AlignIdentMatrix_LET_X_Sorted(i,1:3);
        Row=Row+1;
    end
end
AlignIdent_LET_X_Sorted_90(:,5:7)=sortrows(AlignIdent_LET_X_Sorted_90(:,1:3),1);
% Local Sequence Alignment LET_Y
SeqPath1='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\LET_Y';
SeqPath2=SeqPath1;
AlignMatrix_LET_Y=zeros(16);AlignIdentMatrix_LET_Y=zeros(16);NormAlignIdentMatrix_LET_Y=zeros(16);
Chr=[2,5,6,8,9,10,12,13,14,16];
for ChrNum1=1:10
    disp(ChrNum1);
    TempSeq1=fastaread([SeqPath1,'\Chromosome_',num2str(Chr(ChrNum1)),'.fsa']);
    for ChrNum2=1:10
        disp(ChrNum2);disp('*****');
        TempSeq2=fastaread([SeqPath2,'\Chromosome_',num2str(Chr(ChrNum2)),'.fsa']);
        [Score,Alignment]=swalign(TempSeq1,TempSeq2);
        AlignMatrix_LET_Y(Chr(ChrNum1),Chr(ChrNum2))=Score;
        
        Size=size(Alignment);
        Identical=0;
        for i=1:Size(2)
            if (strcmp(Alignment(2,i),'|'))
                Identical=Identical+1;
            end
        end
        AlignIdentMatrix_LET_Y(Chr(ChrNum1),Chr(ChrNum2))=(Identical*100)/Size(2);
        NormAlignIdentMatrix_LET_Y(Chr(ChrNum1),Chr(ChrNum2))=AlignIdentMatrix_LET_Y(Chr(ChrNum1),Chr(ChrNum2))*(Size(2)/6870);
    end
end

% LET_Y Sorting
AlignIdentMatrix_LET_Y_Sorted=zeros(1,3);
Row=1;
for Align=1:15
    for Align2=Align+1:16
        AlignIdentMatrix_LET_Y_Sorted(Row,1:3)=[Align,Align2,AlignIdentMatrix_LET_Y(Align,Align2)];
        Row=Row+1;
    end
end
AlignIdentMatrix_LET_Y_Sorted=sortrows(AlignIdentMatrix_LET_Y_Sorted,-3);
% 90% cutoff LET_Y
Row=1;AlignIdent_LET_Y_Sorted_90=zeros(1,3);
for i=1:120
    if (AlignIdentMatrix_LET_Y_Sorted(i,3)>90)
        AlignIdent_LET_Y_Sorted_90(Row,:)=AlignIdentMatrix_LET_Y_Sorted(i,1:3);
        Row=Row+1;
    end
end
%
AlignIdent_LET_Y_Sorted_90(:,5:7)=sortrows(AlignIdent_LET_Y_Sorted_90(:,1:3),1);
% LET_XCR 
SeqPath1='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\LET_XCR';
SeqPath2=SeqPath1;
AlignMatrix_LET_XCR=zeros(16);AlignIdentMatrix_LET_XCR=zeros(16);
Chr=[1,2,3,4,6,7,8,9,10,11,12,13,14,15,16];
for ChrNum1=1:15
    disp(ChrNum1);
    TempSeq1=fastaread([SeqPath1,'\Chromosome_',num2str(Chr(ChrNum1)),'.fsa']);
    for ChrNum2=1:15
        disp(ChrNum2);disp('*****');
        TempSeq2=fastaread([SeqPath2,'\Chromosome_',num2str(Chr(ChrNum2)),'.fsa']);
        [Score,Alignment]=swalign(TempSeq1,TempSeq2);
        AlignMatrix_LET_XCR(Chr(ChrNum1),Chr(ChrNum2))=Score;
        
        Size=size(Alignment);
        Identical=0;
        for i=1:Size(2)
            if (strcmp(Alignment(2,i),'|'))
                Identical=Identical+1;
            end
        end
        AlignIdentMatrix_LET_XCR(Chr(ChrNum1),Chr(ChrNum2))=(Identical*100)/Size(2);
    end
end
% LET_XCR Sorting

AlignIdentMatrix_LET_XCR_Sorted=zeros(1,3);
Row=1;
for Align=1:15
    for Align2=Align+1:16
        AlignIdentMatrix_LET_XCR_Sorted(Row,1:3)=[Align,Align2,AlignIdentMatrix_LET_XCR(Align,Align2)];
        Row=Row+1;
    end
end
AlignIdentMatrix_LET_XCR_Sorted=sortrows(AlignIdentMatrix_LET_XCR_Sorted,-3);

% 90% cutoff LET_XCR
Row=1;AlignIdent_LET_XCR_Sorted_90=zeros(1,3);
for i=1:120
    if (AlignIdentMatrix_LET_XCR_Sorted(i,3)>90)
        AlignIdent_LET_XCR_Sorted_90(Row,:)=AlignIdentMatrix_LET_XCR_Sorted(i,1:3);
        Row=Row+1;
    end
end
% Local Sequence Alignment RET_X
SeqPath1='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\RET_X';
SeqPath2=SeqPath1;
AlignMatrix_RET_X=zeros(16);AlignIdentMatrix_RET_X=zeros(16);NormAlignIdentMatrix_RET_X=zeros(16);
for ChrNum1=1:16
    TempSeq1=fastaread([SeqPath1,'\Chromosome_',num2str(ChrNum1),'.fsa']);
    for ChrNum2=1:16
        TempSeq2=fastaread([SeqPath2,'\Chromosome_',num2str(ChrNum2),'.fsa']);
        [Score,Alignment]=swalign(TempSeq1,TempSeq2);
        AlignMatrix_RET_X(ChrNum1,ChrNum2)=Score;
        
        Size=size(Alignment);
        Identical=0;
        for i=1:Size(2)
            if (strcmp(Alignment(2,i),'|'))
                Identical=Identical+1;
            end
        end
        AlignIdentMatrix_RET_X(ChrNum1,ChrNum2)=(Identical*100)/Size(2);
        NormAlignIdentMatrix_RET_X(ChrNum1,ChrNum2)=AlignIdentMatrix_RET_X(ChrNum1,ChrNum2)*(Size(2)/470);

    end
end
% RET_X Sorting
AlignIdentMatrix_RET_X_Sorted=zeros(1,3);
Row=1;
for Align=1:15
    for Align2=Align+1:16
        AlignIdentMatrix_RET_X_Sorted(Row,1:3)=[Align,Align2,AlignIdentMatrix_RET_X(Align,Align2)];
        Row=Row+1;
    end
end

AlignIdentMatrix_RET_X_Sorted=sortrows(AlignIdentMatrix_RET_X_Sorted,-3);
% 90% cutoff RET_X
Row=1;AlignIdent_RET_X_Sorted_90=zeros(1,3);
for i=1:120
    if (AlignIdentMatrix_RET_X_Sorted(i,3)>90)
        AlignIdent_RET_X_Sorted_90(Row,:)=AlignIdentMatrix_RET_X_Sorted(i,1:3);
        Row=Row+1;
    end
end
% Local Sequence Alignment RET_Y
SeqPath1='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\RET_Y';
SeqPath2=SeqPath1;

AlignMatrix_RET_Y=zeros(16);AlignIdentMatrix_RET_Y=zeros(16);NormAlignIdentMatrix_RET_Y=zeros(16);
Chr=[4,5,7,8,12,15,16];
for ChrNum1=1:7
    disp(ChrNum1);
    TempSeq1=fastaread([SeqPath1,'\Chromosome_',num2str(Chr(ChrNum1)),'.fsa']);
    for ChrNum2=1:7
        disp(ChrNum2);disp('*****');
        TempSeq2=fastaread([SeqPath2,'\Chromosome_',num2str(Chr(ChrNum2)),'.fsa']);
        [Score,Alignment]=swalign(TempSeq1,TempSeq2);
        AlignMatrix_RET_Y(Chr(ChrNum1),Chr(ChrNum2))=Score;
        
        Size=size(Alignment);
        Identical=0;
        for i=1:Size(2)
            if (strcmp(Alignment(2,i),'|'))
                Identical=Identical+1;
            end
        end
        AlignIdentMatrix_RET_Y(Chr(ChrNum1),Chr(ChrNum2))=(Identical*100)/Size(2);
        NormAlignIdentMatrix_RET_Y(Chr(ChrNum1),Chr(ChrNum2))=AlignIdentMatrix_RET_Y(Chr(ChrNum1),Chr(ChrNum2))*(Size(2)/6652);

        
    end
end
% RET_Y Sorting
AlignIdentMatrix_RET_Y_Sorted=zeros(1,3);
Row=1;
for Align=1:15
    for Align2=Align+1:16
        AlignIdentMatrix_RET_Y_Sorted(Row,1:3)=[Align,Align2,AlignIdentMatrix_RET_Y(Align,Align2)];
        Row=Row+1;
    end
end
AlignIdentMatrix_RET_Y_Sorted=sortrows(AlignIdentMatrix_RET_Y_Sorted,-3);

% 90% cutoff RET_Y
Row=1;AlignIdent_RET_Y_Sorted_90=zeros(1,3);
for i=1:120
    if (AlignIdentMatrix_RET_Y_Sorted(i,3)>90)
        AlignIdent_RET_Y_Sorted_90(Row,:)=AlignIdentMatrix_RET_Y_Sorted(i,1:3);
        Row=Row+1;
    end
end

% RET_XCR 
SeqPath1='D:\IITD\All\Academics\SEM 9\BED851 (Project)\S.Cerevisiae Sequence\RET_XCR';
SeqPath2=SeqPath1;
AlignMatrix_RET_XCR=zeros(16);AlignIdentMatrix_RET_XCR=zeros(16);
Chr=[2,3,4,5,7,8,9,10,11,12,13,14,15];
for ChrNum1=1:13
    disp(ChrNum1);
    TempSeq1=fastaread([SeqPath1,'\Chromosome_',num2str(Chr(ChrNum1)),'.fsa']);
    for ChrNum2=1:13
        disp(ChrNum2);disp('*****');
        TempSeq2=fastaread([SeqPath2,'\Chromosome_',num2str(Chr(ChrNum2)),'.fsa']);
        [Score,Alignment]=swalign(TempSeq1,TempSeq2);
        AlignMatrix_RET_XCR(Chr(ChrNum1),Chr(ChrNum2))=Score;
        
        Size=size(Alignment);
        Identical=0;
        for i=1:Size(2)
            if (strcmp(Alignment(2,i),'|'))
                Identical=Identical+1;
            end
        end
        AlignIdentMatrix_RET_XCR(Chr(ChrNum1),Chr(ChrNum2))=(Identical*100)/Size(2);
    end
end
% RET_XCR Sorting
AlignIdentMatrix_RET_XCR_Sorted=zeros(1,3);
Row=1;
for Align=1:15
    for Align2=Align+1:16
        AlignIdentMatrix_RET_XCR_Sorted(Row,1:3)=[Align,Align2,AlignIdentMatrix_RET_XCR(Align,Align2)];
        Row=Row+1;
    end
end
AlignIdentMatrix_RET_XCR_Sorted=sortrows(AlignIdentMatrix_RET_XCR_Sorted,-3);
% 90% cutoff RET_XCR
Row=1;AlignIdent_RET_XCR_Sorted_90=zeros(1,3);
for i=1:120
    if (AlignIdentMatrix_RET_XCR_Sorted(i,3)>90)
        AlignIdent_RET_XCR_Sorted_90(Row,:)=AlignIdentMatrix_RET_XCR_Sorted(i,1:3);
        Row=Row+1;
    end
end
TelomereConservedSize=cell(1); 
for i=1:16
    TelomereConservedSize{1}(i,1)=X_element_L_End(i)-X_element_L_Start(i);
    TelomereConservedSize{1}(i,3)=X_element_R_End(i)-X_element_R_Start(i);
    
    TelomereConservedSize{1}(i,5)=X_element_Comb_Repeat_L_End(i)-X_element_Comb_Repeat_L_Start(i);
    TelomereConservedSize{1}(i,7)=X_element_Comb_Repeat_R_End(i)-X_element_Comb_Repeat_R_Start(i);
    
    TelomereConservedSize{1}(i,9)=Y_Prime_element_L_End(i)-Y_Prime_element_L_Start(i);
    TelomereConservedSize{1}(i,11)=Y_Prime_element_R_End(i)-Y_Prime_element_R_Start(i);
end

TelomereConservedSize{1}(:,2)=TelomereConservedSize{1}(:,1)/max(TelomereConservedSize{1}(:,1));
TelomereConservedSize{1}(:,4)=TelomereConservedSize{1}(:,3)/max(TelomereConservedSize{1}(:,3));
TelomereConservedSize{1}(:,6)=TelomereConservedSize{1}(:,5)/max(TelomereConservedSize{1}(:,5));
TelomereConservedSize{1}(:,8)=TelomereConservedSize{1}(:,7)/max(TelomereConservedSize{1}(:,7));
TelomereConservedSize{1}(:,10)=TelomereConservedSize{1}(:,9)/max(TelomereConservedSize{1}(:,9));
TelomereConservedSize{1}(:,12)=TelomereConservedSize{1}(:,11)/max(TelomereConservedSize{1}(:,11));





