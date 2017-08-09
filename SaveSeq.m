function SaveSeq()
    % Organism: Saccharomyces cerevisiae [S288c strain]; Chromosomes: 16 chromosomes; Data Source: GenBank [NCBI].
    % GenBank Accesion Number: {'BK006935',....,'BK006949'};['BK006934' is for Chr VIII (except this, all in proper sequence)]
   
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
        disp(['Saving Chr:',num2str(ChrNum)]);
        disp('Chr Seq');getgenbank(AccessionNum{ChrNum},'SequenceOnly',true,'ToFile',[Path,'CHR\Chr_',num2str(ChrNum),'.fsa']);
        disp('LET Seq');getgenbank(AccessionNum{ChrNum},'SequenceOnly',true,'PARTIALSEQ',[1,LETEnd(ChrNum)],'ToFile',[Path,'LET\Chr_',num2str(ChrNum),'.fsa']);
        disp('RET Seq');getgenbank(AccessionNum{ChrNum},'SequenceOnly',true,'PARTIALSEQ',[RETStart(ChrNum),RETEnd(ChrNum)],'ToFile',[Path,'RET\Chr_',num2str(ChrNum),'.fsa']);
        disp('CEN Seq');getgenbank(AccessionNum{ChrNum},'SequenceOnly',true,'PARTIALSEQ',[CENStart(ChrNum),CENEnd(ChrNum)],'ToFile',[Path,'CEN\Chr_',num2str(ChrNum),'.fsa']);
    end
    clc;disp('ALL SEQUENCE DATA SAVED');
end