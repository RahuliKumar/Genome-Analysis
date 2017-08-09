DataPath='D:\IITD\All\Academics\SEM 9\BED851 (Project)\MATLAB Data';
SortedInfoMatrix=load([DataPath,'\SortedInfoMatrix_LET_LET']);
SortedInfoMatrix=SortedInfoMatrix.SortedInfoMatrix;

Sotrted_New=cell(1);

for i=1:15
    k=zeros(16,3);
    for j=1:16
        k(j,1:3)=SortedInfoMatrix{1,1}(i,1+(j-1)*4:3+(j-1)*4);
     end
    k=sortrows(k,-3);
    Sotrted_New{1}(:,1+(i-1)*3:3+(i-1)*3)=k;
end
