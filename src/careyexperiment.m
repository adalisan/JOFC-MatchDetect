                                                   
                                                                     
                                             
function [ rhovalues,prcntcrrct ] = CareyExperiment(numexper,deltarho,r,s)
% [ rhovalues,prcntcrrct ] = CareyExperiment(numexper,deltarho,r,s)
%  does numexper experiments for rho values 0 to .5 in increments 
% of deltarho;  r vertex correspondences are given and s are estimated.
% ready Jan 9, 2012 *

rhovalues=[0:deltarho:1];
numrho=length(rhovalues);

Results=zeros(numexper,numrho);

for k=1:numexper
    A=round(rand(r+s,r+s));
    A=A-triu(A);
    A=A+A';
    for q=1:numrho 
        for i=1:r+s
            for j=1:i-1
                temp=rand;
                if A(i,j)==1
                    B(i,j)=temp<1-rhovalues(q);
                    B(j,i)=B(i,j);
                else
                    B(i,j)=temp<rhovalues(q);
                    B(j,i)=B(i,j);  
                end
            end
        end
        mixit=[ [1:r]  r+randperm(s)];
        A=A(mixit,mixit);
        [corr,P] = graphmatchHARDSEED( A,B,r );
        Results(k,q)=sum(mixit(r+1:r+s)'==corr(r+1:r+s))/s;
    end
end
prcntcrrct=(ones(1,numexper)*Results)/numexper;
plot(rhovalues,prcntcrrct,'*')
end
