function [tbScalar, tbVec ]=getPerfMetricROCcompute(predictList,trueList,poslabel, plotOption)
         if size(predictList,2)>1
            predictList=reshape(predictList,numel(predictList),1);
            trueList=reshape(trueList,numel(trueList),1);
         end
    Npos = nnz(trueList==poslabel);
    Nneg = nnz(trueList~=poslabel);
    len = size( predictList, 1 ); 
    [~, ord ] = sort(predictList , 1,   'descend' );    
    rank_t = zeros( size(predictList) );          
    rank_t(ord ) = [1: len ]'; 
    predictList =1-rank_t/len;   
         
         low=min(predictList);
         high=max(predictList);
         threshold=linspace(high,low,1000);
         
         Sen=zeros(1,1000); Spe=zeros(1,1000); Pre=zeros(1,1000);
         Acc=zeros(1,1000);
         ind_cut = inf; 
         for I=1:1000
             Vector=zeros(numel(predictList),1);
             v=predictList>=threshold(I);
             Vector(v)=1;
             
             tp=sum(Vector==1&trueList==poslabel);
             tn=sum(Vector==0&trueList~=poslabel);
             np=sum(Vector==1&trueList~=poslabel);
             nn=sum(Vector==0&trueList==poslabel);
             
             if tp+nn==0
                Sen(I)=1;
             else
                Sen(I)=tp/(tp+nn);
             end
% %            [Npos,  (tp+nn) ] 
             if tn+np==0
                Spe(I)=1;
             else
                Spe(I)=tn/(tn+np);
             end
          
             if tp+np==0
                Pre(I)=1;
             else
                Pre(I)=tp/(tp+np);
             end
            
             Acc(I)=(tn+tp)/(tn+tp+np+nn);
             
             if ind_cut>abs( Npos-(tp+np) ); ind_cut = I; end 
         end
         
% %          Sen=[0,Sen];Spe=[1,Spe];Pre=[1,Pre];
         Sen=[0,Sen];Spe=[1,Spe];Pre=[Pre(1),Pre];
         Acc=[sum(trueList~=poslabel)/length(trueList) Acc];
         ind_cut = ind_cut + 1;
         
         AUC=abs(trapz(1-Spe,Sen)); 
         
         tbScalar = table;  
         tbScalar.AUROC   = AUC;
         %
         tbScalar.Acc     =  (Acc(ind_cut));
         tbScalar.Sen     =  (Sen(ind_cut));
         tbScalar.Spe     =  (Spe(ind_cut));
         tbScalar.Pre     =  (Pre(ind_cut));
         %
         tbScalar.mAcc     = mean(Acc(:));
         tbScalar.mSen     = mean(Sen(:));
         tbScalar.mSpe     = mean(Spe(:));
         tbScalar.mPre     = mean(Pre(:));
         
     
        tbVec=table;
        tbVec{'Acc',1} = Acc;
        tbVec{'Pre',1} = Pre;
        tbVec{'TPR_Rec_Sen',1} = Sen;   
        tbVec{'FPR',1} = 1-Spe;
        tbVec{'Spe',1} = Spe; 
         

         if plotOption==1
            plot(1-Spe,Sen);
            axis([-0.01 1.00 0 1.01]);
            xlabel('1-Spe');
            ylabel('Sen');
            
            figure;
            plot(Sen,Pre);
            axis([0 1.01 0 1.01]);
            xlabel('Sen');
            ylabel('Pre');
         end
         
         
