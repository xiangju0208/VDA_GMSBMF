function [CVdata] = getKfoldCrossValidMatIndSet(Wdv, nfold,CVtype,NegType, seed )
% % Wdv = rand(10,5)>0.5; nfold  = 2; CVtype = 'CVc'; seed = [] ;   
% % if ~exist('CVtype','var') || isempty( CVtype )
% %     CVtype = 'CVa'; 
% % %     CVtype = 'CVr'
% % %     CVtype = 'CVc'
% % end
% % if exist('NegType','var') && ~isempty( NegType )
% %     NegType = 'Unlabel'  ; 
% % end
if exist('seed','var') && ~isempty( seed )
    rng('default')
    rng(seed) 
end

IndSet_pos_test = repmat({[]},nfold,1  ); 
IndSet_neg_test = repmat({[]},nfold,1  ); 
IndSet_pos_train = repmat({[]},nfold,1  ); 
IndSet_neg_train = repmat({[]},nfold,1  ); 
Wdv_mark = zeros( size(Wdv) ); 
Wdv_mark(Wdv~=0 )= 1; 
Wdv_mark(Wdv==0 )= -1;     
switch CVtype
    case  'CVa'  % split all positive and negative pairs 
        Ind_pos = find(Wdv_mark==1);  
        Ind_neg = find(Wdv_mark==-1);  
        Ipos    = crossvalind('Kfold',length(Ind_pos),nfold );     %%%% ��Ҫ�޸� 
        Ineg    = crossvalind('Kfold',length(Ind_neg),nfold );
        for i_fold = 1:nfold
            IndSet_pos_test{i_fold} = Ind_pos(Ipos==i_fold);  
            IndSet_pos_train{i_fold} = Ind_pos(Ipos~=i_fold);  
            if  strcmpi( NegType, 'Unlabel'  )
                IndSet_neg_test{i_fold} = Ind_neg ; 
                % % % ȷ�����Լ����Ժ�����������ƽ�Դ�,��˽����Լ�����������Ҳ��������ѵ����,
                IndSet_neg_train{i_fold} = union(IndSet_neg_test{i_fold}, IndSet_pos_test{i_fold} ) ;                
            else
                IndSet_neg_test{i_fold} = Ind_neg(Ineg==i_fold);   
                IndSet_neg_train{i_fold} = Ind_neg(Ineg~=i_fold);              
            end
            %
        end 
    

    otherwise; error(['There is on cross-validation definition of : ' , CVtype  ] ); 
end 
CVdata = []; 
CVdata.MatIndSet_pos_test = IndSet_pos_test; 
CVdata.MatIndSet_neg_test = IndSet_neg_test; 
CVdata.MatIndSet_pos_train =IndSet_pos_train; 
CVdata.MatIndSet_neg_train = IndSet_neg_train; 

end
