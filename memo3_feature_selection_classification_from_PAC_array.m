%%% Figure7 and Table 3 from ''A novel CFC-based biomarker of amnestic MCI''
%
% 1) loading the array of TV_PAC-measurements for ALL subjects
% 2) bootstrap-procedure for selecting features
% 3) training-testing classifier
% 4) representing the selected features, i.e. the PAC-biomarker

load  FV_from_TV_PAC_analysis  % ERP_PAC_array [#subjects x #bands x #bands x #segments]
t=-100:923;times=10:20:990;% the latencies that correspond to the 50 segments we used to sample the PAC-estimates from
Nsegments=size(ERP_PAC_array,4);
 

% PART 2)
%%%%% BootStrap : Random Sampling with Replacement
                                        % ammendments todo: stratified sampling
 FVs=reshape(ERP_PAC_array,40,7*7*Nsegments);
%_________ perform Feture-Ranking over multiple samples and averaging the Wilcoxon_Score 
%%% if you are patient run the next loop for 1000 times  
repeatZ=[]; for ri=1:100; ri
   rlist=(ceil(40*rand(1,40)));X_Train=FVs(rlist,:); Train_labels= subj_label(rlist);
   [IDX, Z] = rankfeatures(X_Train',Train_labels,'criterion', 'wilcoxon');
   repeatZ(ri,:)=Z;   
 end   
   Feat_Score=mean(repeatZ)./std(repeatZ) ; % estimating a score for each feature     
   [ranked_Score,slist]=sort(Feat_Score); % ranking the features
    % selecting the 30=60/2 most discriminative features
  sel=slist(end-59:2:end);FVs=FVs(:,sel); Str=ranked_Score(end-59:2:end);% the strength of the Score  

 % these are computed for defined the sign of the difference (used below in Fig.7)
 NI_CFC_level=median(FVs(subj_label==0,:)); MCI_CFC_level=median(FVs(subj_label==1,:)); 

 % we may visually check the potential of the selected features 
 % via  the follwing scatter plot derived via MDS
 M=cmdscale(pdist(FVs(:,:))) ;figure(10),clf, plot(M(1:15,1),M(1:15,2),'bo',M(16:40,1),M(16:40,2),'ro' )
 for i=1:40, text(M(i,1),M(i,2),num2str(i)),end   % note: #10 / #36  can be considered as typical subjects for NI/MCI group
 
 
 % PART 3)
 %%%%%%%%%%%%%%%%%%%%%%%%555
 %___SVM classifier ________  
 %%%%  multiple repetitions of 2-fold cross-validation scheme

 for i=1:200;
rlist=randperm(40); 
Xtrain=FVs(rlist(1:35),:); train_labels=subj_label(rlist(1:35));
Xtest=FVs(rlist(36:40),:); test_labels=subj_label(rlist(36:40));
svmStruct = svmtrain(Xtrain,train_labels,'kernel_function','linear');
Pred_labels = svmclassify(svmStruct,Xtest);
Error(i)=sum(abs(Pred_labels'-test_labels))/numel(test_labels);
end
mean(Error) % ~6.40


%%%%  or Leave-one-out subject -validation for the SVM Classifier
Error=[]; for i=1:40   %subj_label
      list=setdiff([1:40],i);
      Xtest=FVs(i,:); test_label=subj_label(i);
      Xtrain=FVs(list,:); train_labels= subj_label(list);      
svmStruct = svmtrain(Xtrain,train_labels,'kernel_function','linear');
Pred_label = svmclassify(svmStruct,Xtest);
Error(i)=Pred_label-test_label ;
end
mean(abs(Error)) %~2.5

% PART 4)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visualizing the most significant PAC-features 

REL_CHANGE_sel=(MCI_CFC_level-NI_CFC_level)./NI_CFC_level;
Q2=zeros(1,7*7*50); Q2(sel)=REL_CHANGE_sel;
TV_REL_CHANGE=reshape(Q2,7,7,50);
TV_Score=reshape(Feat_Score,7,7,50); % from long vector bring to TV-format
FF=[];for i=1:50;FF(i)=max(max(TV_Score(:,:,i))); end
[pks,locs] = findpeaks(FF,'minpeakheight',2.5); % used in Fig.7A

%%%% Paper-Figure-7_upper_panel
figure(3),clf,subplot(2,2,1),plot(t(times),FF,t(times(locs)),pks,'ro'), ylabel('Wscore*') ,xlabel('time')
Amax=max(TV_Score,[],3);lim1=min(Amax(:));lim2=max(Amax(:));%Asum=max(TV_Score,3);lim1=min(Asum(:));lim2=max(Asum(:));
[sources,destinations,Strengths]=find(triu(Amax));
colors=summer(64);N=7; z=exp(j*2*pi*[0:1/N:1-1/N]);   z=[z(5); z(6); z(7);z(1); z(2); z(3); z(4);z(5)]';
subplot(2,2,2), plot(z,'k.'),hold
for i=1:length(sources);
color_index=max(1,(round(((Strengths(i)-lim1)/(lim2-lim1))*64)));
quiver(real(z(sources(i))),imag(z(sources(i))),real(z(destinations(i))-z(sources(i))),imag(z(destinations(i))-z(sources(i))),...
'AutoScaleFactor',0.9, 'linewidth',1,'color',colors(color_index,:))
end,
 plot(z,'o','markersize',25,'markerfacecolor',[0.6 0.6 0.6]),hold, axis ('off')
FR_names=['ä '; 'è ';'á1';'á2';'â1';'â2';'ã1'];
 for ii=1:7; text(real(z(ii)),imag(z(ii)),FR_names(ii,:),'color','r'),end
title('CFC-relational graph based on maximal PAC')
 
 %____ Paper-Figure-7-lower  panel
figure(4),clf,
colors=jet(64);
limit1=min(REL_CHANGE_sel);limit2=max(REL_CHANGE_sel);       %limit1=min(Feat_Score); limit2=max(Feat_Score);
for i=1:numel(locs);
A=TV_REL_CHANGE(:,:,locs(i))';[sources,destinations,Strengths]=find(triu(A));
subplot(2,ceil(numel(locs)/2),i),plot(z,'o','markersize',5,'markerfacecolor',[0.7 0.7 0.7],'markeredgecolor',[0 0 0]),hold
for i2=1:length(sources);
color_index=max(1,(round(((Strengths(i2)-limit1)/(limit2-limit1))*64)));
quiver(real(z(sources(i2))),imag(z(sources(i2))),real(z(destinations(i2))-z(sources(i2))),imag(z(destinations(i2))-z(sources(i2))),...
    'AutoScaleFactor',0.95, 'linewidth',2,'color',colors(color_index,:))
end,
 hold, axis ('off')
 title(num2str( t(times(locs(i)))))         
end


