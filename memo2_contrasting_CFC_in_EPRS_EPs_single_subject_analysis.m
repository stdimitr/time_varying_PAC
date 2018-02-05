%%% Figure4  from ''A novel CFC-based biomarker of amnestic MCI''
%
% 1) loading the ERP and EP single-trials of a non-impaired subject
% 2) computing the TV_PAC-estimated among the 7 frequency bands for different temporal segmends
% 3) deriving ''summarizing''  profiles and comparing them between ERP and EP responses 
%
load single_trial_data EP_STdata ERP_STdata
Fs=1024; % sampling frequency
t=-100:923; % the time axis
FR_Bands=[1 4; 4 7.5; 8 10; 10 13; 13 20; 20 30; 30 45]
FR_names=['delta '; 'theta ';'alpha1';'alpha2';'beta1 ';'beta2 ';'gamma1'];

%% TV_PAC_patterns derived across-trials  [#fr-bands x #fr-bands x #temporal segments]  
Nsegments=250;
%ERP
STdata=ERP_STdata;Ntrials=size(STdata,1); 
ERP_PAC_matrix=[]; for i1=1:size(FR_Bands,1)-1; for i2=i1+1:size(FR_Bands,1);
   Pf1=FR_Bands(i1,1);  Pf2=FR_Bands(i1,2);Af1=FR_Bands(i2,1);  Af2=FR_Bands(i2,2);
      [tPAC,times]=moving_multitrial_pac(STdata,Fs,Pf1,Pf2,Af1,Af2,Nsegments);
      ERP_PAC_matrix(i1,i2,:)=tPAC; ERP_PAC_matrix(i2,i1,:)=tPAC; % technically symmetrize
    end,end,

%%% EPs
STdata=EP_STdata;Ntrials=size(STdata,1); 
EP_PAC_matrix=[]; for i1=1:size(FR_Bands,1)-1; for i2=i1+1:size(FR_Bands,1);
   Pf1=FR_Bands(i1,1);  Pf2=FR_Bands(i1,2);Af1=FR_Bands(i2,1);  Af2=FR_Bands(i2,2);
      %plv=multitrial_pac(STdata,1024,Pf1,Pf2,Af1,Af2,t1,t2);
      [tPAC,times]=moving_multitrial_pac(STdata,1024,Pf1,Pf2,Af1,Af2,Nsegments);
      EP_PAC_matrix(i1,i2,:)=tPAC; EP_PAC_matrix(i2,i1,:)=tPAC; % technically symmetrize
    end,end,


%% computing the ''summarizing'' profiles
for i=1:Nsegments, FF_ERPs(i)=mean(mean(ERP_PAC_matrix(:,:,i))); max_FF_ERPs(i)=max(max(ERP_PAC_matrix(:,:,i))),end
for i=1:Nsegments, FF_EPs(i)=mean(mean(EP_PAC_matrix(:,:,i))); max_FF_EPs(i)=max(max(EP_PAC_matrix(:,:,i))),end

figure(2),clf,subplot(3,1,1),plot(t,mean(ERP_STdata),'b',t,mean(EP_STdata),'k',t,zeros(1,1024),':k','linewidth',2), 
axis([-100 1000 -2000 2000]), legend('ERP','EP'),title('averaged-responses')
subplot(3,1,2), plot(t(times),FF_ERPs,'b',t(times),FF_EPs,'k','linewidth',2),ylabel('aver.PAC-value')%legend('ERP','EP')
axis([-100 1000 0 0.5]),title('TV-  PAC averaged- across fr.pairs - profiles')
subplot(3,1,3), plot(t(times),max_FF_ERPs,'b',t(times),max_FF_EPs,'k','linewidth',2),ylabel('max PAC-value')%legend('ERP','EP')
axis([-100 1000 0 0.5]),title('TV-  PAC maximum-across fr.pairs  - profiles')















load ensemble_averages_8fr_bands
subj_label=[zeros(1,15), ones(1,25)];t=-100:923; % NonImpaired (NI)label-->0 ,Mild Cognitive Impaired (MCI)label-->1
 
FR_Bands=[1 45;1 4; 4 8; 8 10;10 12; 13 20; 20 30; 30 45];
FR_names=['wide  ';'delta '; 'theta ';'alpha1';'alpha2';'beta1 ';'beta2 ';'gamma1'];

%____ For each group we derive and Grand-Average 
NI_ERP=mean(ALL_subject_ERPave(:,:,subj_label==0),3);MCI_ERP=mean(ALL_subject_ERPave(:,:,subj_label==1),3);

freq_SCORES=[];for i=1:size(FR_Bands,1);
 NI_AVEs=squeeze(ALL_subject_ERPave(i,:,subj_label==0))'; MCI_AVEs=squeeze(ALL_subject_ERPave(i,:,subj_label==1))';
[IDX,Z] = rankfeatures([NI_AVEs;MCI_AVEs]',[zeros(1,15),ones(1,25)],'criterion', 'wilcoxon');
 freq_SCORES(i,:)=Z; 
end 
YSCORES=stack_plot(freq_SCORES,4);
YNI_ERP=stack_plot(NI_ERP,520);YMCI_ERP=stack_plot(MCI_ERP,520);



figure(1),clf,subplot(1,3,1), plot(t,YNI_ERP,'linewidth',2),title('NI: ERP Grand-Average'),grid,legend(FR_names)

plot(t,mean(EP_STdata),t,mean(ERP_STdata))
subplot(3,1,1),plot(t,mean(ERP_STdata),'b',t,mean(EP_STdata),'k',t,zeros(1,1024),':k','linewidth',2), legend('ERP','EP')


subplot(1,3,2), plot(t,YMCI_ERP,'linewidth',2),title('MCI: ERP Grand-Average'),grid
subplot(1,3,3), plot(t,YSCORES),title('Wilcoxon Score'),grid, xlabel('#sample')

