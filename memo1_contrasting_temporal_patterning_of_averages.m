%%% Figure3 from ''A novel CFC-based biomarker of amnestic MCI''
%
% loading the ensemble_averages from all subjects (and frequency bands)
%       and contrasting them by means of wilcoxon-score
%

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
subplot(1,3,2), plot(t,YMCI_ERP,'linewidth',2),title('MCI: ERP Grand-Average'),grid
subplot(1,3,3), plot(t,YSCORES),title('Wilcoxon Score'),grid, xlabel('#sample')

