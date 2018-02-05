function [temporal_plv,Times]=moving_multitrial_pac(signals,srate,Pf1,Pf2,Af1,Af2,Nsegments)
%estimation of phase-amplitude cross-frequency coupling
%[temporal_plv,Times]=multitrial_pac(signals,srate,Pf1,Pf2,Af1,Af2,Nsegments)
% Inputs:
% signals = Single_trial- time series
% srate = frequency range
% Pf1,Pf2 = frequency range of the low-frequency
% Af1,Af2 = frequency range of the high-frequency
%    Nsegments : # temporal segments ....
%               
% Outputs:
% temporal_plv = Phase-Amplitude Coupling array over the different temporal segments
%Times -> the middle-latency of each segment


%adapted from Cohen MX."Assessing transient cross-frequency coupling in EEG data"
% Journal of Neuroscience Methods 168 (2008) 494–499

%W.D. Penny et al., "Testing for nested oscillation"
%Journal of Neuroscience Methods 174 (2008) 50–61

% DIMITRIADIS STAVROS 9/2017  
%     http://neuroinformatics.gr/ 
%     http://users.auth.gr/~stdimitr/index.html



%get the phase of the low-frequency band
%PhaseFreq=eegfilt(signal,srate,Pf1,Pf2); % this is just filtering 
[bb_p,aa_p]=butter(3,[Pf1/(srate/2),Pf2/(srate/2)]);
low_filtered_signals=filtfilt(bb_p,aa_p,signals')';
Phase_signals=angle(hilbert(low_filtered_signals'))'; % this is getting the phase time series
Ntrials=size(signals,1);


%STEP 1: Filtering the original signal in the range of high-frequency range
[bb,aa]=butter(3,[Af1/(srate/2),Af2/(srate/2)]);
high_filtered_signals=filtfilt(bb,aa,signals')';

%STEP 2:Get the analytic signal 
Env_high_filtered_signals=abs(hilbert(high_filtered_signals'))'; % getting the amplitude envelope

%STEP 3: Filtering the amplitude of analytic signal of the high-frequency range within the
%    frequency range of the low-frequency band

lowfromhigh=filtfilt(bb_p,aa_p,Env_high_filtered_signals')'; 
low_Env_high_filtered_signals=lowfromhigh-repmat(mean(lowfromhigh),Ntrials,1);

%STEP 4:Get the phase
Amp_phase_signals=angle(hilbert(low_Env_high_filtered_signals'))';

Ntime=size(signals,2);
Nstep=floor(Ntime/Nsegments);
for ii=1:Nsegments  % repear over different segments
    tt=[1+(ii-1)*Nstep:ii*Nstep]; Times(ii)=tt(round(Nstep/2)); % the center-latency of each window
plv=abs(sum(sum(exp(1i*(Phase_signals(:,tt)-Amp_phase_signals(:,tt))))))/(Ntrials*(Nstep));
temporal_plv(ii)=plv; end

end