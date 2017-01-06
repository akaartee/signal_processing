function [ ecgPreFiltered, highpassed9HzEcg, highpassed3HzEcg, lowpassed5HzEcg ] = prefilter( ecg, fs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Calculate filter coefficients (fs/2 = Nyquist frequency. 1Hz/(Nyquist)=Normalized freq.)
[bBs,aBs] = butter(2, [ 57/(fs/2),63/(fs/2)],'stop');
[bBp,aBp] = butter(1, [ 1/(fs/2),37/(fs/2)],'bandpass');
[b2,a2]   = butter(1, 9.0/(fs/2),'high');
% [b2,a2]   = butter(1, [9.0/(fs/2),40/(fs/2)],'bandpass');
[b3,a3]   = butter(2, 3.5/(fs/2),'high');
% [b3,a3]   = butter(2, [3.5/(fs/2),40/(fs/2)],'bandpass');
[bLp,aLp] = butter(2, [1.5/(fs/2),5/(fs/2)],'bandpass');

% Remove 60Hz power line noise.
ecgPreFiltered = filtfilt(bBs,aBs,ecg);

% Bandpass ECG to remove high frequency noise and some baseline wondering.
ecgPreFiltered = filtfilt(bBp,aBp,ecgPreFiltered);

% High passed signal with 9Hz cutoff frequency for better damping of low frequency artefacts.
highpassed9HzEcg = filtfilt(b2,a2,ecgPreFiltered);

% Alternative high passed signal with 3.5Hz cutoff frequency for better detection of abnormal beats.
highpassed3HzEcg = filtfilt(b3,a3,ecgPreFiltered);

% Low passed signal with 5Hz cutoff for identifying T-waves and false R-peak detections.
lowpassed5HzEcg = filtfilt(bLp,aLp,ecgPreFiltered);


end

