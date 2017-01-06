
clearvars;
clc;

data = dlmread('207.txt',' ');
% 100 ok.  101 ok.  102 ok.  103 ok.
% 104 a few false detections due to noisy signal.
% 105 too many beats. smaller threshold divider helps.
% 106 ok.  107 quite ok.  108 extra beats.  109 ok.
% 111 a few extra beats.
% 112 OK.  113 OK.  114 ok.  115 OK.  116 pretty ok.  117 OK.  118 ok.  119 pretty ok.
% 121 OK.  122 OKKK.  123 OK.  124 ok.
% 200 pretty ok.  201 pretty ok.
% 202 PVC criteria results in missclassification.
% 203 PVC criteria results in missclassification.
% 205 ok.  207 ok, challenging signal.
% 208 some wide abnormal R-peaks are missed.
% 209 ok, correctly detected beats missclassified as ventricular flutter.
% 210 ok, but PVC criteria results in missclassification.
% 212 ok. 213 ok. 214 ok. 215 ok. 217 ok. 
% 219 ok, pretty ok. 220 OK. 221 ok, missing some abnormal beats.
% 222 ok, but PVC criteria results in missclassification.
% 223 ok.  228 ok.  230 ok.
% 231 ok, but heart block criteria isn't working.
% 232 ok, heart blocks found.
% 233 ok.  234 ok.

time = data(1:end,1);
ecg = data(:,2);

fs = 1/mean(diff(time))


% Due to poor resolution of the time samples it is best to calculate a new
% time vector using accurate sampling frequency.
timeNew = (0:length(ecg)-1)'./fs;

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

figure(1)
plot(timeNew,ecgPreFiltered,'k')
hold on
plot(timeNew,highpassed9HzEcg,'b')
plot(timeNew,highpassed3HzEcg,'m')
plot(timeNew,lowpassed5HzEcg,'g')
hold off
grid on


%% R-peak detection
threshold = mean(abs(highpassed9HzEcg))*4;
candidatePeaks(1) = threshold;

rIndeces = [];
rAmplitudes = [];

kk = 1;

prevPeak = 1;
prevPeakAmp = 0;

adaptiveThreshold = zeros(size(highpassed9HzEcg));

for i = 14:length(highpassed9HzEcg)-1
    % Detect a peak in highpassed9HzEcg.
    if highpassed9HzEcg(i) > highpassed9HzEcg(i+1) && highpassed9HzEcg(i) > highpassed9HzEcg(i-1) && highpassed9HzEcg(i) > 0
        if highpassed9HzEcg(i) > 0.05 && length(rIndeces) > 1
            candidatePeaks(end+1) = highpassed9HzEcg(i); %#ok<*SAGROW>
        end
        % Check amplitude threshold criteria.
        if highpassed9HzEcg(i) > threshold || (highpassed3HzEcg(i) > 1.2*threshold && highpassed9HzEcg(i) > threshold*0.7)
            % Check that the slopes of the peak are steep enough.
            if highpassed9HzEcg(i) - highpassed9HzEcg(i-8) > highpassed9HzEcg(i)*0.26 ...
                    || highpassed9HzEcg(i) - highpassed9HzEcg(i-13) > highpassed9HzEcg(i)*0.38
                if (length(rIndeces) > 1)
                    % Ignore a peak when its amplitude is lower than the
                    % previous one which occured very recently (< 200ms)
                    if ((i - rIndeces(end))/fs < 0.2)
                        if rAmplitudes(end) < highpassed9HzEcg(i)
                            kk = kk - 1;
                        else
                            adaptiveThreshold(i) = threshold;
                            continue
                        end
                    % Ignore T-wave peaks.
                    elseif (rAmplitudes(end) > highpassed9HzEcg(i)*1.5 && lowpassed5HzEcg(i) > highpassed9HzEcg(i) && rAmplitudes(end) > threshold*0.5)
                        if ecgPreFiltered(rIndeces(kk-1)) < rAmplitudes(kk-1)*2
                            adaptiveThreshold(i) = threshold;
                            continue
                        end
                    end
                end
                for iii = 1:10
                    candidatePeaks(end+1) = highpassed9HzEcg(i);
                end
                if length(candidatePeaks) > 18
                    threshold = mean(candidatePeaks(end-18:end))/2.1;
                else
                    threshold = mean(candidatePeaks(end-length(candidatePeaks)+1:end))/2.1;
                end
                rIndeces(kk) = i;
                rAmplitudes(kk) = highpassed9HzEcg(i);
                kk = kk + 1;
                adaptiveThreshold(i) = threshold;
                continue
            end
        end
        if length(candidatePeaks) > 18
            threshold = mean(candidatePeaks(end-18:end))/2.1;
        else
            threshold = mean(candidatePeaks(end-length(candidatePeaks)+1:end))/2.1;
        end
    end
    adaptiveThreshold(i) = threshold;
    % Detection for W-shaped low amplitude (abnormal) R-peaks.
    % Detect a peak in highpassed3HzEcg.
    if (highpassed3HzEcg(i) > highpassed3HzEcg(i+1) && highpassed3HzEcg(i) > highpassed3HzEcg(i-1) && highpassed3HzEcg(i) > 0)
        if (highpassed3HzEcg(i) > threshold*0.25) &&  ((i - prevPeak)/fs < 0.1) && (prevPeakAmp > highpassed3HzEcg(i)*0.7)
            notch1 = 0;
            notch2 = 0;
            if i <= length(highpassed3HzEcg)-11
                % Check that both notches of W-shape are found.
                for kkk = 3:10
                    if (highpassed3HzEcg(i-kkk) < -highpassed3HzEcg(i)*1.3 && highpassed3HzEcg(i+kkk) < -highpassed3HzEcg(i)*1.3)
                        if highpassed3HzEcg(i-kkk) < highpassed3HzEcg(i-kkk-1) && highpassed3HzEcg(i-kkk) < highpassed3HzEcg(i-kkk+1)
                            notch1 = 1;
                        end
                        if highpassed3HzEcg(i+kkk) < highpassed3HzEcg(i+kkk-1) && highpassed3HzEcg(i+kkk) < highpassed3HzEcg(i+kkk+1)
                            notch2 = 1;
                        end
                    end
                end
            end
            if notch1 + notch2 == 2
                rIndeces(kk) = i;
                rAmplitudes(kk) = highpassed3HzEcg(i);
                kk = kk + 1;
            end
        end
        prevPeak = i;
        prevPeakAmp = highpassed3HzEcg(i);
    % Detection for inverted R-peak (deep notch without a prominent peak).
    elseif (highpassed3HzEcg(i) < highpassed3HzEcg(i+1) && highpassed3HzEcg(i) < highpassed3HzEcg(i-1) ...
            && highpassed3HzEcg(i) < 0 && not(isempty(rIndeces)) )
        if ((-0.7*highpassed3HzEcg(i) > prevPeakAmp) && (-highpassed3HzEcg(i) > threshold*1.1) && (i-rIndeces(end))/fs > 0.2)
            rIndeces(kk) = i;
            rAmplitudes(kk) = highpassed3HzEcg(i);
            kk = kk + 1;
        end
    end
end

% Number of detected R-peaks
Number_of_detected_R_peaks = kk-1

figure(1)
hold on
plot(timeNew(rIndeces),rAmplitudes,'r*','linewidth',2)
plot(timeNew,adaptiveThreshold,'k')
hold off
grid on
leg = legend('ecgPreFiltered','highpassed9HzEcg','highpassed3HzEcg','lowpassed5HzEcg','R-peaks','adaptive threshold');
set(leg,'fontsize',15)


%% RR-interval calculation

rPeakTimeStamps = timeNew(rIndeces);
rrIntervalsInMs = diff(rPeakTimeStamps).*1000;

% Correction for short double beat detections
for i = 2:length(rrIntervalsInMs)-1
    if rrIntervalsInMs(i) < 100 && rrIntervalsInMs(i-1) > 150 && rrIntervalsInMs(i+1) > 150
        rrIntervalsInMs(i) = rrIntervalsInMs(i-1) + rrIntervalsInMs(i);
        rrIntervalsInMs(i-1) = 0;
    end
end
rPeakTimeStamps(rrIntervalsInMs == 0) = [];
rrIntervalsInMs(rrIntervalsInMs == 0) = [];


% Additional correction for remaining short duration intervals
for ind = 1:8
    for i = 2:length(rrIntervalsInMs)-1
        if rrIntervalsInMs(i) < 200
            rrIntervalsInMs(i) = rrIntervalsInMs(i-1) + rrIntervalsInMs(i);
            rrIntervalsInMs(i-1) = 0;
        end
    end
    rPeakTimeStamps(rrIntervalsInMs == 0) = [];
    rrIntervalsInMs(rrIntervalsInMs == 0) = [];
end

Number_of_approved_R_peaks = length(rrIntervalsInMs)+1;

%% Classification of RR-intervals

prevPVC = -2;
rrCategory = ones(size(rrIntervalsInMs));
skip = 0;

for i = 2:length(rrIntervalsInMs)-1
    if skip < 1
        rrCategory(i) = 1;

        rr1 = rrIntervalsInMs(i-1);
        rr2 = rrIntervalsInMs(i);
        rr3 = rrIntervalsInMs(i+1);

        % Rule 1, ventricular flutter beats
        if rr2 < 600 && 1.8*rr2 < rr1
            rrCategory(i) = 3;
            for k = i+1:length(rrIntervalsInMs)-1
                rr1k = rrIntervalsInMs(k-1);
                rr2k = rrIntervalsInMs(k);
                rr3k = rrIntervalsInMs(k+1);
                if (rr1k < 700 && rr2k < 700 && rr3k < 700) || (rr1k + rr2k + rr3k < 1700)
                    rrCategory(k) = 3;
                else
                    if k <= i+3
                       rrCategory(i:k) = 1;
                    else
                        skip = k-i-1;
                    end
                    break
                end

            end
        end
        
        % Rule 2, premature ventricular contractions
        if ((1.15*rr2 < rr1) && (1.15*rr2 < rr3)) || ...
                ((abs(rr1-rr2) < 300) && ((rr1 < 800) && ( rr2 < 800)) ...
                && (rr3 > 1.2*mean([rr1,rr2])) || ...
               ((abs(rr2-rr3) < 300) && ((rr2 < 800) && ( rr3 < 800)) ...
                && (rr1 > 1.2*mean([rr2,rr3])) ) )
            rrCategory(i) = 2;
            prevPVC = i;
        end
        
        % Rule 3, second degree heart block
        if (2200 < rr2 && rr2 < 3000) && ((abs(rr1-rr2) < 200) || (abs(rr2-rr3) < 200))
            rrCategory(i) = 4;
        end
    else
        skip = skip - 1;
    end
end

number_of_normal_beats = size(rrCategory(rrCategory == 1), 1);
number_of_premature_ventricular_contractions = size(rrCategory(rrCategory == 2), 1);
number_of_ventricular_flutter_beats = size(rrCategory(rrCategory == 3), 1);
number_of_second_degree_heart_block_beats = size(rrCategory(rrCategory == 4), 1);

sprintf('Number of approved R-peaks: %d\nNumber of normal beats: %d\nNumber of premature ventricular contractions: %d\nNumber of ventricular flutter beats: %d\nNumber of 2nd degree heart block beats: %d',...
    Number_of_approved_R_peaks,number_of_normal_beats,number_of_premature_ventricular_contractions,number_of_ventricular_flutter_beats,number_of_second_degree_heart_block_beats)


%% Additional features
hrInBpm = 60./(rrIntervalsInMs./1000);

% Heart rate variability (HRV) features
RMSSD = sqrt( mean( (rrIntervalsInMs(2:end)-rrIntervalsInMs(1:end-1)).^2 ) ) %#ok<*NOPTS>
SDRR = std(rrIntervalsInMs)
SDSD = std( rrIntervalsInMs(1:end-1)-rrIntervalsInMs(2:end) )
SD1 = (1/sqrt(2))*SDSD
SD2 = sqrt( 2*SDRR^2 - 0.5*SDSD^2 )


figure(2)
subplot(3,1,1)
plot(rPeakTimeStamps(2:end),rrIntervalsInMs,'r')
grid on
leg = legend('RR-Intervals');
set(leg,'fontsize',16)
ylabel('RRI [ms]','fontsize',16)

subplot(3,1,2)
plot(rPeakTimeStamps(2:end),hrInBpm ,'b')
grid on
leg = legend('HR');
set(leg,'fontsize',16)
ylabel('HR [bpm]','fontsize',16)

subplot(3,1,3)
plot(rPeakTimeStamps(2:end), rrCategory)
axis([0,length(rrCategory),0,5])
grid on
leg = legend('Arrhythmia category');
set(leg,'fontsize',16)
ylabel('Category','fontsize',16)

figure(3)
for pp=1:length(rrIntervalsInMs)-1
    if rrCategory(pp) == 1
        plot(rrIntervalsInMs(pp),rrIntervalsInMs(pp+1),'k.')
    elseif rrCategory(pp) == 2
        plot(rrIntervalsInMs(pp),rrIntervalsInMs(pp+1),'r.')
    elseif rrCategory(pp) == 3
        plot(rrIntervalsInMs(pp),rrIntervalsInMs(pp+1),'m.')
    elseif rrCategory(pp) == 4
        plot(rrIntervalsInMs(pp),rrIntervalsInMs(pp+1),'b.')
    end
    if pp == 1
        hold on
    end
end
hold off
grid on
title(sprintf('Poincare plot. RMSSD: %.1f, SDRR: %.1f, SDSD: %.1f, SD1: %.1f, SD2: %.1f.',RMSSD,SDRR,SDSD,SD1,SD2))






