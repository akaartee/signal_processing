function [ rIndeces, rAmplitudes, adaptiveThreshold ] = detectRPeaks( timeNew, ecgPreFiltered,highpassed9HzEcg, highpassed3HzEcg, lowpassed5HzEcg, fs )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% R-peak detection
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



end

