function [ rrIntervalsInMs, rPeakTimeStamps ] = calculateRrIntervals( timeNew, rIndeces )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% RR-interval calculation

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


end

