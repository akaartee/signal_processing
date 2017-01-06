function [ rrCategory ] = classifyRrIntervals( rrIntervalsInMs )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

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



end

