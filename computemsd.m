function [data] = computemsd(tRange,xy,sqrtmsd,sqrtt)

% compute mean-squared-displacement given a trajectory and timestep range 
% INPUTS
% tRange : range of time values 
% xy : x,y position 
% sqrtmsd : 
% sqrtt :
% OUTPUTS
% data: struct with xy,msd,meanval,stddev,tRange
% 
% 

if sqrtt
    tau = sqrt(tRange);
else
    tau = tRange;
end
    
for tau = 1:20
    % MSD
    msd{tau} = sum((xy(end:-1:1+tau,:) - xy(end-tau:-1:1,:)).^2,2);
    
    meanval(tau) = mean(msd{tau});
    stddev(tau) = std(msd{tau});
    
end

if sqrtmsd
        % Root MSD
        msd = sqrt(msd);
end

data.xy = xy;
data.msd = msd;
data.meanval = meanval;
data.stddev = stddev;
data.tRange = tRange;


