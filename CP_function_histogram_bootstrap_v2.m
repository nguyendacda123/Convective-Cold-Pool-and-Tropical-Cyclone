function [freq, CI] = CP_function_histogram_bootstrap_v2(dat1,dat2,bins,Msize)
    % dat1 - Cold pool time (minute)
    % dat2 - instrument time (days)
    % compute median and confidence interval of median using bootstrapping
    % method
    
    rng(1);
    nbin = length(bins)-1;
    N1 = length(dat1);
    N2 = length(dat2);
    CI = nan(nbin,2);
    count = nan(nbin,Msize);
    parfor k = 1:Msize        
        dats1 = datasample(dat1,N1);
        dats2 = datasample(dat2,N2);
        for i = 1:nbin
            count1 = nansum(dats1 >= bins(i) & dats1 < bins(i+1))*10; % minute
            count2 = nansum(dats2 >= bins(i) & dats2 < bins(i+1))*600/86400; % day
            count(i,k) = count1/count2;
        end         
    end
    
    freq = nanmedian(count,2);
    for i = 1:nbin        
        CI(i,1) = prctile(count(i,:),2.5);
        CI(i,2) = prctile(count(i,:),97.5);
    end

    
    return
end


        













