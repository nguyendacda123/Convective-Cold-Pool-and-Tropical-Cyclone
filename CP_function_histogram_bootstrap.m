function [freq, CI] = CP_function_histogram_bootstrap(dat,bins,Msize)
    % compute median and confidence interval of median using bootstrapping
    % method
    
    rng(1);
    nbin = length(bins)-1;
    N = length(dat);
    CI = nan(nbin,2);
    count = nan(nbin,Msize);
    parfor k = 1:Msize        
        dat2 = datasample(dat,N);
        for i = 1:nbin
            count(i,k) = nansum(dat2 >= bins(i) & dat2 < bins(i+1));
        end         
    end
    
    freq = nanmedian(count,2);
    for i = 1:nbin        
        CI(i,1) = prctile(count(i,:),2.5);
        CI(i,2) = prctile(count(i,:),97.5);
    end

    
    return
end


        













