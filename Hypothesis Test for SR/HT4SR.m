function p_value = HT4SR(r1, r2, type)
    if length(r1) ~= length(r2)
        p_value = 'Different Length!';
    else
        t = length(r1);
        mu1 = mean(r1);
        mu2 = mean(r2);
        sd1 = std(r1);
        sd2 = std(r2);
        corr = corrcoef(r1, r2);
        sr1 = mu1 / sd1;
        sr2 = mu2 / sd2;
        
        if strcmp(type, 'Memmel03')
            theta = 1/t * (2 * (1 - corr) + 1/2 * (sr1^2 + sr2^2 - 2 * sr1 * sr2 * corr^2));
            z = (sr1 - sr2) / sqrt(theta);
        elseif strcmp(type, 'JK81')
            theta = 1/t * (2 * sd1^2 * sd2^2 * (1 - corr) + 1/2 * mu1^2 * sd1^2 + 1/2 * mu2^2 * sd2^2 - mu1 * mu2 / (2 * sd1 * sd2) * (sd1^2 * sd2^2 * (1 + corr^2)));
            z = (mu1 * sd2 - mu2 * sd1) / sqrt(theta);
        end
        
        p_value = 1 - tcdf(z, t - 1);
    end
end