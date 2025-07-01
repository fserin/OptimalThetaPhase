function pvalFisher=pval_Fisher(pval,u)
%input
%pval:  a Vxn matrix of contrast p-values
pvalFisher =  NaN*ones(size(pval,1),1);
pval(pval==0) = 1.0000e-016;
for v=1:size(pval,1)
    p=squeeze(pval(v,:));
    n = sum(~isnan(p));
    if n>(u-1)
        p=sort(p(~isnan(p)));
        pvalFisher(v) = 1-chi2cdf(-2*sum(log(p(u:n))), 2*(n-u+1));
    end
end