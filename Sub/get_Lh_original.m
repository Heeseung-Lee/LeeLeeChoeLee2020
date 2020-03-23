function [Lh,sum_log_Lh,pL,mu_p_s,sigma_p_s,mu_p_c,sigma_p_c,sigma_s,sigma_c,dof] = get_Lh_original(fitConds)
% Lh:           log-likelihood for each trial
% sum_log_Lh:   sum of the log-likelihoods
% pL:           proportions of large choice for each trial
% mu_p_s:       mean of the sampling distribution of s
% sigma_p_s:    standard deviation of the sampling distribution of s
% mu_p_c:       mean of the sampling distribution of c
% sigma_p_c:    standard deviation of the sampling distribution of c
% sigma_s:      standard deviation of s
% sigma_c:      standard deviation of c
% dof:          degree of freedom

Chc     = fitConds.Chc;
RT      = fitConds.RT;
thrRT   = fitConds.thrRT;
iname   = fitConds.imodelname;
%
if strcmp(iname(1),'c') == 1 % for the constant criterion model
    sigma_p_c       = zeros(size(Chc));
    mu_p_c          = fitConds.mu_0*ones(size(Chc));
    mu_p_s          = fitConds.Z;
    sigma_p_s       = fitConds.sigma_m*ones(size(Chc));
    sigma_s         = sigma_p_s;
    sigma_c         = sigma_p_c;
elseif strcmp(iname(1),'N') == 1 % for the criterion inference model
    [mu_p_c,sigma_p_c,sigma_c]   = get_Cri_original(fitConds);
    sigma_m         = fitConds.sigma_m;
    mu_0            = fitConds.mu_0;
    sigma_0         = fitConds.sigma_0;
    Z               = fitConds.Z;
    Wt_m            = sigma_0^2/(sigma_0^2 + sigma_m^2);
    Wt_0            = sigma_m^2/(sigma_0^2 + sigma_m^2);
    mu_p_s          = Z*Wt_m + mu_0*Wt_0;
    sigma_p_s       = ones(size(mu_p_s))*sigma_m*Wt_m;
    sigma_s         = sqrt((ones(size(mu_p_s))*sigma_m*Wt_m).^2 + sigma_0.^2*Wt_0.^2);
end
pL                  = normcdf((mu_p_s-mu_p_c)./sqrt(sigma_p_s.^2+sigma_p_c.^2));
%
Lh              = NaN(size(Chc));
Lh(Chc == -1)   = 1-pL(Chc == -1);
Lh(Chc == 1)    = pL(Chc == 1);
Lh(RT<thrRT)    = NaN;
Lh(1,:)         = NaN;
sum_log_Lh      = sum(sum(log(Lh(isnan(Lh)==0))));
dof             = sum(isnan(Lh(:))==0);
end
