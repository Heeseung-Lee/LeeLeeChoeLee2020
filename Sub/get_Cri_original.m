function [mu_p_c,sigma_p_c,sigma_c] = get_Cri_original(fitConds)
% parameters = {sigma_m, mu_0, sigma_0, kappa_WMdecay}
sigma_m         = fitConds.sigma_m;
mu_0            = fitConds.mu_0;
sigma_0         = fitConds.sigma_0;
kappa_WMdecay   = fitConds.kappa_WMdecay;
%
Z           = fitConds.Z;
maxTB       = fitConds.maxTB;
nT          = size(Z,1);
nR          = size(Z,2);
mu_p_c      = NaN(size(Z));
sigma_p_c   = NaN(size(Z));
sigma_c     = NaN(size(Z));
%
for iT = 2:nT
    % exponential working memory decay
    if iT >= 2 && iT <= (maxTB+1)
        pre_sigma_m     = sigma_m*fliplr((1+kappa_WMdecay).^(1:(iT-1)));
    else
        pre_sigma_m     = sigma_m*fliplr((1+kappa_WMdecay).^(1:maxTB));
    end
    %
    if iT >= 2 && iT <= (maxTB+1)
        pre_Z       = Z(1:iT-1,:);
    else
        pre_Z       = Z((iT-maxTB):iT-1,:);
    end
    Sigmas          = [sigma_0 pre_sigma_m];
    Wts             = Sigmas.^-2/sum(Sigmas.^-2);
    imu_p_c         = Wts*[mu_0*ones(1,nR); pre_Z];
    isigma_p_c      = sqrt(sum(Wts(2:end).^2.*pre_sigma_m.^2));
    isigma_c        = sqrt(sum(Sigmas.^2.*Wts.^2));
    mu_p_c(iT,:)    = imu_p_c;
    sigma_p_c(iT,:) = isigma_p_c;
    sigma_c(iT,:)   = isigma_c;
end