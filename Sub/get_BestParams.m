function Params = get_BestParams(dir0,modelname)

idir = [dir0 '/figures/Fig3/modelfitting/' modelname '/SecondStage/'];

Params.model    = modelname;
sigma_m         = NaN(18,1);
mu_0            = NaN(18,1);
sigma_0         = NaN(18,1);
kappa_WMdecay   = NaN(18,1);
for iSub = 1:18
    load([idir '/' num2str(iSub)])
    sigma_m(iSub,1)        = fitResults.fit_sigma_m(1);
    mu_0(iSub,1)           = fitResults.fit_mu_0(1);
    sigma_0(iSub,1)        = fitResults.fit_sigma_0(1);
    kappa_WMdecay(iSub,1)  = fitResults.fit_kappa_WMdecay(1);
end
Params.sigma_m         = sigma_m;
Params.mu_0            = mu_0;
Params.sigma_0         = sigma_0;
Params.kappa_WMdecay   = kappa_WMdecay;
end

