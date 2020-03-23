function [fitParams,minus_sum_log,guessIn] = fitModel(fitConds)
bnds        = fitConds.bnds;
InitRange   = fitConds.InitRange;
if fitConds.iStage == 1 % initial parameter randomization
    guessIn = NaN(1,4);
    for iParam = 1:4
        guessIn(iParam) = rand*(InitRange{iParam}(2)-InitRange{iParam}(1)) + InitRange{iParam}(1);
    end
    [fitParams,minus_sum_log] = fminsearchbnd(@(guessIn) mleCriterion(guessIn,fitConds),...
        guessIn,[bnds{1}(1) bnds{2}(1) bnds{3}(1) bnds{4}(1)],[bnds{1}(2) bnds{2}(2) bnds{3}(2) bnds{4}(2)],fitConds.options);
else
    for iCall = 1:fitConds.nCallfminsearch
        switch iCall
            case 1
                guessIn = fitConds.BestParams;
                fitParams_1Call = fminsearchbnd(@(guessIn) mleCriterion(guessIn,fitConds),...
                    guessIn,[bnds{1}(1) bnds{2}(1) bnds{3}(1) bnds{4}(1)],[bnds{1}(2) bnds{2}(2) bnds{3}(2) bnds{4}(2)],fitConds.options);
            case 2                
                [fitParams,minus_sum_log] = fminsearchbnd(@(fitParams_1Call) mleCriterion(fitParams_1Call,fitConds),...
                    fitParams_1Call,[bnds{1}(1) bnds{2}(1) bnds{3}(1) bnds{4}(1)],[bnds{1}(2) bnds{2}(2) bnds{3}(2) bnds{4}(2)],fitConds.options);
        end
    end
end
end

function value = mleCriterion(guessIn,fitConds)
% parameters = {sigma_m, mu_0, sigma_0, kappa_WMdecay}
sigma_m         = guessIn(1);
mu_0            = guessIn(2);
sigma_0         = guessIn(3);
kappa_WMdecay   = guessIn(4);
maxTB           = fitConds.maxTB;
%
fitConds.sigma_m        = sigma_m;
fitConds.mu_0           = mu_0;
fitConds.sigma_0        = sigma_0;
fitConds.kappa_WMdecay  = kappa_WMdecay;
%
[~,sum_log] = get_Lh_original(fitConds);
%
value = -sum_log;
fprintf([fitConds.imodelname ' model, iStage: %d, iIter: %d/%d, iSub: %d, maxTB: %d, ' ...
    '-LL: %.4f, sigma_m: %.4f, mu_0: %.4f, sigma_0: %.4f, kappa_WMdecay: %.4f\n'],...
    fitConds.iStage, fitConds.iIter, fitConds.nIter, fitConds.iSub, maxTB, value, sigma_m, mu_0, sigma_0, kappa_WMdecay)

end