clear all; close all; clc;
modelnames          = {'NBMC'}; % ["NBMC", or "constant"]
subjects            = 1:18;  
stages              = [1 2]; % stage 1: initial parameter radomization, stage2: fitting until parameters are converged
%%
nModel              = length(modelnames);
allbnds             = {{[10^-5 5],[-5 5],[10^-5 100],[10^-5 5]},...
                       {[10^-5 5],[-5 5],[0 0],[0 0]}};  
% parameters = {sigma_m, mu_0, sigma_0, kappa_WMdecay}   
bnds = cell(nModel,1);
for iModel = 1:nModel
    iname = modelnames{iModel};
    if strcmp(iname(1),'N') == 1 % for NBMC 
        bnds{iModel} = allbnds{1};
    elseif strcmp(iname(1),'c') == 1 % for the constant criterion model
        bnds{iModel} = allbnds{2};
    end
end
maxTrialBack        = 7;
InitRange           = bnds;
thrRT               = 0.3;
nIter4InitRandomize = 20;
nBest4SecondStage   = 20;
dir0 = pwd;
addpath([dir0 '/Sub'])
idir = [dir0 '/results/'];
if isempty(dir(idir)) == 1
    mkdir(idir)
end
cd(idir)
load([dir0 '/maindata.mat'])
fitConds.thrRT = thrRT;
for iStage = stages
    fitConds.iStage = iStage;
    switch iStage
        case 1
            stagename                   = 'First';
            fitConds.options            = optimset('MaxFunEvals',50,'MaxIter',50);
            fitConds.nCallfminsearch    = 1;
        case 2
            stagename                   = 'Second';
            fitConds.stagename          = 'SecondStage';
            fitConds.options            = optimset('MaxFunEvals',10^5,'MaxIter',10^5,'TolFun',10^-7,'TolX',10^-7);
            fitConds.nCallfminsearch    = 2;
    end
    for imodel = 1:nModel
        switch iStage
            case 1
                fitConds.InitRange           = InitRange{imodel};
                nIter                        = nIter4InitRandomize;
                fitConds.nIter               = nIter;
                fitConds.nIter4InitRandomize = nIter;
            case 2
                nIter                       = nBest4SecondStage;
                fitConds.nIter              = nIter;
                fitConds.nBest4SecondStage  = nIter;
        end
        imodelname          = modelnames{imodel};
        fitConds.imodel     = imodel;
        fitConds.bnds       = bnds{imodel};
        fitConds.maxTB      = maxTrialBack;
        fitConds.imodelname = imodelname;
        idir = ['./' imodelname '/' stagename 'Stage/' ];
        jdir = ['./' imodelname '/FirstStage/'];
        if isempty(dir(idir)) == 1
            mkdir(idir)
        end
        for iSub = subjects
            imat1 = isempty(dir(['./' idir '/' num2str(iSub) '.mat']));
            imat2 = isempty(dir(['./' idir '/' num2str(iSub) '_computing.mat']));
            if imat1*imat2 == 1
                save([idir '/' num2str(iSub) '_computing.mat'],'iSub')
                if iStage == 2
                    load([jdir '/' num2str(iSub) '.mat'],'fitResults')
                end
                fitConds.iSub   = iSub;
                fitConds.Z      = Stm{iSub};
                fitConds.Chc    = Chc{iSub};
                fitConds.RT     = RT{iSub};
                %
                fit_sigma_m         = NaN(nIter,1);
                fit_mu_0            = NaN(nIter,1);
                fit_sigma_0         = NaN(nIter,1);
                fit_kappa_WMdecay   = NaN(nIter,1);
                %
                init_sigma_m        = NaN(nIter,1);
                init_mu_0           = NaN(nIter,1);
                init_sigma_0        = NaN(nIter,1);
                init_kappa_WMdecay  = NaN(nIter,1);
                %
                minus_sum_log_Lh = NaN(nIter,1);
                for iIter = 1:nIter
                    fitConds.iIter = iIter;
                    %
                    if iStage == 2
                        fitConds.BestParams = ...
                            [fitResults.fit_sigma_m(iIter) fitResults.fit_mu_0(iIter) ...
                            fitResults.fit_sigma_0(iIter) fitResults.fit_kappa_WMdecay(iIter)];
                    end
                    [ifitParams,iminus_sum_log_Lh,guessIn] = fitModel(fitConds);
                    % parameters = {sigma_m, mu_0, sigma_0, kappa_WMdecay}
                    fit_sigma_m(iIter)          = ifitParams(1);
                    fit_mu_0(iIter)             = ifitParams(2);
                    fit_sigma_0(iIter)          = ifitParams(3);
                    fit_kappa_WMdecay(iIter)    = ifitParams(4);
                    %
                    init_sigma_m(iIter)         = guessIn(1);
                    init_mu_0(iIter)            = guessIn(2);
                    init_sigma_0(iIter)         = guessIn(3);
                    init_kappa_WMdecay(iIter)   = guessIn(4);
                    %
                    minus_sum_log_Lh(iIter) = iminus_sum_log_Lh;
                end
                %
                [minus_sum_log_Lh,sInd]         = sort(minus_sum_log_Lh);
                fitResults                      = [];
                fitResults.minus_sum_log_Lh     = minus_sum_log_Lh;
                fitResults.fit_sigma_m          = fit_sigma_m(sInd);
                fitResults.fit_mu_0             = fit_mu_0(sInd);
                fitResults.fit_sigma_0          = fit_sigma_0(sInd);
                fitResults.fit_kappa_WMdecay    = fit_kappa_WMdecay(sInd);
                %
                fitConds.init_sigma_m           = init_sigma_m(sInd);
                fitConds.init_mu_0              = init_mu_0(sInd);
                fitConds.init_sigma_0           = init_sigma_0(sInd);
                fitConds.init_kappa_WMdecay     = init_kappa_WMdecay(sInd);
                %
                plot_BestParams(fitConds,fitResults,idir)
                %
                save([idir '/' num2str(iSub) '.mat'],'fitConds','fitResults')
                delete([idir '/' num2str(iSub) '_computing.mat'])
            end
        end
    end
end