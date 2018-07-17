%% analyzeLFContrast
%
% This script calls function in order to analyze the data for the
% LFContrast experiment.

%% Convenience variables
projectName  = 'LFContrastAnalysis';
flywheelName = 'LFContrast';
subjID       = 'sub-HEROgka1';
session      = 'ses-0411181853PM';
sessionLabel = '04/11/18 18:53 PM';

%% Analysis labels that we are going to go and get
fmriprepLabel   = 'fmriprep 04/12/2018 15:16:06';
neuropythyLabel = 'retinotopy-templates 04/13/2018 16:46:22';
fwInfo          = getAnalysisFromFlywheel(flywheelName,fmriprepLabel,fullfile(getpref('LFContrastAnalysis','projectRootDir'),'fmriprep'),'nodownload',true);
sessionDir      = fullfile(getpref('LFContrastAnalysis','projectRootDir'),[fwInfo.subject,'_', fwInfo.timestamp(1:10)]);
if ~isfolder(sessionDir)
    mkdir(sessionDir);
end
if ~isfolder(fullfile(sessionDir,'fmriprep'))
    fwInfo          = getAnalysisFromFlywheel(flywheelName,fmriprepLabel,fullfile(sessionDir,'fmriprep'));
    fwInfoRetino    = getAnalysisFromFlywheel(flywheelName,neuropythyLabel,fullfile(sessionDir,'neuropythy'));
end
%% Relevant Nifti names for analysis

% functional runs
 functionalRuns = {'sub-HEROgka1_ses-0411181853PM_task-tfMRIFLASHAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastAP_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'};

% brain mask of function run for the reference volume in ANTs step
refFileName  = 'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';

% output files of Neuropythy (retinotopy template)
retinoFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz',};

% warp file name (product of running fmriprep)
warpFileName = 'sub-HEROgka1_ses-0411181853PM_T1w_space-MNI152NLin2009cAsym_target-T1w_warp.h5';

% Set up paths to nifti and .h5 files
retinoPath     = fullfile(sessionDir,'neuropythy');
functionalPath = fullfile(sessionDir, 'fmriprep', [fwInfo.subject '_' fwInfo.analysis_id], fwInfo.analysis_id, 'fmriprep', subjID, session, 'func');
warpFilePath   = fullfile(sessionDir, 'fmriprep', [fwInfo.subject '_' fwInfo.analysis_id], fwInfo.analysis_id, 'fmriprep', subjID, session, 'anat');

%% Create restricted V1 mask

% load ecc nifti file
eccenPos       = find(~cellfun(@isempty,strfind(retinoFiles,'eccen')));
[~,tempName,~] = fileparts(retinoFiles{eccenPos});
[~,outName,~]  = fileparts(tempName);
eccenFileName  = fullfile(retinoPath,[outName '.nii.gz']);
eccen          = MRIread(eccenFileName);

% load areas nifti file
areasPos       = find(~cellfun(@isempty,strfind(retinoFiles,'areas')));
[~,tempName,~] = fileparts(retinoFiles{areasPos});
[~,outName,~]  = fileparts(tempName);
areasFileName  = fullfile(retinoPath,[outName,'.nii.gz']);
areas          = MRIread(areasFileName);

% make mask from the area and eccentricity maps
areaNum     = 1;
eccenRange  = [3 20];
[~,maskSaveName] = makeMaskFromRetino(eccen,areas,areaNum,eccenRange,retinoPath);

%% Apply the warp to the mask and T1 files using ANTs

files2warp = {'HERO_gka1_T1_MNI_resampled.nii.gz'};
for ii = 1:length(files2warp)
    % input file
    inFile = fullfile(retinoPath,files2warp{ii});
    
    % output file
    [~,tempName,~] = fileparts(inFile);
    [~,outName,~] = fileparts(tempName);
    outFile = fullfile(retinoPath,'HERO_gka1_T1.nii.gz');
    
    % reference file
    refFile = fullfile(functionalPath,refFileName);
    
    % warp file
    warpFile = fullfile(warpFilePath,warpFileName);
    if ~exist(outFile)
        applyANTsWarpToData(inFile, outFile, warpFile, refFile);
    else
        [~,fileName,~] = fileparts(outFile);
        display(sprintf('%s already exist in the specified directory',fileName));
    end
end

%% Extract Signal from voxels
% Load mask nifti
maskPos       = find(~cellfun(@isempty,strfind(files2warp,'mask')));
[~,tempName,~] = fileparts(files2warp{maskPos});
[~,tmpName,~] = fileparts(maskSaveName);
[~,outName,~] = fileparts(tmpName);
maskOutFileName = fullfile(retinoPath,[outName '_MNI_resampled.nii.gz']);
mask = MRIread(maskOutFileName);
maskVol = mask.vol;

% make full file path to functional runs
functionalRuns = fullfile(functionalPath,functionalRuns);

% extract the mean signal from voxels
[voxelTimeSeries, voxelIndex] = extractTimeSeriesFromMask(functionalRuns,maskVol);
meanSignal = squeeze(mean(voxelTimeSeries,1));
% convert to percent signal change relative to the mean
meanMat = repmat(mean(meanSignal,1),[size(meanSignal,1),1]);
PSC = 100*((meanSignal - meanMat)./meanMat);


%% Get trial order info:
trialOrderDir = '~/Documents/MELAStuff/session_1';
trialOrderFiles = {'CRF_session_1_scan0.mat', 'CRF_session_1_scan1.mat','CRF_session_1_scan2.mat', 'CRF_session_1_scan3.mat', 'CRF_session_1_scan4.mat', 'CRF_session_1_scan5.mat'};


%% Create a cell of stimulusStruct (one struct per run)  
for jj = 1:length(trialOrderFiles)
    dataParamFile = fullfile(trialOrderDir,trialOrderFiles{jj});
    TR = 0.800;
    expParams = getExpParams(dataParamFile,TR,'hrfOffset', false, 'stripInitialTRs', false);
    [avgPerCond(:,jj), blockOrder(:,jj)] = sortDataByConditions(PSC(:,jj),expParams);   
    
    % make stimulus timebase
    protocolParams =  load(dataParamFile, 'protocolParams');
    totalTime = protocolParams.protocolParams.nTrials * protocolParams.protocolParams.trialDuration * 1000;
    deltaT = 800;
    stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
    
    % make stimulus values 
    % Stim coding: 80% = 1, 40% = 2, 20% = 3, 10% = 4, 5% = 5, 0% = 6; 
    stimulusStruct.values = zeros(totalTime/deltaT,length(unique(expParams(:,3))));
    for kk = 1:size(expParams,1) 
         stimulusStruct.values(expParams(kk,1):expParams(kk,2), expParams(kk,3)) = 1;
    end  
end

%% create contrast response function
zeroCondMean = repmat(avgPerCond(6,:),[size(avgPerCond,2),1]);
CRF = avgPerCond -zeroCondMean;

%%plot stuff

% load the last param file (doesnt matter because all were the same
% data = load(dataParamFile);
% plotTimeCourse(meanSignal,data.block,data.responseStruct);

%plot CRF
figure; hold on
plot([.8,.4,.2,.1,.05,0],mean(CRF,2),'k','LineWidth',2)
plot([.8,.4,.2,.1,.05,0],CRF,'--o')
ylabel('Percent Signal Change (diff from zero cond)')
xlabel('Contrast Level')
legend('Mean of Runs', 'Run 1', 'Run 2', 'Run 3', 'Run 4', 'Run 5', 'Run 6')
axis square
set(gca, 'YGrid', 'on', 'XGrid', 'off')
title('Contrast Response Function')