function createRetinotopicMaps(flywheelName, subjID, session, pulseOxDir, stimulusFile, varargin)
% Script that downloads the analysis products from flywheel into the
% directory specified in the local hook file and creates eccentricity and
% polar angle maps for a the average of functional runs.
%
% Syntax:
%  createRetinotopicMaps(flywheelName, subjID, session, pulseOxDir, stimulusFile)
%
% Description:
%   This routine connects to a Flywheel server and downloads fmriprep and
%   neuropythy analysis data. It then runs the TFE IAMP and TFE RPRF on it
%   to create temporal fits of the RETINO data to create eccentricity and polar angle retinotopic maps. 
%
% Inputs:
%   flywheelName            - Define the project. This is the Flywheel project you
%                             are interested in extracting retinotopic
%                             mapping data from.
%
%   subjID                  - Define the subject you are mapping data for.
%                             Must be in the format "sub-'subjectLabel'".
%
%   session                 - Define the session for this subject that you
%                             are interested in. Must be in the format
%                             "ses-'sessionLabel'".
%
%   pulseOxDir              - Location of the directory containing PulseOx
%                             files for all functional runs for this 
%                             subject's session. These are puls.mat files.
%
%   stimulusFile              - Location of the .mat stimulus file.
%
%
% Outputs:
%   none
%
% Optional key/value pairs:
%   'ANTSwarp'                - Logical, default true
%
% Examples:
%{
    createRetinotopicMaps('tome', 'sub-TOME3016', 'ses-Session2', '~/Documents/PulseOx_TOME_3016_Sess2/Dicom', '~/Documents/flywheel/retinotopyTOMEAnalysis/pRFpacket_10x10.mat')
%}

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = false;

p.addRequired('flywheelName', @ischar);
p.addRequired('subjID', @ischar);
p.addRequired('session', @ischar);
p.addRequired('pulseOxDir', @ischar);
p.addRequired('stimulusFile', @ischar);

p.addParameter('ANTSwarp',true, @islogical);

p.parse(flywheelName, subjID, session, pulseOxDir, stimulusFile, varargin{:});

% Distribute into variables
flywheelName = p.Results.flywheelName;
subjID = p.Results.subjID;
session = p.Results.session;
pulseOxDir = p.Results.pulseOxDir;
stimulusFile = p.Results.stimulusFile;
ANTSwarp = p.Results.ANTSwarp;

%% Open flywheel object
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

%% Analysis that we are using data from 
fmriprepLabel = 'fmriprep 08/01/2018 18:47:45';
neuropythyLabel = 'retinotopy-templates 03/15/2018 11:34:56';

%% Check for/make the project level directory
projectDir = getpref('retinotopyTOMEAnalysis','projectRootDir');
if ~exist(projectDir)
    mkdir(projectDir)
end

%% Set up fmriprep directory and download
fmriprepDir = fullfile(projectDir,'fmriprep');
if ~exist(fmriprepDir)
    mkdir(fmriprepDir)
end

% download the data
[fwInfo] = getAnalysisFromFlywheel(flywheelName,fmriprepLabel,fmriprepDir, 'verbose', false, 'searchDir', projectDir);

%% make session dir
sessionDate = fwInfo.timestamp(1:10);
sessionFolder = [fwInfo.subject,'_',sessionDate];
sessionDir = fullfile(projectDir,sessionFolder);
if ~exist(sessionDir)
    mkdir(sessionDir)
end

%% Move contents of fmriprep into session directory
dataLocation = fullfile(fmriprepDir, [fwInfo.subject, '_', fwInfo.analysis_id],fwInfo.analysis_id);

% move fmriprep output
if exist(fullfile(dataLocation,'fmriprep'))
    movefile(fullfile(dataLocation,'fmriprep'), sessionDir)
else
    warning('retinotopyTOMEAnalysis:analysisDirectoryMissing','WARNING: fmriprep directory missing from %s folder',fwInfo.analysis_id);
end

% move freesurfer output
if exist(fullfile(dataLocation,'freesurfer'))
    movefile(fullfile(dataLocation,'freesurfer'), sessionDir)
else
    warning('retinotopyTOMEAnalysis:analysisDirectoryMissing','WARNING: freesurfer directory missing from %s folder',fwInfo.analysis_id);
end

%% Make bookkeping files

% create bookkeeping file so that the download function can find the
% analysis id after the clean up at the end of this script. This allows
% for download skip

% fmriprep file
cmd = ['echo -e "' 'This is a bookkeeping file.\nDo not delete unless you delete the whole analysis directory.\nIf you want to download the data again and you got both frmiprep and freesurfer from the same gear, then you need to dete both analysis directories" > ' ...
    fullfile(sessionDir,'fmriprep',[fwInfo.analysis_id '.txt'])];
system(cmd);

% freesurfer file
cmd = ['echo -e "' 'This is a bookkeeping file.\nDo not delete unless you delete the whole analysis directory.\nIf you want to download the data again and you got both frmiprep and freesurfer from the same gear, then you need to dete both analysis directories" > ' ...
    fullfile(sessionDir,'freesurfer',[fwInfo.analysis_id '.txt'])];
system(cmd);


%% Make Benson Gear (Neuropythy) directory
bensonDir = fullfile(sessionDir,'neuropythy');
if ~exist(bensonDir)
    mkdir(bensonDir)
end

%% Download Benson Gear (Neuropythy) output
clear fwInfo
[fwInfo] = getAnalysisFromFlywheel(flywheelName,neuropythyLabel,bensonDir, 'verbose', false, 'searchDir', projectDir);

%% Benson bookkeeping file
%
% freesurfer file
cmd = ['echo -e "' 'This is a bookkeeping file.\nDo not delete unless you delete the whole analysis directory." > ' ...
    fullfile(sessionDir,'neuropythy',[fwInfo.analysis_id '.txt'])];
system(cmd);


%% Clean Up
% remove the unused files from the fmriprep download
cmd = ['rm -rf ' fmriprepDir];
system(cmd);

%% Now to identify the RETINO runs

% Set up necessary path to functional runs
functionalPath = fullfile(sessionDir, 'fmriprep', subjID, session, 'func');

% Get all files and filenames
allFiles = dir(strcat(functionalPath,'/*'));
fileNames = {allFiles.name};

% Set up variables for later
functionalRuns = {};
functionalMasks = {};
confoundFiles = {};

% Now only filter out only RETINO files
for ii = 1:length(fileNames)
    if contains(fileNames{ii},'RETINO')
        if contains(fileNames{ii},'preproc')
            functionalRuns{end+1} = fileNames{ii};
        elseif contains(fileNames{ii},'brainmask')
            functionalMasks{end+1} = fileNames{ii};
        elseif contains(fileNames{ii},'.tsv')
            confoundFiles{end+1} = fileNames{ii};
        end
    end
end

fullFileConfounds = fullfile(functionalPath,confoundFiles);

%% Set up GM volume and V1 Mask for filtering data

% We have to warp the GM mask to MNI EPI space
% Native GM brainmask to warp to MNI EPI space
inFileName = strcat(subjID,'_',session,'_T1w_class-GM_probtissue.nii.gz');

% Reference functional brainmask in MNI EPI space
refFileName = functionalMasks{1};

% Warp file name (product of running fmriprep)
warpFileName = strcat(subjID,'_',session,'_T1w_target-MNI152NLin2009cAsym_warp.h5');

% Path to warp file and GM file
GMFilePath   = fullfile(sessionDir, 'fmriprep', subjID, session, 'anat');
retinoPath = fullfile(sessionDir,'neuropythy');

inFile = fullfile(GMFilePath,inFileName);
outFile = fullfile(GMFilePath,'GM_brainmask_MNI_EPI.nii.gz');
warpFile = fullfile(GMFilePath,warpFileName);
refFile = fullfile(functionalPath,refFileName);

% Apply the ANTS transform
if ANTSwarp
    applyANTsWarpToData(inFile,outFile,warpFile,refFile);
end

%% Extract Signal from Voxels

% make full file path to functional runs
functionalRunFiles = fullfile(functionalPath,functionalRuns);

% extract the mean signal from voxels
mask = MRIread(outFile);
maskVol = mask.vol;
maskVol(round(size(maskVol,1)/2):end,:,:) = 0;
[voxelTimeSeries, voxelIndex] = extractTimeSeriesFromMask(functionalRunFiles,maskVol,'threshold', 0.5);
dim1 = voxelIndex(:,1);
dim2 = voxelIndex(:,2);
dim3 = voxelIndex(:,3);

% Clip first two data points from the time series, as the signal is not yet
% steady state. We need to do more to either prevent these volumes from
% being saved, or to automatically detect this condition.
voxelTimeSeries = voxelTimeSeries(:,1:end,:);

%% IAMP TFE Confound Regression

% Construct model object
temporalFit = tfeIAMP('verbosity','none');

% Define the TR
    TR = 0.800;

        % make stimulus timebase
    deltaT = 800;
    totalTime = 420*deltaT;
    stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
    responseStruct.timebase = stimulusStruct.timebase;
    thePacket.kernel = [];
    thePacket.metaData = [];

%% Create a cell of stimulusStruct (one struct per run)

cleanFunctionalRuns = {};

for jj = 1:length(functionalRuns)
    
    % get confound regressors 
    confoundRegressors = getConfoundRegressors(fullFileConfounds{jj});
    load(strcat(pulseOxDir,'/run',num2str(jj),'/puls.mat'));
    confoundRegressors(:,end+1:(end+size(output.all,2)-4)) = output.all(:,5:8);

    % normalize the regressors
    confoundRegressors = confoundRegressors - nanmean(confoundRegressors);
    confoundRegressors = confoundRegressors ./ nanstd(confoundRegressors);

    % This is a hack to set the initial value of the stdVar (or perhaps
    % framewise displacement vector) to zero instead of nan.
    confoundRegressors(1,3)=0;
    
    thePacket.stimulus.values = confoundRegressors';
    
    defaultParamsInfo.nInstances = size(thePacket.stimulus.values,1);

    % get the data for all masked voxel in a run 
    runData = voxelTimeSeries(:,:,jj);
    
    % convert to percent signal change relative to the mean
    voxelMeanVec = mean(runData,2);
    PSC = 100*((runData - voxelMeanVec)./voxelMeanVec);

    % timebase will be the same for every voxel
    thePacket.response.timebase = stimulusStruct.timebase;
    thePacket.stimulus.timebase = stimulusStruct.timebase;

    % loop over voxels --> returns a "cleaned" time series
    tic
    for vxl = 1:size(PSC,1)
        % place time series from this voxel into the packet
        thePacket.response.values = PSC(vxl,:);
        
        % TFE linear regression here
%         disp(jj);
%         disp(vxl);
        [paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(thePacket,...
            'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression','errorType','1-r2');
        
        confoundBetas(:,vxl) = paramsFit.paramMainMatrix;
        cleanRunData(vxl,:) = thePacket.response.values - modelResponseStruct.values;
        IAMPfval(vxl) = 1-fVal;
    end
    toc
    cleanVolTemp = cleanRunData(:,:);
    IAMPfvalVol = nan(size(maskVol));
    cleanRun = MRIread(fullfile(functionalPath,functionalRuns{jj}));
    cleanVol = cleanRun.vol;
    IAMPfile = mask;
    for aa = 1:length(dim1)
        cleanVol(dim1(aa),dim2(aa),dim3(aa),:) = cleanVolTemp(aa);
        IAMPfvalVol(dim1(aa), dim2(aa), dim3(aa)) = IAMPfval(aa);
    end
    
    cleanRun.vol = cleanVol;
    cleanRun.fspec = strcat(erase(functionalRuns{jj},'.nii.gz'),'_clean.nii.gz');
    MRIwrite(cleanRun, fullfile(functionalPath,cleanRun.fspec));
    
    IAMPfile.vol = IAMPfvalVol;
    IAMPfile.fspec = strcat(erase(functionalRuns{jj},'.nii.gz'),'_fval.nii.gz');
    MRIwrite(IAMPfile,fullfile(functionalPath,IAMPfile.fspec));
    
    cleanFunctionalRuns{end+1} = cleanRun.fspec;

end

 %% Construct the model object
tfeHandle = tfeRPRF('verbosity','none');

%% Temporal definition of the stimulus and response
deltaT = 800; % in msecs
totalTime = 420*deltaT; % in msecs

% Load stimulus
load(stimulusFile);

%% Set up TFE

% Define a kernelStruct. In this case, a double gamma HRF
hrfParams.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in secs)
hrfParams.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in secs)
hrfParams.gammaScale = 10; % scaling factor between the positive and negative gamma componenets
kernelStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);

% The timebase is converted to seconds within the function, as the gamma
% parameters are defined in seconds.
hrf = gampdf(kernelStruct.timebase/1000, hrfParams.gamma1, 1) - ...
gampdf(kernelStruct.timebase/1000, hrfParams.gamma2, 1)/hrfParams.gammaScale;
kernelStruct.values=hrf;

% prepare this kernelStruct for use in convolution as a BOLD HRF
kernelStruct.values=kernelStruct.values-kernelStruct.values(1);
kernelStruct=normalizeKernelArea(kernelStruct);

% Get the default forward model parameter
defaultParamsInfo.nInstances=1;
params0 = tfeHandle.defaultParams('defaultParamsInfo', defaultParamsInfo);

% start the packet assembly
thePacket.stimulus = stimulus;
thePacket.kernel = kernelStruct;
thePacket.metaData = [];
thePacket.response.timebase = thePacket.kernel.timebase;

%% Loop through RETINO runs and perform TFE commands

% Load donor
donorFilePath = fullfile(functionalPath, functionalMasks{1});
donorX = MRIread(donorFilePath);
donorY = donorX;
donorSigma = donorX;
donorFval = donorX;
donorA = donorX;
donorFit = donorX;

% File to write median volume to
writeFileNameX = 'RETINO_average_map_x_new.nii.gz';
writeFileNameY = 'RETINO_average_map_y_new.nii.gz';
writeFileNameSigma = 'RETINO_average_map_sigma_new.nii.gz';
writeFileNameFval = 'RETINO_average_map_fval_new.nii.gz';
writeFileNameA = 'RETINO_average_map_amplitude_new.nii.gz';
writeFileNameFit = 'RETINO_average_map_fit_new.nii.gz';

writeFilePathX = fullfile(functionalPath,writeFileNameX);
writeFilePathY = fullfile(functionalPath,writeFileNameY);
writeFilePathSigma = fullfile(functionalPath,writeFileNameSigma);
writeFilePathFval = fullfile(functionalPath,writeFileNameFval);
writeFilePathA = fullfile(functionalPath,writeFileNameA);
writeFilePathFit = fullfile(functionalPath,writeFileNameFit);

% Load GM Brainmask
niiBrainmask = MRIread(outFile);
niiBrainmaskVol = niiBrainmask.vol;

% Load and Create Secondary Constriction Brainmask
niiBrainmaskV1 = MRIread(outFile);
niiBrainmaskV1.vol = nan(size(niiBrainmaskVol));
niiBrainmaskV1.vol(1:58,:,:) = 1;
niiBrainmask.fspec = 'V1_surroundings_brainmask_MNI_EPI.nii.gz';
MRIwrite(niiBrainmaskV1, fullfile(retinoPath,'V1_surroundings_brainmask_MNI_EPI.nii.gz')); 
niiBrainmaskV1Vol = niiBrainmaskV1.vol;


% Load temporary initial Functional Run
filePath = fullfile(functionalPath, cleanFunctionalRuns{1});
    
% Load the file 
nii = MRIread(filePath);
fmri = nii.vol;

% Loop through desired files
for ff = 1:length(cleanFunctionalRuns)
    filePath = fullfile(functionalPath, cleanFunctionalRuns{ff});
    
    % Load the file 
    nii = MRIread(filePath);
    fmri = nii.vol;
    
    if ff==1
        avgVol = fmri;
        % subplot(1,4,ff)
%         plot(squeeze(avgVol(42,25,41,:)));
    else
        avgVol = avgVol+fmri;
        % subplot(1,4,ff)
%         plot(squeeze(avgVol(42,25,41,:)));
    end
end
avgVol = avgVol ./ length(cleanFunctionalRuns);

% Identify which parts of the GM brainmask are the "brain"
indices = find(niiBrainmaskVol > 0.1 & niiBrainmaskV1Vol == 1);
[x,y,z] = ind2sub(size(niiBrainmaskVol),indices); 

% Create empty volume to fill
newVolX = nan(size(niiBrainmaskVol));
newVolY = nan(size(niiBrainmaskVol));
newVolSigma = nan(size(niiBrainmaskVol));
newVolFval = nan(size(niiBrainmaskVol));
newVolA = nan(size(niiBrainmaskVol));
FitVol = nan(size(avgVol));

% Time to fill it in with a median map
tic
for idx = 1:1
    disp(idx);
    % This is the step where you grab the time series for a voxel
    
%     if idx==1
%         avgTS = avgVol(x(idx),y(idx),z(idx),:);
%     else
%         avgTS = avgTS+avgVol(x(idx),y(idx),z(idx),:);
%     end
%     
%     
%     avgTS = avgTS ./ length(indices);
%     avgTS = squeeze(avgTS); 
%     plot(avgTS);

    ts = avgVol(x(idx),y(idx),z(idx),:);
    ts = squeeze(ts)';
    ts = ts ./ mean(ts);
    ts = ts - mean(ts);
    thePacket.response.values = ts;
    
    % This is the call to the fitting engine
    [paramsFit,fVal,modelResponseStruct] = ...
    tfeHandle.fitResponse(thePacket,...
    'defaultParamsInfo', defaultParamsInfo,'errorType','1-r2');

%    Store the parameter value
    newVolX(x(idx),y(idx),z(idx))= paramsFit.paramMainMatrix(1);
    newVolY(x(idx),y(idx),z(idx))= paramsFit.paramMainMatrix(2);
    newVolSigma(x(idx),y(idx),z(idx))= paramsFit.paramMainMatrix(3);
    newVolFval(x(idx),y(idx),z(idx))= 1-fVal;
    newVolA(x(idx),y(idx),z(idx))= paramsFit.paramMainMatrix(4);
    FitVol(x(idx),y(idx),z(idx),:)=modelResponseStruct.values;
    
end
toc

% Write the new maps to the writeFiles 50,20,43 (21, 51, 44) | 
donorX.vol = newVolX;
donorX.fspec = writeFilePathX;
MRIwrite(donorX,writeFilePathX);

donorY.vol = newVolY;
donorY.fspec = writeFilePathY;
MRIwrite(donorY,writeFilePathY);

donorSigma.vol = newVolSigma;
donorSigma.fspec = writeFilePathSigma;
MRIwrite(donorSigma,writeFilePathSigma);

donorFval.vol = newVolFval;
donorFval.fspec = writeFilePathFval;
MRIwrite(donorFval,writeFilePathFval);

donorA.vol = newVolA;
donorA.fspec = writeFilePathA;
MRIwrite(donorA,writeFilePathA);

donorFit.vol = FitVol;
donorFit.fspec = writeFilePathFit;
MRIwrite(donorFit,writeFilePathFit);

%% ANTS Warp these files to Native T1 space

inFile = writeFilePathFval;
outFile = fullfile(functionalPath,'RETINO_average_map_fval_native.nii.gz');
warpFile = fullfile(GMFilePath,strcat(subjID,'_',session,'_T1w_space-MNI152NLin2009cAsym_target-T1w_warp.h5'));
refFile = fullfile(GMFilePath, strcat(subjID,'_',session,'_T1w_class-GM_probtissue.nii.gz'));

applyANTsWarpToData(inFile,outFile, warpFile, refFile);

inFile = writeFilePathX;
outFile = fullfile(functionalPath,'RETINO_average_map_fval_native.nii.gz');
warpFile = fullfile(GMFilePath,strcat(subjID,'_',session,'_T1w_space-MNI152NLin2009cAsym_target-T1w_warp.h5'));
refFile = fullfile(GMFilePath, strcat(subjID,'_',session,'_T1w_class-GM_probtissue.nii.gz'));

applyANTsWarpToData(inFile,outFile, warpFile, refFile);

inFile = writeFilePathY;
outFile = fullfile(functionalPath,'RETINO_average_map_fval_native.nii.gz');
warpFile = fullfile(GMFilePath,strcat(subjID,'_',session,'_T1w_space-MNI152NLin2009cAsym_target-T1w_warp.h5'));
refFile = fullfile(GMFilePath, strcat(subjID,'_',session,'_T1w_class-GM_probtissue.nii.gz'));

applyANTsWarpToData(inFile,outFile, warpFile, refFile);

inFile = writeFilePathSigma;
outFile = fullfile(functionalPath,'RETINO_average_map_fval_native.nii.gz');
warpFile = fullfile(GMFilePath,strcat(subjID,'_',session,'_T1w_space-MNI152NLin2009cAsym_target-T1w_warp.h5'));
refFile = fullfile(GMFilePath, strcat(subjID,'_',session,'_T1w_class-GM_probtissue.nii.gz'));
applyANTsWarpToData(inFile,outFile, warpFile, refFile);

inFile = writeFilePathA;
outFile = fullfile(functionalPath,'RETINO_average_map_fval_native.nii.gz');
warpFile = fullfile(GMFilePath,strcat(subjID,'_',session,'_T1w_space-MNI152NLin2009cAsym_target-T1w_warp.h5'));
refFile = fullfile(GMFilePath, strcat(subjID,'_',session,'_T1w_class-GM_probtissue.nii.gz'));

applyANTsWarpToData(inFile,outFile, warpFile, refFile);

fvalVol = MRIread(fullfile(functionalPath,'RETINO_average_map_fval_native.nii.gz'));
xVol = MRIread(fullfile(functionalPath,'RETINO_average_map_x_native.nii.gz'));
yVol = MRIread(fullfile(functionalPath,'RETINO_average_map_y_native.nii.gz'));
sigmaVol = MRIread(fullfile(functionalPath,'RETINO_average_map_sigma_native.nii.gz'));
ampVol = MRIread(fullfile(functionalPath,'RETINO_average_map_amplitude_native.nii.gz'));

writeFileNameEccen = 'RETINO_average_map_eccen_new.nii.gz';
writeFileNamePA = 'RETINO_average_map_angle_new.nii.gz';



EccenVol = nan(size(xVol.vol));
PAVol = nan(size(xVol.vol));

testIndices = find(fvalVol.vol < 0.00);
[x1,y1,z1] = ind2sub(size(fvalVol.vol),testIndices); 
for iii = 1:length(testIndices)
    xVol.vol(x1(iii),y1(iii),z1(iii)) = 0;
    yVol.vol(x1(iii),y1(iii),z1(iii)) = 0;
    sigmaVol.vol(x1(iii),y1(iii),z1(iii)) = 0;
    ampVol.vol(x1(iii),y1(iii),z1(iii)) = 0;
end

testIndices2 = find(fvalVol.vol > 0);
[x2,y2,z2] = ind2sub(size(fvalVol.vol),testIndices2); 
for jj = 1:length(testIndices2)
    xDist = xVol.vol(x2(jj), y2(jj), z2(jj))-5;
    yDist = yVol.vol(x2(jj), y2(jj), z2(jj))-5;
    
    EccenVol(x2(jj), y2(jj), z2(jj)) = sqrt((xDist)^2 + (yDist)^2);
    
    PAVol(x2(jj), y2(jj), z2(jj)) = rad2deg(atan(abs(xDist)/abs(yDist)));
    if yDist < 0 && xDist > 0
        PAVol(x2(jj), y2(jj), z2(jj)) = 180 - PAVol(x2(jj), y2(jj), z2(jj));
    elseif yDist < 0 && xDist < 0
        PAVol(x2(jj), y2(jj), z2(jj)) = -180 + PAVol(x2(jj), y2(jj), z2(jj));
    elseif yDist > 0 && xDist > 0 
        PAVol(x2(jj), y2(jj), z2(jj)) = PAVol(x2(jj), y2(jj), z2(jj));
    else
        PAVol(x2(jj), y2(jj), z2(jj)) = -90 + PAVol(x2(jj), y2(jj), z2(jj));
    end
end

donorEccen = xVol;
donorPA = xVol;

donorEccen.vol = EccenVol;
donorPA.vol = PAVol;

writeFilePathEccen = fullfile(functionalPath,writeFileNameEccen);
writeFilePathPA = fullfile(functionalPath,writeFileNamePA);

MRIwrite(donorEccen,writeFilePathEccen);
MRIwrite(donorPA,writeFilePathPA);

MRIwrite(xVol,fullfile(functionalPath,'RETINO_average_map_x_native.nii.gz'));
MRIwrite(yVol,fullfile(functionalPath,'RETINO_average_map_y_native.nii.gz'));
MRIwrite(sigmaVol,fullfile(functionalPath,'RETINO_average_map_sigma_native.nii.gz'));
MRIwrite(ampVol,fullfile(functionalPath,'RETINO_average_map_amplitude_native.nii.gz'));