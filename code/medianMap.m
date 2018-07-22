%% medianMap
%
% Script that downloads the analysis products from flywheel into the
% directory specified in the local hook file and creates a median 3D volume 
% for a 4D functional run

%% %% Convenience variables
projectName = 'tome';
flywheelName = 'tome';
subjID = 'sub-TOME3016';
session = 'ses-Session2';

%% Analysis that we are using data from 
fmriprepLabel = 'fmriprep 03/14/2018 13:00:10';
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
    warning('LFContrastAnalysis:analysisDirectoryMissing','WARNING: fmriprep directory missing from %s folder',fwInfo.analysis_id);
end

% move freesurfer output
if exist(fullfile(dataLocation,'freesurfer'))
    movefile(fullfile(dataLocation,'freesurfer'), sessionDir)
else
    warning('LFContrastAnalysis:analysisDirectoryMissing','WARNING: freesurfer directory missing from %s folder',fwInfo.analysis_id);
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
allFiles = dir(strcat(functionalPath,'/*.nii.gz'));
fileNames = {allFiles.name};

% Set up variables for later
functionalRuns = {};
functionalMasks = {};

% Now only filter out only RETINO files
for ii = 1:length(fileNames)
    if contains(fileNames{ii},'FLASH')
        if contains(fileNames{ii},'preproc')
            functionalRuns{end+1} = fileNames{ii};
        elseif contains(fileNames{ii},'brainmask')
            functionalMasks{end+1} = fileNames{ii};
        end
    end
end

%% Create restricted V1 mask

retinoFiles = {'tome_3016_native.template_angle.nii.gz','tome_3016_native.template_areas.nii.gz','tome_3016_native.template_eccen.nii.gz',};
retinoPath = fullfile(sessionDir,'neuropythy');

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

%% Set up GM volume and V1 Mask for filtering data

% We have to warp the GM mask to MNI EPI space
% Native GM brainmask to warp to MNI EPI space
inFileName = 'sub-TOME3016_ses-Session2_T1w_class-GM_probtissue.nii.gz';

% Reference functional brainmask in MNI EPI space
refFileName = functionalMasks{1};

% Warp file name (product of running fmriprep)
warpFileName = 'sub-TOME3016_ses-Session2_T1w_target-MNI152NLin2009cAsym_warp.h5';

% Path to warp file and GM file
GMFilePath   = fullfile(sessionDir, 'fmriprep', subjID, session, 'anat');

inFile = fullfile(GMFilePath,inFileName);
outFile = fullfile(GMFilePath,'GM_brainmask_MNI_EPI.nii.gz');
warpFile = fullfile(GMFilePath,warpFileName);
refFile = fullfile(functionalPath,refFileName);

% Apply the ANTS transform
applyANTsWarpToData(inFile,outFile,warpFile,refFile);

% Apply it to V1 Mask
V1InFile = fullfile(retinoPath,maskSaveName);
V1OutFile = fullfile(retinoPath,'V1_brainmask_MNI_EPI.nii.gz');
applyANTsWarpToData(V1InFile,V1OutFile,warpFile,refFile);

 %% Construct the model object
tfeHandle = tfeRPRF('verbosity','none');

%% Temporal definition of the stimulus and response
deltaT = 800; % in msecs
totalTime = 420*deltaT; % in msecs

% Spatial definition of the stimulus; we'll replace this later
load('/Users/radhika/Documents/flywheel/retinotopyTOMEAnalysis/flashStimForRPRF.mat');
stimulus.timebase = stimulus.timebase.*1000;


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

% File to write median volume to
writeFileNameX = 'RETINO_average_map_x.nii.gz';
writeFileNameY = 'RETINO_average_map_y.nii.gz';
writeFileNameSigma = 'RETINO_average_map_sigma.nii.gz';
writeFileNameFval = 'RETINO_average_map_fval.nii.gz';
writeFileNameA = 'RETINO_average_map_amplitude.nii.gz';

writeFilePathX = fullfile(functionalPath,writeFileNameX);
writeFilePathY = fullfile(functionalPath,writeFileNameY);
writeFilePathSigma = fullfile(functionalPath,writeFileNameSigma);
writeFilePathFval = fullfile(functionalPath,writeFileNameFval);
writeFilePathA = fullfile(functionalPath,writeFileNameA);

% Load GM Brainmask
niiBrainmask = MRIread(outFile);
niiBrainmaskVol = niiBrainmask.vol;

% Load V1 Brainmask
niiBrainmaskV1 = MRIread(V1OutFile);
niiBrainmaskV1Vol = niiBrainmaskV1.vol;

% Load temporary initial Functional Run
filePath = fullfile(functionalPath, functionalRuns{1});
    
% Load the file 
nii = MRIread(filePath);
fmri = nii.vol;

% Loop through desired files
for ff = 1:length(functionalRuns)
    filePath = fullfile(functionalPath, functionalRuns{ff});
    
    % Load the file 
    nii = MRIread(filePath);
    fmri = nii.vol;
    
    if ff==1
        avgVol = fmri;
        subplot(1,4,ff)
%         plot(squeeze(avgVol(42,25,41,:)));
    else
        avgVol = avgVol+fmri;
        subplot(1,4,ff)
%         plot(squeeze(avgVol(42,25,41,:)));
    end
end
avgVol = avgVol ./ length(functionalRuns);

% Identify which parts of the GM brainmask are the "brain"
indices = find(niiBrainmaskVol > 0.1 & niiBrainmaskV1Vol > 0.1);
[x,y,z] = ind2sub(size(niiBrainmaskVol),indices); 

% Create empty volume to fill
newVolX = nan(size(niiBrainmaskVol));
newVolY = nan(size(niiBrainmaskVol));
newVolSigma = nan(size(niiBrainmaskVol));
newVolFval = nan(size(niiBrainmaskVol));
newVolA = nan(size(niiBrainmaskVol));

% Time to fill it in with a median map
tic
for idx = 1:length(indices)
    % This is the step where you grab the time series for a voxel
    
    if idx==1
        avgTS = fmri(x(idx),y(idx),z(idx),:);
    else
        avgTS = avgTS+fmri(x(idx),y(idx),z(idx),:);
    end
    
    ts = squeeze(avgVol(x(idx),y(idx),z(idx),:))';
    ts = ts ./ mean(ts);
    ts = ts - mean(ts);
    thePacket.response.values = ts;

    % This is the call to the fitting engine
    [paramsFit,fVal,modelResponseStruct] = ...
    tfeHandle.fitResponse(thePacket,...
    'defaultParamsInfo', defaultParamsInfo,'errorType','1-r2');

    % Store the parameter value
%     newVolX(x(idx),y(idx),z(idx))= paramsFit.paramMainMatrix(1);
%     newVolY(x(idx),y(idx),z(idx))= paramsFit.paramMainMatrix(2);
%     newVolSigma(x(idx),y(idx),z(idx))= paramsFit.paramMainMatrix(3);
    newVolFval(x(idx),y(idx),z(idx))= fVal;
    newVolA(x(idx),y(idx),z(idx))= paramsFit.paramMainMatrix(4);
end
toc

avgTS = avgTS ./ length(indices);
plot(avgTS);
avgTS = squeeze(avgTS);    

% Write the new maps to the writeFiles
% donorX.vol = newVolX;
% donorX.fspec = writeFilePathX;
% MRIwrite(donorX,writeFilePathX);
% 
% donorY.vol = newVolY;
% donorY.fspec = writeFilePathY;
% MRIwrite(donorY,writeFilePathY);
% 
% donorSigma.vol = newVolSigma;
% donorSigma.fspec = writeFilePathSigma;
% MRIwrite(donorSigma,writeFilePathSigma);
% 
donorFval.vol = newVolFval;
donorFval.fspec = writeFilePathFval;
MRIwrite(donorFval,writeFilePathFval);

donorA.vol = newVolA;
donorA.fspec = writeFilePathA;
MRIwrite(donorA,writeFilePathA);