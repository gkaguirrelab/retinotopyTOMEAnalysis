% %% createSurfaceMap
% %
% % This script creates a median MNI time-series from multiple runs for a
% % session. It then converts this MNI volume into a cortical surface
% % representation.
% 
% %% Convenience variables
% projectName  = 'LFContrastAnalysis';
% flywheelName = 'LFContrast';
% subjID       = 'sub-HEROgka1';
% session      = 'ses-0411181853PM';
% sessionLabel = '04/11/18 18:53 PM';
% 
% %% Analysis labels that we are going to go and get
% fmriprepLabel   = 'fmriprep 04/12/2018 15:16:06';
% neuropythyLabel = 'retinotopy-templates 04/13/2018 16:46:22';
% fwInfo          = getAnalysisFromFlywheel(flywheelName,fmriprepLabel,fullfile(getpref('LFContrastAnalysis','projectRootDir'),'fmriprep'),'nodownload',true);
% sessionDir      = fullfile(getpref('LFContrastAnalysis','projectRootDir'),[fwInfo.subject,'_', fwInfo.timestamp(1:10)]);
% if ~isfolder(sessionDir)
%     mkdir(sessionDir);
% end
% if ~isfolder(fullfile(sessionDir,'fmriprep'))
%     fwInfo          = getAnalysisFromFlywheel(flywheelName,fmriprepLabel,fullfile(sessionDir,'fmriprep'));
%     fwInfoRetino    = getAnalysisFromFlywheel(flywheelName,neuropythyLabel,fullfile(sessionDir,'neuropythy'));
% end
% 
% 
% %% Relevant Nifti names for analysis
% 
% % functional runs
% functionalRuns = {'sub-HEROgka1_ses-0411181853PM_task-tfMRIFLASHAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
%     'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
%     'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastAP_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
%     'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
%     'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
%     'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'};
% 
% % functional masks
% functionalMasks = {'sub-HEROgka1_ses-0411181853PM_task-tfMRIFLASHAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz', ...
%     'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz', ...
%     'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastAP_run-2_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz', ...
%     'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz', ...
%     'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-2_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz', ...
%     'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-3_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz'};
% 
% % output files of Neuropythy (retinotopy template)
% retinoFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz',};
% 
% % brain mask of function run for the reference volume in ANTs step
% refFileName  =  'sub-HEROgka1_ses-0411181853PM_T1w_preproc.nii.gz';
% 
% % warp file name (product of running fmriprep)
% warpFileName = 'sub-HEROgka1_ses-0411181853PM_T1w_space-MNI152NLin2009cAsym_target-T1w_warp.h5';
% 
% % Set up paths to nifti and .h5 files
% retinoPath     = fullfile(sessionDir,'neuropythy');
% functionalPath = fullfile(sessionDir, 'fmriprep', subjID, session, 'func');
% warpFilePath   = fullfile(sessionDir, 'fmriprep', subjID, session, 'anat');
% 
% 
% %% Create restricted V1 mask
% 
% % load ecc nifti file
% eccenPos       = find(~cellfun(@isempty,strfind(retinoFiles,'eccen')));
% [~,tempName,~] = fileparts(retinoFiles{eccenPos});
% [~,outName,~]  = fileparts(tempName);
% eccenFileName  = fullfile(retinoPath,[outName '.nii.gz']);
% eccen          = MRIread(eccenFileName);
% 
% % load areas nifti file
% areasPos       = find(~cellfun(@isempty,strfind(retinoFiles,'areas')));
% [~,tempName,~] = fileparts(retinoFiles{areasPos});
% [~,outName,~]  = fileparts(tempName);
% areasFileName  = fullfile(retinoPath,[outName,'.nii.gz']);
% areas          = MRIread(areasFileName);
% 
% % make mask from the area and eccentricity maps
% areaNum     = 1;
% eccenRange  = [3 20];
% [~,maskSaveName] = makeMaskFromRetino(eccen,areas,areaNum,eccenRange,retinoPath);
% 
% %% Apply the warp to the mask and T1 files using ANTs
% 
% arrayForMed = [];
% files2warp = functionalRuns;
% for ii = 1:length(files2warp)
%     
%     % donor file
%     donor = MRIread(fullfile(functionalPath, functionalMasks{ii}));
%     
%     % input file
%     inFileTemp = fullfile(functionalPath,files2warp{ii});
%     inFile = fullfile(functionalPath,strcat(erase(files2warp{ii},'.nii.gz'),'_3D.nii.gz'));
%     nii = MRIread(inFileTemp);
%     timeSeries = nii.vol;
%     medianTimeSeries = median(timeSeries,4);
%     donor.vol = medianTimeSeries;
% %     donor.nframes = 1;
%     donor.fspec = inFile;
% %     donor.volres = nii.volres;
%     MRIwrite(donor,inFile);
%     % output file
%     [~,tempName,~] = fileparts(inFile);
%     [~,outName,~] = fileparts(tempName);
%     outFile = fullfile(retinoPath, strcat(erase(outName,'-MNI152NLin2009cAsym_preproc'),'_native.nii.gz'));
%     
%     % reference file
%     refFile = fullfile(warpFilePath,refFileName);
%     
%     % warp file
%     warpFile = fullfile(warpFilePath,warpFileName);
%     if ~exist(outFile)
%         applyANTsWarpToData(inFile, outFile, warpFile, refFile);
%     else
%         [~,fileName,~] = fileparts(outFile);
%         display(sprintf('%s already exist in the specified directory',fileName));
%     end
%     test = MRIread(outFile);
%     
%      arrayForMed = cat(4, arrayForMed, test.vol);
% end
% 
% %% Take a median of all the native functional median files
% 
% MedianVolAcrossRuns = median(arrayForMed,4);
% resultFile = fullfile(retinoPath,'averageForFunctionalRuns_median_native.nii.gz');
% donor1 = MRIread(outFile);
% donor1.vol = MedianVolAcrossRuns;
% donor1.fspec = resultFile;
% MRIwrite(donor1,resultFile);

%% Convert Volume Map to Cortical Surface Maps

% Set the SUBJECTS_DIR environment variable
setenv('SUBJECTS_DIR', fullfile(sessionDir,'freesurfer'));

% Create variable for freesurfer fsaverage directory
fsaveragePath = fullfile(sessionDir,'freesurfer','fsaverage');

% Convert fsaverage surface files into GIFTI surface files
% First the WM files
system(['mris_convert ' fullfile(fsaveragePath,'surf','lh.inflated') ' ' fullfile(fsaveragePath,'surf','WM_fsaverage.L.surf.gii')]);
system(['mris_convert ' fullfile(fsaveragePath,'surf','rh.inflated') ' ' fullfile(fsaveragePath,'surf','WM_fsaverage.R.surf.gii')]);

% Next the sphere files
system(['mris_convert ' fullfile(fsaveragePath,'surf','lh.sphere') ' ' fullfile(fsaveragePath,'surf','sphere_fsaverage.L.surf.gii')]);
system(['mris_convert ' fullfile(fsaveragePath,'surf','rh.sphere') ' ' fullfile(fsaveragePath,'surf','sphere_fsaverage.R.surf.gii')]);

% Convert fsaverage curvature files into GIFTI metric shape files
system(['mri_convert ' fullfile(fsaveragePath,'surf','lh.sphere') ' ' fullfile(fsaveragePath,'surf','curvature_fsaverage.L.shape.gii')]);
system(['mri_convert ' fullfile(fsaveragePath,'surf','rh.sphere') ' ' fullfile(fsaveragePath,'surf','curvature_fsaverage.R.shape.gii')]);

% Register NIFTI volume to GIFTI fsaverage files to produce fsaverage
% cortical maps
cmd = [fullfile(getpref('flywheelMRSupport','binWorkbench'),'wb_command') ' -volume-to-surface-mapping ' resultFile ' ' fullfile(fsaveragePath,'surf','WM_fsaverage.L.surf.gii') ' ' fullfile(fsaveragePath,'surf','MedianCorticalMap_fsaverage.L.shape.gii') ' -cubic'];
system(cmd);
cmd = [fullfile(getpref('flywheelMRSupport','binWorkbench'),'wb_command') ' -volume-to-surface-mapping ' resultFile ' ' fullfile(fsaveragePath,'surf','WM_fsaverage.R.surf.gii') ' ' fullfile(fsaveragePath,'surf','MedianCorticalMap_fsaverage.R.shape.gii') ' -cubic'];
system(cmd);

% Create variable for freesurfer subject native directory
nativePath = fullfile(sessionDir,'freesurfer',subjID);

% Now to convert freesurfer native files into GIFTI surface format
% First the WM files
system(['mris_convert ' fullfile(nativePath,'surf','lh.inflated') ' ' fullfile(nativePath,'surf','WM_native.L.surf.gii')]);
system(['mris_convert ' fullfile(nativePath,'surf','rh.inflated') ' ' fullfile(nativePath,'surf','WM_native.R.surf.gii')]);

% Next the sphere files
system(['mris_convert ' fullfile(nativePath,'surf','lh.sphere') ' ' fullfile(nativePath,'surf','sphere_native.L.surf.gii')]);
system(['mris_convert ' fullfile(nativePath,'surf','rh.sphere') ' ' fullfile(nativePath,'surf','sphere_native.R.surf.gii')]);

% Convert native curvature files into GIFTI metric shape files
system(['mri_convert ' fullfile(nativePath,'surf','lh.sphere') ' ' fullfile(nativePath,'surf','curvature_fsaverage.L.shape.gii')]);
system(['mri_convert ' fullfile(nativePath,'surf','rh.sphere') ' ' fullfile(nativePath,'surf','curvature_fsaverage.R.shape.gii')]);

% Convert fsaverage cortical maps to subject native cortical maps
system(['mri_surf2surf --srcsubject fsaverage --sval ' fullfile(fsaveragePath,'surf','MedianCorticalMap_fsaverage.L.shape.gii') ' --trgsubject ' subjID ' --tval ' fullfile(nativePath,'surf','MedianCorticalMap_native.L.shape.gii')]);
system(['mri_surf2surf --srcsubject fsaverage --sval ' fullfile(fsaveragePath,'surf','MedianCorticalMap_fsaverage.R.shape.gii') ' --trgsubject ' subjID ' --tval ' fullfile(nativePath,'surf','MedianCorticalMap_native.R.shape.gii')]);

%% Done