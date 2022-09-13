% Extracts values from an image data file based on regions/clusters of
% voxels in mask marsbar or image file(s) or voxel.
%
% Usage: [Ym R info] = extract_voxel_values(mask,data,sp,t)
%
% Ym    - matrix of output mean of values across voxels per image (row) per
%         region (column).
% R     - structure containing voxel info per region (R), per image (I):
%         R.I.Ya  - output vector of values of all voxels.
%         R.I.xyz - matrix indices of voxels in mask region/cluster.
%         R.I.mni - mni coordinates of voxels in mask region/cluster.
%         E.g. to get vector of value for region r, image i, after 
%         extraction, type R(r).I(i).Ya
% info  - structure containing region and image files used:
%         info.regions - string array of region files used.
%         info.images - string array of image files used.
% mask  - string array of mask image file name(s).
% data  - string array of image data file name(s). Can be any analyze or
%         nifti image that contains voxel values of interest, which can be
%         beta values, contrast values, t values etc.
%         e.g. ['beta_0005.img';'beta_0001.img'];
% sp    - space specification - Currently only 'mni' or 'subj'.
% t     - masking threshold if nii or img used during subj space mode.
%
% If input not complete, GUI invoked.
%
% Example usage(s):
% [Ym R info] = extract_voxel_values('L_Fusiform_roi.mat','beta_0005.img','mni');
%
% [Ym R info] = extract_voxel_values;
%
% Dependencies: SPM; Marsbar (if *_roi.mat ROI mask files used).
% 
% Josh Goh 24 May 2013 - Modified from extract_voxel_values_old.m.
% Josh Goh 28 May 2014 - Added facility for subject space.
% Josh Goh 05 Jun 2014 - Resolved bug for creating temp mask file.
% Josh Goh 19 Aug 2014 - Resolved bug for deleting temp mask file + NaN in mask files (for
%                        when there are -ve values in images).
% Josh Goh 25 Apr 2015 - Added masking threshold for nii/img during subj space mode.

function [Ym R info] = extract_voxel_values(mask,data,sp,t)

% Check input
if nargin<3
    mask = spm_select(Inf,'any','Select mask ROI files',[],pwd);
    data = spm_select(Inf,'image','Select data file (*.img or *.nii)',[],pwd);
    sp   = questdlg('Specify brain image space','Image Space','mni','subj','mni');
    if strcmp(sp,'subj')
    	t = str2num(cell2mat(inputdlg('Specify masking threshold as proportion','Masking Threshold',1,{'0.7'})));
    end
end

info.regions = mask;
info.images  = data;
currdir      = pwd;

% Loop image
for i = 1:size(data,1)
    
    % Read image header
    V = spm_vol(deblank(data(i,:)));
    
    % Loop regions
    for r = 1:size(mask,1)
        
        % Get voxels from mask
        [~,~,e] = fileparts(deblank(mask(r,:)));

		% Set space
		switch sp
			case {'mni'}
				switch e
					case {'.mat'}
						roi = maroi(deblank(mask(r,:)));
						xyz = voxpts(roi,deblank(data(i,:)));
					case {'.nii','.img'}
						roi = maroi_image(deblank(mask(r,:)));
						xyz = voxpts(roi,deblank(data(i,:)));
						%maskdata = spm_read_vols(spm_vol(deblank(mask(r,:))));
						%[x,y,z]  = ind2sub(size(maskdata),find(maskdata));
						%xyz      = [x y z]';
				end
		
				% Extract data
				Ya = spm_sample_vol(V,xyz(1,:),xyz(2,:),xyz(3,:),0);

				% Compute MNI coordinates      
				R(r).I(i).mni = vox2mni(V.mat,xyz);

			case {'subj'}
				switch e   
					case {'.mat'}
						error('Only .nii or .img masks allowed for subj space for now');
					case {'.nii','.img'}
						[p n ex]  = fileparts(V.fname);
						cd(p);
						%spm_mask(deblank(mask(r,:)),deblank(data(i,:)));
						spm_mask(deblank(mask(r,:)),deblank(data(i,:)),t);
						Vmaskfile = [pwd '/m' n ex];
						Vmaskdata = spm_vol(Vmaskfile);
						maskdata  = spm_read_vols(Vmaskdata);
						if length(size(maskdata))>3;
							maskdata = maskdata(:,:,:,size(maskdata,4));
						end
						% Check zero or NaN mask (when there are -ve responses)
						if min(min(min(maskdata))) == 0;	
							I = find(maskdata>min(min(min(maskdata))));
						else
							I = find(~isnan(maskdata(:)));
						end
						[x,y,z]   = ind2sub(size(maskdata),I);
						xyz       = [x y z]';
						delete(Vmaskfile);
						if ~isempty(strfind(ex,'.img'));
							Vmaskfilehdr = [pwd '/m' n '.hdr'];
							delete(Vmaskfilehdr);
						end
						cd(currdir);
				end

				% Extract data
				Ya = maskdata(I)';

				% No MNI coordinates
				R(r).I(i).mni = 'No mni for subj space';

		end
		
		% Log voxel data
		R(r).I(i).Ya = Ya;

		% Compute mean across voxels
		Ym(i,r) = mean(Ya(~isnan(Ya)));
		
		% Compute system coordinates       
		R(r).I(i).xyz = xyz;
        
    end
end
