function [DLI] = DynLaterIdx(data,atlas_path,l,r,w,s)
% Input:
%         data: a m * n array, m = number of rois, n = number of time points.
%         atlas_path: the path of the parcellation nifti file. 
%         l: index of roi of left brain.
%         r: index of roi of right brain.
%         w: window size.
%         s: step length.
% Depended on SPM package (https://www.fil.ion.ucl.ac.uk/spm/).
% Output:
%         Dynamic Laterality Index (DLI).
ws = fix((size(data,1) - w)/s) + 1;
DLI = zeros(ws,size(data,2));
[gs_L] = roi2gs(data(:,l,:),roi_size(l));
[gs_R] = roi2gs(data(:,r,:),roi_size(r));
for j = 1:ws
    time_idx = (((j - 1)*s + 1):((j - 1)*s + w));
    gL = gs_L(time_idx);gR = gs_R(time_idx);
    [LL] = corr(data(time_idx,:),gL);
    [RR] = corr(data(time_idx,:),gR);
    % -----------------------------------------%
    DLI(j,:) = (fisher_r2z(LL) - fisher_r2z(RR));  
end
disp('DLI finish.')
end

function [roi_size] = atlas_roi_size(atlas_path)
Va = spm_vol(atlas_path);
[Ya,~] = spm_read_vols(Va);
Atl = unique(Ya);Atl(1) = [];
roi_size = zeros(length(Atl),1);
for i = 1:length(Atl)
    roi_size(i,1) = sum(sum(sum(Ya == Atl(i))));
end
