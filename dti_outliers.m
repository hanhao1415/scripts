function [b_slice_init ] = dti_outliers( dti_path, mask_path, bvals_path )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%read in the bvals
fid = fopen(bvals_path); bvals = fscanf(fid,'%d');

%read in the dti sequence and mask
dti = MRIread(dti_path); mask = MRIread(mask_path); 

%get logical mask
mask_log = logical(mask.vol); 

%loop through every frame, calculating mean intensity and natural log of
%the ratio of intensity to b0 intensity
dims = size(dti.vol); b_init = zeros(dims(4),1); b_ln = b_init;
b_slice_init = zeros(dims(4),dims(3)); b_slice_ln = b_slice_init;
for f=1:dims(4)
    dim = dti.vol(:,:,:,f);
    b_init(f,1) = mean(dim(mask_log));
    b_ln(f,1) = log(b_init(f,1)/b_init(1,1));
    for t=1:dims(3)
        b_slice = dti.vol(:,:,t,f);
        b_slice_init(f,t) = mean(b_slice(mask_log(:,:,t)));
        b_slice_ln(f,t) = log(b_slice_init(f,t)/b_slice_init(f,1));
    end
    scatter(bvals,b_slice_ln(:,f))
end

b_slice_init
b_slice_ln