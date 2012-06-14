function img2mask(Mimg)

% Load image first
V = spm_vol(Mimg);
Y = spm_read_vols(V);

% First round the image...
Y = round(Y);

% Anything above 0 is 1
Y(Y>0) = 1;
% Anything below 0 is 0
Y(Y<0) = 0;

% Any NaNs are 0
Y(isnan(Y)) = 0;

% Adjust V properties to ensure mask works fine...
V.dt(1) = 2;
V.pinfo(1) = 0; 

% Write image back...
spm_write_vol(V,Y);