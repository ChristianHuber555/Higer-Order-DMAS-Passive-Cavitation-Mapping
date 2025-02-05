clear all;
addpath("data")
load("RData_Real_MI_3_fs_1.mat");
RF = RData.data;
fs = RData.fs;

param.xarray = (-31.5:1:31.5)*0.3e-3;
param.c = 1540;


N = size(RF,2);

%--  field of view for the reconstruction
xdim = 51;
zdim = 51;
x = linspace(-2,2,xdim)*1e-3;
z = linspace(30,55,zdim)*1e-3;
[X, Z] = meshgrid(x, z);

param.t = 0:1/fs:size(RF,1)*(1/fs)-1/fs;
Cav_Map = zeros(size(X));
tic
RF = sign(RF).*nthroot(abs(RF),3);

parfor ix = 1:xdim
    xp = x(ix);
    for iz = 1:zdim                 
        zp = z(iz);  
        Cav_Map(iz, ix) = compute_PAM(RF, param, xp ,zp);
    end
    fprintf("Column %d from %d is finished", ix,xdim);
    fprintf('\n')
end

toc

%%
Cav_map_norm = Cav_Map./max(Cav_Map(:));
Cav_map_norm = 10*log10(Cav_map_norm);

figure()
imagesc(x*1000,z*1000,Cav_map_norm, [-10 0])
colorbar
set(gca,'DataAspectRatio',[1 1 1])

function pam_value = compute_PAM(RF_rcb, param, xp , zp)
    dx = sqrt((param.xarray - xp).^2 + zp.^2);
    RF_b = zeros(size(RF_rcb));
    for iarray = 1:64
        t_new = param.t + dx(iarray)/param.c;
        RF_b(:,iarray) = interp1(param.t,RF_rcb(:,iarray),t_new,'linear',0);
    end
    e1 = sum(RF_b,2);
    e2 = (e1*e1-e2)/2;
    pam_value = sum(e1.*e1);
end


