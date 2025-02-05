clear all;
addpath("data")

param.xarray = (-31.5:1:31.5)*0.3e-3;
param.c = 1540;

xdim = 51;
zdim = 51;
x = linspace(-2,2,xdim)*1e-3;
z = linspace(30,55,zdim)*1e-3;
[X, Z] = meshgrid(x, z);

Cav_Map = zeros(5,size(X,1),size(X,2));

tic
fs_i = [1 2 3 5 10];
for i = 1:5
    load(strcat("RData_Real_MI_3_fs_",num2str(fs_i(i)),".mat"));
    RF = RData.data;
    RF = sign(RF).*nthroot(abs(RF),fs_i(i));
    fs = RData.fs;
    param.t = 0:1/fs:size(RF,1)*(1/fs)-1/fs;
    parfor ix = 1:xdim
        xp = x(ix);
        for iz = 1:zdim                 
            zp = z(iz);  
            if i == 1
                Cav_Map(i,iz, ix) = compute_PAM_DAS(RF, param, xp ,zp);
            elseif i == 2
                Cav_Map(i,iz, ix) = compute_PAM_DMAS(RF, param, xp ,zp);
            elseif i == 3 
                Cav_Map(i,iz, ix) = compute_PAM_DMAS3(RF, param, xp ,zp);
            elseif i == 4
                Cav_Map(i,iz, ix) = compute_PAM_DMAS5(RF, param, xp ,zp);
            elseif i == 5
                Cav_Map(i,iz, ix) = compute_PAM_DMAS10(RF, param, xp ,zp);
            end
        end
        fprintf("Column %d from %d is finished", ix,xdim);
        fprintf('\n')
    end
    Cav_Map(i,:,:) = Cav_Map(i,:,:)./max(Cav_Map(i,:,:),[],'all');
    Cav_Map(i,:,:) = 10*log10(Cav_Map(i,:,:));
    toc
end


%%

figure()
subplot(1,5,1)
imagesc(x*1000,z*1000,squeeze(Cav_Map(1,:,:)), [-10 0])
colorbar
set(gca,'DataAspectRatio',[1 1 1])
subplot(1,5,2)
imagesc(x*1000,z*1000,squeeze(Cav_Map(2,:,:)), [-10 0])
colorbar
set(gca,'DataAspectRatio',[1 1 1])
subplot(1,5,3)
imagesc(x*1000,z*1000,squeeze(Cav_Map(3,:,:)), [-10 0])
colorbar
set(gca,'DataAspectRatio',[1 1 1])
subplot(1,5,4)
imagesc(x*1000,z*1000,squeeze(Cav_Map(4,:,:)), [-10 0])
colorbar
set(gca,'DataAspectRatio',[1 1 1])
subplot(1,5,5)
imagesc(x*1000,z*1000,squeeze(Cav_Map(5,:,:)), [-10 0])
colorbar
set(gca,'DataAspectRatio',[1 1 1])

function pam_value = compute_PAM_DAS(RF_rcb, param, xp , zp)
    dx = sqrt((param.xarray - xp).^2 + zp.^2);
    RF_b = zeros(size(RF_rcb));
    for iarray = 1:64
        t_new = param.t + dx(iarray)/param.c;
        RF_b(:,iarray) = interp1(param.t,RF_rcb(:,iarray),t_new,'linear',0);
    end
    e1 = sum(RF_b,2);
    pam_value = sum(e1.*e1);
end

function pam_value = compute_PAM_DMAS(RF_rcb, param, xp , zp)
    dx = sqrt((param.xarray - xp).^2 + zp.^2);
    RF_b = zeros(size(RF_rcb));
    for iarray = 1:64
        t_new = param.t + dx(iarray)/param.c;
        RF_b(:,iarray) = interp1(param.t,RF_rcb(:,iarray),t_new,'linear',0);
    end
    e1 = sum(RF_b,2);
    e2 = sum(RF_b.*RF_b,2);
    e1 = (e1.*e1-e2)/2;
    pam_value = sum(e1.*e1);
end

function pam_value = compute_PAM_DMAS3(RF_rcb, param, xp , zp)
    dx = sqrt((param.xarray - xp).^2 + zp.^2);
    RF_b = zeros(size(RF_rcb));
    for iarray = 1:64
        t_new = param.t + dx(iarray)/param.c;
        RF_b(:,iarray) = interp1(param.t,RF_rcb(:,iarray),t_new,'linear',0);
    end
    e1 = sum(RF_b,2);
    e2 = sum(RF_b.*RF_b,2);
    e3 = sum(RF_b.*RF_b.*RF_b,2);
    e1 = (e1.*e1.*e1- 3*e1.*e2 - 2*e3)/6; 
    pam_value = sum(e1.*e1);
end

function pam_value = compute_PAM_DMAS5(RF_rcb, param, xp , zp)
    dx = sqrt((param.xarray - xp).^2 + zp.^2);
    RF_b = zeros(size(RF_rcb));
    for iarray = 1:64
        t_new = param.t + dx(iarray)/param.c;
        RF_b(:,iarray) = interp1(param.t,RF_rcb(:,iarray),t_new,'linear',0);
    end
    e1 = sum(RF_b,2);
    e2 = sum(RF_b.*RF_b,2);
    e3 = sum(RF_b.*RF_b.*RF_b,2);
    e4 = sum(RF_b.*RF_b.*RF_b.*RF_b,2);
    e5 = sum(RF_b.*RF_b.*RF_b.*RF_b.*RF_b,2);
    e1 = (e1.*e1.*e1.*e1.*e1 - 10*e1.*e1.*e1.*e2 + 20*e3.*e1.*e1 + 15*e1.*e2.*e2 - 30*e4.*e1 - 20*e3.*e2 + 24.*e5)/120;
    pam_value = sum(e1.*e1);
end

function pam_value = compute_PAM_DMAS10(RF_rcb, param, xp , zp)
    dx = sqrt((param.xarray - xp).^2 + zp.^2);
    RF_b = zeros(size(RF_rcb));
    for iarray = 1:64
        t_new = param.t + dx(iarray)/param.c;
        RF_b(:,iarray) = interp1(param.t,RF_rcb(:,iarray),t_new,'linear',0);
    end
    e1 = sum(RF_b,2);
    e2 = sum(RF_b.*RF_b,2);
    e3 = sum(RF_b.*RF_b.*RF_b,2);
    e4 = sum(RF_b.*RF_b.*RF_b.*RF_b,2);
    e5 = sum(RF_b.*RF_b.*RF_b.*RF_b.*RF_b,2);
    e6 = sum(RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b,2);
    e7 = sum(RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b,2);
    e8 = sum(RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b,2);
    e9 = sum(RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b,2);
    e10 = sum(RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b.*RF_b,2);
    e1 = e1.*e1.*e1.*e1.*e1.*e1.*e1.*e1.*e1.*e1 - 45.*e1.*e1.*e1.*e1.*e1.*e1.*e1.*e1.*e2 +...
        240.*e1.*e1.*e1.*e1.*e1.*e1.*e1.*e3 + 630.*e1.*e1.*e1.*e1.*e1.*e1.*e2.*e2 - ...
        1260.*e1.*e1.*e1.*e1.*e1.*e1.*e4 - 5040.*e1.*e1.*e1.*e1.*e1.*e2.*e3 + ...
        6048.*e1.*e1.*e1.*e1.*e1.*e5 - 3150.*e1.*e1.*e1.*e1.*e2.*e2.*e2 + 18900.*e1.*e1.*e1.*e1.*e2.*e4 +...
        8400.*e1.*e1.*e1.*e1.*e3.*e3 - 25200.*e6.*e1.*e1.*e1.*e1 + 25200.*e1.*e1.*e1.*e2.*e2.*e3 - ...
        60480.*e1.*e1.*e1.*e2.*e5 - 50400.*e1.*e1.*e1.*e3.*e4 + 86400.*e7.*e1.*e1.*e1 + 4725.*e1.*e1.*e2.*e2.*e2.*e2 - ...
        56700.*e1.*e1.*e2.*e2.*e4 - 50400.*e1.*e1.*e2.*e3.*e3 + 151200.*e6.*e1.*e1.*e2 + 120960.*e1.*e1.*e3.*e5 + ...
        56700.*e1.*e1.*e4.*e4 - 226800.*e8.*e1.*e1 - 25200.*e1.*e2.*e2.*e2.*e3 + 90720.*e1.*e2.*e2.*e5 + ...
        151200.*e1.*e2.*e3.*e4 - 259200.*e7.*e1.*e2 + 22400.*e1.*e3.*e3.*e3 - 201600.*e6.*e1.*e3 - 181440.*e1.*e4.*e5 + ...
        403200.*e9.*e1 - 945.*e2.*e2.*e2.*e2.*e2 + 18900.*e2.*e2.*e2.*e4 + 25200.*e2.*e2.*e3.*e3 - ...
        75600.*e6.*e2.*e2 - 120960.*e2.*e3.*e5 - 56700.*e2.*e4.*e4 + 226800.*e8.*e2 - 50400.*e3.*e3.*e4 + ...
        172800.*e7.*e3 + 151200.*e6.*e4 + 72576.*e5.*e5 - 362880.*e10;
    e1 = e1./3628800;
    pam_value = sum(e1.*e1);
end


