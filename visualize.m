close all
clear all

%load files
load calcium
load currents
load dyad
load cacyt
load cajsr

Nx = 100;
Ny = 1;
Nz = 1;

figure
subplot(5,2,1);
plot(calcium(:,1),calcium(:,3),'LineWidth',2);      %cadyad
title('Cadyad (uM)');
subplot(5,2,3);
plot(calcium(:,1),calcium(:,4),'LineWidth',2)       %casub
title('Casub (mM)');
subplot(5,2,5);
plot(calcium(:,1),calcium(:,5),'LineWidth',2)       %cai
title('Cai (uM)');
subplot(5,2,7);
plot(calcium(:,1),calcium(:,6),'LineWidth',2);      %nsr
title('NSR (mM)');
subplot(5,2,9);
plot(calcium(:,1),calcium(:,7),'LineWidth',2);      %jsr
title('JSR (mM)');

subplot(5,2,2);
plot(dyad(:,1),dyad(:,2),'LineWidth',2);            %cads
title('Ca_{ds} (uM)');
subplot(5,2,4);
plot(dyad(:,1),dyad(:,3),'LineWidth',2)             %cass
title('Ca_{ss} (mM)');
subplot(5,2,6);
plot(dyad(:,1),dyad(:,4),'LineWidth',2)             %cacyt
title('Ca_{cyt} (uM)');
subplot(5,2,8);
plot(dyad(:,1),dyad(:,5),'LineWidth',2);            %cansr
title('Ca_{NSR} (mM)');
subplot(5,2,10);
plot(dyad(:,1),dyad(:,6),'LineWidth',2);            %cajsr
title('Ca_{JSR} (mM)');



figure
subplot(5,2,1);
plot(calcium(:,1), calcium(:,2),'LineWidth',2);       %v
title('Vm (mV)');
subplot(5,2,2);
plot(currents(:,1),currents(:,2),'LineWidth',2);      %ina
title('ina (uA/uF)');
subplot(5,2,3);
plot(currents(:,1),currents(:,3),'LineWidth',2)       %ilca
title('ilca (uA/uF)');
subplot(5,2,4);
plot(currents(:,1),currents(:,4),'LineWidth',2)       %Ito
title('Ito (uA/uF)');
subplot(5,2,5);
plot(currents(:,1),currents(:,5),'LineWidth',2);      %IKs
title('IKs (uA/uF)');
subplot(5,2,6);
plot(currents(:,1),currents(:,6),'LineWidth',2);      %IKr
title('IKr (uA/uF)');
subplot(5,2,7);
plot(currents(:,1),currents(:,7),'LineWidth',2);      %IK1
title('IK1 (uA/uF)');
subplot(5,2,8);
plot(currents(:,1),currents(:,8),'LineWidth',2);      %INaK
title('INaK (uA/uF)');
subplot(5,2,9);
plot(currents(:,1),currents(:,9),'LineWidth',2);      %inaca + inacass
title('inaca (uA/uF)');
subplot(5,2,10);
plot(currents(:,1),currents(:,10),'LineWidth',2);     %jrel
title('jrel (uM/ms)');


ca_cyt = zeros(size(cacyt(:,1)),Nz,Ny,Nx);
ca_jsr = zeros(size(cajsr(:,1)),Nz,Ny,Nx);

for t = 1:length(cacyt(:,1))
    for iz = 1:Nz
        for iy=1:Ny
            for ix = 1:Nx
                ca_cyt(t,iz,iy,ix) = cacyt(t, (iz-1)*Nx*Ny + (iy-1)*Nx + ix);
                ca_jsr(t,iz,iy,ix) = cajsr(t, (iz-1)*Nx*Ny + (iy-1)*Nx + ix);
            end
        end
    end
end

%figure
for t =1:size(ca_cyt(:,1,1,1));
%     subplot(2,1,1);
%     imagesc(squeeze(ca_cyt(t,Nz/2,:,:))); axis off; axis image; colorbar; caxis([0 0.8e-3]);
%     title('Ca_{myo} (mM)');
%     subplot(2,1,2);
%     imagesc(squeeze(ca_jsr(t,Nz/2,:,:))); axis off; axis image; colorbar; caxis([0 2]);
%     title('Ca_{JSR} (mM)');
%     set(gcf,'color','w');
%     F(t) = getframe(gcf);
    cacyt_linescan(t,:) = squeeze(ca_cyt(t,1,1,:));
    cajsr_linescan(t,:) = squeeze(ca_jsr(t,1,1,:));
end

% movie2avi(F,'calcium.avi');

figure
subplot(2,1,1);
imagesc(cacyt_linescan');  colorbar;  axis off;
title('Ca_{myo} (mM)');
subplot(2,1,2);
imagesc(cajsr_linescan');  colorbar;  axis off;
title('Ca_{JSR} (mM)');





