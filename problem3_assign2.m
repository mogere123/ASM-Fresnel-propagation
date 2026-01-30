clear; clc; close all;
% Author: Mogere Mogere
% --------------------  parameters --------------------
lambda = 632e-9;          
D = 5e-3;                
R = D/2;                 
L = 10e-3;               
z_mm  = [10 20 50 100];  
z_list = z_mm * 1e-3;     

% -------------------- Spatial grid ------------------------
N  = 1024;                
dx = L/N;
x  = -L/2 : dx : L/2-dx;
y  = x;
[X,Y] = meshgrid(x,y);

% -------------------- Circular aperture ------------------
psi0 = double(sqrt(X.^2 + Y.^2) <= R);

% -------------------- Frequency grid ---------------------
fx = -1/(2*dx) : 1/L : 1/(2*dx)-(1/L);    
fy = fx;
[Fx,Fy] = meshgrid(fx,fy);

ASM = cell(1,4);
TF  = cell(1,4);
IR  = cell(1,4);

for iz = 1:4
    z = z_list(iz);

    psi_ASM = ASM_propagation(psi0, z, Fx, Fy, lambda);
    ASM{iz} = abs(psi_ASM);
  
    psi_TF = Fres_TF(psi0, z, Fx, Fy, lambda);
    TF{iz} = abs(psi_TF);
  

    psi_IR = Fres_IR(psi0, z, lambda, X, Y);
    IR{iz} = abs(psi_IR);
end


% -------------- figure --------------------------
figAll = figure('Color','w','Position',[50 50 1800 1100]);
tiledlayout(3,4,'Padding','loose','TileSpacing','loose');
sgtitle('Circular aperture propagation','FontSize',18,'FontWeight','bold');

% -------- ASM --------
for iz = 1:4
    nexttile
    imagesc(x*1e3,y*1e3,ASM{iz});
    axis image; axis xy; colormap pink
    title(sprintf('ASM  z = %d mm',z_mm(iz)),'FontSize',14)
    xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13)
end

% -------- Fresnel TF --------
for iz = 1:4
    nexttile
    imagesc(x*1e3,y*1e3,TF{iz});
    axis image; axis xy; colormap pink
    title(sprintf('Fres-TF  z = %d mm',z_mm(iz)),'FontSize',14)
    xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13)
end

% -------- Fresnel IR --------
for iz = 1:4
    nexttile
    imagesc(x*1e3,y*1e3,IR{iz});
    axis image; axis xy;colormap pink
    title(sprintf('Fres-IR  z = %d mm',z_mm(iz)),'FontSize',14)
    xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13)
end

% ---- colorbar ----
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Intensity';
cb.FontSize = 14;

drawnow
% ---------------saving----------------------------
exportgraphics(figAll,'figures.png','Resolution',300);
