clear; clc; close all;
% Author: Mogere Mogere

% -------------------- parameters --------------------
lambda = 632e-9;          % [m]
D      = 5e-3;            % [m] aperture diameter
R      = D/2;             % [m]
L      = 6.25e-3;         % [m] window size
f      = 100e-3;          % [m] lens focal length
k      = 2*pi/lambda;

z_mm   = [10 20 50 100];  % [mm]
z_list = z_mm * 1e-3;      % [m]

% -------------------- Spatial grid ------------------------
N  = 1024;
dx = L/N;
x  = -L/2 : dx : L/2-dx;
y  = x;
[X,Y] = meshgrid(x,y);

% -------------------- Circular aperture + lens ------------
aperture   = double(sqrt(X.^2 + Y.^2) <= R);
lens_phase = exp(-1i*k*(X.^2 + Y.^2)/(2*f));      % lens phase
psi0       = aperture .* lens_phase;

% -------------------- Frequency ----------------------
fx = -1/(2*dx) : 1/L : 1/(2*dx)-(1/L);                           % cycles/m
fy = fx;
[Fx,Fy] = meshgrid(fx,fy);

% -------------------- Output folder -----------------------
outDir = "Problem4_figures";
if ~exist(outDir,'dir'); mkdir(outDir); end
%_--------------intensities-------------------------
ASM = cell(1,numel(z_list));
TF  = cell(1,numel(z_list));
IR  = cell(1,numel(z_list));

for iz = 1:numel(z_list)
    z = z_list(iz);

    psi_ASM = ASM_propagation(psi0, z, Fx, Fy, lambda);
    ASM{iz} = abs(psi_ASM);

    psi_TF = Fres_TF(psi0, z, Fx, Fy, lambda);
    TF{iz}=abs(psi_TF);

    psi_IR = Fres_IR(psi0, z, lambda, X, Y);
    IR{iz}= abs(psi_IR);
 end
  

%----------------FFT2 aperture ------------------------
Psi0_F = fftshift(fft2(ifftshift(psi0)));
FFT_intensity = abs(Psi0_F);
FFT_phase = angle(Psi0_F);

% ---------- Propagation at z = 100 mm ---------------
z100 = 100e-3;
psi_ASM_100 = ASM_propagation(psi0, z100, Fx, Fy, lambda);
psi_TF_100  = Fres_TF(psi0, z100, Fx, Fy, lambda);
psi_IR_100  = Fres_IR(psi0, z100, lambda, X, Y);

ASM_100 = abs(psi_ASM_100);
TF_100  = abs(psi_TF_100 ); 
IR_100  = abs(psi_IR_100 ); 

P_ASM_100 = angle(psi_ASM_100);
P_TF_100  = angle(psi_TF_100);
P_IR_100  = angle(psi_IR_100);


% FIGURE A-(12): (INTENSITy)(mm)

figA = figure('Color','w','Position',[50 50 2200 1400]);
tiledlayout(3,4,'Padding','loose','TileSpacing','loose');
sgtitle(sprintf('Intensity propagation | L=%.2f mm, D=%.2f mm, f=%.0f mm, N=%d', ...
    L*1e3, D*1e3, f*1e3, N), 'FontSize',18,'FontWeight','bold');

%---------- ASM-----------------
for iz = 1:4
    nexttile
    imagesc(x*1e3, y*1e3, ASM{iz}); axis image; axis xy;colormap pink
    title(sprintf('ASM, z=%d mm', z_mm(iz)), 'FontSize',14);
    xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13);
end

% ---------Fres-TF--------------------
for iz = 1:4
    nexttile
    imagesc(x*1e3, y*1e3, TF{iz}); axis image; axis xy;colormap pink
    title(sprintf('Fres TF, z=%d mm', z_mm(iz)), 'FontSize',14);
    xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13);
end

% ----------------Fres-IR----------------
for iz = 1:4
    nexttile
    imagesc(x*1e3, y*1e3, IR{iz}); axis image; axis xy; colormap pink
    title(sprintf('Fres IR, z=%d mm', z_mm(iz)), 'FontSize',14);
    xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13);
end

cbA = colorbar; cbA.Layout.Tile = 'east';
cbA.Label.String = 'intensity';
cbA.FontSize = 14;

exportgraphics(figA, fullfile(outDir,"intensity.png"), 'Resolution', 300);



% FIGURE B: FFT

figB = figure('Color','w','Position',[80 80 1800 650]);

subplot(1,2,1)
imagesc(fx/1e3, fy/1e3, FFT_intensity); axis image; axis xy; colormap pink; colorbar;
xlabel('(cycles/mm)'); ylabel('(cycles/mm)');
title('FFT2','FontSize',14);set(gca,'FontSize',13);
 
subplot(1,2,2)
imagesc(fx/1e3, fy/1e3, FFT_phase); axis image; axis xy;
colormap pink; colorbar;
xlabel('(cycles/mm)'); ylabel('(cycles/mm)');
title('phase angle','FontSize',14); set(gca,'FontSize',13);

sgtitle('Fourier transform of input field (aperture + lens)','FontSize',16,'FontWeight','bold');
exportgraphics(figB, fullfile(outDir,"B_FFT2_phase.png"), 'Resolution', 300);


% ---------FIGURE C: z=100 mm comparison — intensity + phase-------------

figC = figure('Color','w','Position',[100 60 2400 1100]);
tiledlayout(2,4,'Padding','loose','TileSpacing','loose');
% Row 1: Intensities
nexttile; imagesc(x*1e3, y*1e3,ASM_100); axis image; axis xy; colormap pink; colorbar;
title('ASM'); xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13);

nexttile; imagesc(x*1e3, y*1e3,TF_100); axis image; axis xy; colormap pink; colorbar;
title('Fres TF'); xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13);

nexttile; imagesc(x*1e3, y*1e3,IR_100); axis image; axis xy; colormap pink; colorbar;
title('Fres IR'); xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13);

nexttile; imagesc(fx/1e3,fy/1e3,FFT_intensity); axis image; axis xy; colormap pink; colorbar;
title('FFT'); xlabel('(cycles/mm)'); ylabel('(cycles/mm)'); set(gca,'FontSize',13);

% Row 2: Phases
nexttile; imagesc(x*1e3, y*1e3,P_ASM_100); axis image; axis xy; colormap pink; colorbar;
title('ASM phase'); xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13);

nexttile; imagesc(x*1e3, y*1e3,P_TF_100); axis image; axis xy; colormap pink; colorbar;
title('Fres TF phase'); xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13);

nexttile; imagesc(x*1e3, y*1e3,P_IR_100); axis image; axis xy; colormap pink; colorbar;
title('Fres IR phase'); xlabel('(mm)'); ylabel('(mm)'); set(gca,'FontSize',13);

nexttile; imagesc(fx/1e3,fy/1e3,FFT_phase); axis image; axis xy; colormap pink; colorbar;
title('FFT phase'); xlabel('(cycles/mm)'); ylabel('(cycles/mm)'); set(gca,'FontSize',13);

exportgraphics(figC, fullfile(outDir,"C_z100_compare_phase.png"), 'Resolution', 300);

disp("Folder: " + outDir);
