%Author: Mogere Mogere
clear;clc;

%------------Parameters------------------
lambda =632.8e-9;                   % [m]red light
k =2*pi/lambda;                    
D =2.5e-03;                         % [m]aperture diameter

%----------spatial domain--------------------
N =1024;                            % grid size
L =5e-03;                           % [m] window [m] 
dx =L/N;                            % sampling interval
x = -L/2 : dx :L/2-dx; 
y = x;
[X, Y] = meshgrid(x,y);

% -----------Fourier domain-----------------------
fx = -1/(2*dx) : 1/L : 1/(2*dx)-(1/L);    %frequency coordinates
fy = fx;
[Fx,Fy] =meshgrid(fx,fy);

%----------apertures---------------
input_field1 = sqrt(X.^2 +Y.^2) < D/2;    % circular aperture
input_field2 = squares_apeture(N);        % square aperture

%----------z list---------------
z_list = 0e-03:1e-03:100e-03;
Nz = numel(z_list);

% axes in mm 
xmm = x*1e3; 
ymm = y*1e3;

% ----------------figure-----------------
figure('Color','w');

ax1 = subplot(131); axis(ax1,'image'); colormap(ax1,'pink');
xlabel(ax1,'mm'); ylabel(ax1,'mm'); set(ax1,'FontSize',18);

ax2 = subplot(132); axis(ax2,'image'); colormap(ax2,'pink');
xlabel(ax2,'mm'); ylabel(ax2,'mm'); set(ax2,'FontSize',18);

ax3 = subplot(133); axis(ax3,'image'); colormap(ax3,'pink');
xlabel(ax3,'mm'); ylabel(ax3,'mm'); set(ax3,'FontSize',18);

% -------------------- circular-------------------------------
sgtitle('Propagation of |E(x,y)| (Circular aperture) using Fresnel TF, Fresnel IR and ASM');

% for iz = 1:Nz
%     z = z_list(iz);
% 
%     %-----------------Circular field-----------------------
%     psi_ir  = Fres_IR(input_field1, z,lambda,X,Y);
%     psi_tf  = Fres_TF(input_field1, z, Fx, Fy, lambda);
%     psi_asm = ASM_propagation(input_field1, z, Fx, Fy, lambda);
% 
%     %-----------------Display------------------
%     subplot(131)
%     imagesc(xmm,ymm,abs(psi_tf)); axis image;colormap pink;
%     title("Fresnel TF |E(x,y)| at z="+z*1e3+"mm",'FontSize',12);
% 
%     subplot(132)
%     imagesc(xmm,ymm,abs(psi_ir));axis image;colormap pink;
%     title("Fresnel IR |E(x,y)| at z="+z*1e3+"mm",'FontSize',12);
% 
%     subplot(133)
%     imagesc(xmm,ymm,abs(psi_asm));axis image;colormap pink;
%     title("ASM Propagation |E(x,y)| at z="+z*1e3+"mm",'FontSize',12);
% 
%     drawnow limitrate
% end

% ------------------ square-------------------------------
sgtitle('Propagation of |E(x,y)| (Square aperture) using Fresnel TF, Fresnel IR and ASM');

for iz = 1:Nz
    z = z_list(iz);

    %-----------------Square field-----------------------
    psi_ir  = Fres_IR(input_field2, z,lambda,X,Y);
    psi_tf  = Fres_TF(input_field2, z, Fx, Fy, lambda);
    psi_asm = ASM_propagation(input_field2, z, Fx, Fy, lambda);

    %-----------------Display inside loop------------------
    subplot(131)
    imagesc(xmm,ymm,abs(psi_tf));axis image
    title("Fresnel TF |E(x,y)| at z="+z*1e3+"mm",'FontSize',12);colormap pink;

    subplot(132)
    imagesc(xmm,ymm,abs(psi_ir));axis image;
    title("Fresnel IR |E(x,y)| at z="+z*1e3+"mm",'FontSize',12);colormap pink;

    subplot(133)
    imagesc(xmm,ymm,abs(psi_asm));axis image; colormap pink;
    title("ASM Propagation |E(x,y)| at z="+z*1e3+"mm",'FontSize',12);

    drawnow limitrate
    pause(0.1)
end
