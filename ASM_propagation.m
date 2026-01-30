% Author: Surya Kamal

function [psi_out] = ASM_propagation(psi_in,z,Fx,Fy,lambda,X,Y)
% Takes input field and propagates using the ASM propagator 
% to the output plane at distance 'z
% Transfer Function for Propagation

k = 2*pi/lambda;

H = exp(1i* k * z * sqrt(1 - (lambda*Fx).^2 - (lambda*Fy).^2 ) )...;
    .*(sqrt((lambda*Fx).^2 + (lambda*Fy).^2)<(1/lambda));

angular_spectrum_phi_in = fftshift(fft2(ifftshift(psi_in)));

angular_spectrum_phi_out = angular_spectrum_phi_in.*H;

psi_out = fftshift(ifft2(ifftshift(angular_spectrum_phi_out)));

end