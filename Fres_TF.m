%Author : Mogere Mogere

function [psi_out] = Fres_TF(psi_in,z,Fx,Fy,lambda)
% Takes input field and propagates using the ASM propagator 
% to the output plane at distance 'z'
% Transfer Function for Propagation

k = 2*pi/lambda;

H = exp(1i*k*z)*exp(-1i*pi*lambda*z*(Fx.^2+Fy.^2)); % Transfer function

FresTF_in = fftshift(fft2(ifftshift(psi_in)));      % FT of the field

FresTF_out = FresTF_in.*H;

psi_out = fftshift(ifft2(FresTF_out));             % IFT of the output

end