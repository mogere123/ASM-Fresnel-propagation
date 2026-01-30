%Author: Mogere Mogere
function [psi_out] = Fres_IR(psi_in,z,lambda,X,Y)
% Takes input field and propagates using the Fred_IR propagator 
% to the output plane at distance 'z'
% Transfer Function for Propagation
k = 2*pi/lambda;
h = (exp(1i*k*z) ./ (1i*lambda*z)) .* exp(1i*k*(X.^2 + Y.^2) ./ (2*z)); % Impulse Response(IR)
H= fftshift(fft2(ifftshift(h)));                                        %FT                                    
FresIR_in = fftshift(fft2(ifftshift(psi_in)));

FresTF_out = FresIR_in.*H;                                              % multply in fourier domain
psi_out = fftshift(ifft2(ifftshift(FresTF_out)));                       % IFT of the output to spatial domain

end