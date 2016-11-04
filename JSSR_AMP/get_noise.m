function w = get_noise(m,sigma_Z)
%generates the noise signal

w = sigma_Z*randn(m,1);