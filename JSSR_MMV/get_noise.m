function w = get_noise(m,sigma_Z,l)
%generates the noise signal

w = repmat(sigma_Z,m,1).*randn(m,l);