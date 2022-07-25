
H = [];
iter = 3125;

for i=1:iter
    channel=(randn(1,Ntap)+1i*randn(1,Ntap)).*sqrt(Power/2);	
    h22=zeros(1,Lch); 
    h22(Delay+1)=channel; % cir: channel impulse response
    h22r2=h22;
    H22r2=fft([h22r2 zeros(1,Nfft-Lch)]); % Channel frequency response
    H = [H H22r2];
end
H=H';
size(H)

chan = [real(H) imag(H)];
csvwrite('myFile.csv',chan);
