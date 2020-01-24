clc;
close all;
clear all;
Nd = 1000;                  % number of data symbols per frame
Ns = 1000;                  % Number of data frames
Nt = Nd*Ns;                 % total
nnumber of data symbols
EbNo = 0:2:20;              % SNR in dB
ModOrd = 1;                 % Modulation
Order (Number of bits per symbol)
j = sqrt(-1);                               % Imaginary Component
Map = modem.pskmod('M', 2^ModOrd, 'PhaseOffset', pi,'SymbolOrder','binary', 'InputType', 'bit');      % PSK Modulation syntax
Demap = modem.pskdemod(Map);            % PSKDemodulation syntax
t = 8;                      % No of Transmitting Users
r = 8;                      % No ofreceiving antennas
berzf = zeros(1,Ns);
bermm = zeros(1,Ns);
berzflas = zeros(1,Ns);
bermmlas = zeros(1,Ns);
Errzf = zeros(1,length(EbNo));
Errmm =zeros(1,length(EbNo));
Errzflas = zeros(1,length(EbNo));
Errmmlas = zeros(1,length(EbNo));
for idx = 1:length(EbNo)
    snr = EbNo(idx)+ 10*log10(ModOrd)
    sigmaS = t/(10^(0.1*snr));
    for i = 1:Ns
        s = randi([0 1],ModOrd*t,Nd);       % generation of data bits
        x = modulate(Map, s);                   % generation of data symbols
        H = (randn(r,t)+j*randn(r,t))/sqrt(2);
            % Rayleigh fat fading channel
        y = awgn(H*x,snr);
            % Received signal vector
        xzf = inv(H'*H)*H'*y;
            % Estimated signal vector using ZF receiver
        xmmse = inv(H'*H+sigmaS*eye(t))*H'*y;
            % Estimated signal vector using MMSE receiver
        Dzf = demodulate(Demap,xzf);
            % Estimated data bits
        Dmmse = demodulate(Demap,xmmse);
            % Estimated data bits
        xinzf = real(modulate(Map, Dzf));
            % generation of data symbols
        xinmm = real(modulate(Map, Dmmse));
            % generation of data symbols
        for m = 1:Nd
            Y = y(:,m);xz=xinzf(:,m);xm=xinmm(:,m);xz1=xz;xm1=xm;
            for n =1:t
                xz1(n)=-xz(n);xm1(n)=-xm(n);
                if ((norm(Y-H*xz1))^2)<((norm(Y-H*xz))^2)
                xz(n)=xz1(n);
                end
                if ((norm(Y-H*xm1))^2)<((norm(Y-H*xm))^2)
                xm(n)=xm1(n);
                end
            end
            xinzf(:,m)=xz;xinmm(:,m)=xm;
        end
        Dzflas = demodulate(Demap,xinzf);
        Dmmlas = demodulate(Demap,xinmm);
        berzf(i)= biterr(s, Dzf);
    % BER per data frame
        bermm(i)= biterr(s, Dmmse);
    % BER per data frame
        berzflas(i)= biterr(s, Dzflas);
    % BER per data frame
        bermmlas(i)= biterr(s, Dmmlas);
    % BER per data frame
    end
    Errzf(idx) = sum(berzf)/(Nt*t);
% Total BER per SNR
    Errmm(idx) = sum(bermm)/(Nt*t);
% Total BER per SNR
    Errzflas(idx) = sum(berzflas)/(Nt*t*idx);
% Total BER per SNR
    Errmmlas(idx) = sum(bermmlas)/(Nt*t*idx);
% Total BER per SNR
end

semilogy(EbNo,Errzf,'-ro','linewidth',2);
hold on;
semilogy(EbNo,Errmm,'-ko','linewidth',2);
hold on;
semilogy(EbNo,Errzflas,'-bo','linewidth',2);
hold on;
semilogy(EbNo,Errmmlas,'-go','linewidth',2);
grid on;