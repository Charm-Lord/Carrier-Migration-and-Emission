% test func_CME
clear all
clc

Data = csvread('10_0 PDOS_E0.csv');
Ep = Data(:,5);
dos = Data(:,6);

% 插值
Ep_1 = min(Ep):0.1:max(Ep)+3.1;
nE = length(Ep_1);
dos_1 = interp1(Ep,dos,Ep_1,'spline');
dE = Ep_1(2)-Ep_1(1);
Ep_1 = Ep_1';
dos_1 = dos_1';

% 费米分布
kT = 0.026;
F = @(x)1./(exp((x)./kT)+1);
F1 = F(Ep_1);

fi = dos_1.*F1;
fi0 = dos_1.*F1;

t = -20:.2:20;

tau = 3;
T = 2;
w = 2*pi/T;
fai = 0;

Tau_carrier = 10000;
npulse = 10;
pulse = 12500;

nn = -5:0.1:-2;
II = 10.^nn;
Wk = zeros(length(II),1);

parfor kk = 1:length(II)
    kk/length(II)
    I01 = II(kk);
    I02 = 0;
    t0 = 10;
    E1 = sqrt(I01).*exp(-2*log(2)*t.^2./tau.^2).*cos(w.*t+fai);
    E2 = sqrt(I02).*exp(-2*log(2)*(t-t0).^2./tau.^2).*cos(2*w.*(t-t0)+2*fai);
    fi1 = fi;
    for nn = 1:npulse
        [wk,fm1] = func_CME(Ep_1,dos_1,fi1,fi0,t,E1,E2,Tau_carrier);
        delt_fm = fm1-fi0;
        delt_fm = exp(-pulse/Tau_carrier)*delt_fm;
        fi1 = fi0+delt_fm;
    end
    Wk(kk) = wk;
end

figure()
loglog(II,Wk)
bb = gradient(log(Wk))./gradient(log(II'));
figure()
plot(bb)

