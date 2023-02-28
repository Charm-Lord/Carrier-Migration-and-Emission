% 载流子迁移发射函数
% 输入态密度 dos,初始电子分布 fi,平衡态下电子分布 fi0
% 输出发射电子 Wk,末态电子分布 fm
function [Wk,fm] = func_CME(Ep,dos,fi,fi0,t,E1,E2,Tau_carrier)
    nE = length(Ep);
    dE = Ep(2)-Ep(1);
    dt = t(2)-t(1);
    WF = 4.4;
    Evm = -12;
    hw = 1.55;
    ff = zeros(nE,length(t));
    WW = zeros(length(t),1);
    ff(:,1) = fi;
    ki0 = dos-fi;
    kki = ki0/(sum(ki0)*dE);
    W = 0;

    for tt = 1:length(t)
        Qi = zeros(nE,1);
        fi_temp = ff(:,tt);
        E11 = E1(tt);
        E22 = E2(tt);
        E12 = E11+E22;
        E12 = abs(E12);
        P1 = E11.^2*dt*10;
        P2 = E22.^2*dt*10;
        P12 = E12.^2*dt*10;
        P12(P12>1) = 1;
        for ii = 1:nE
            E0 = Ep(ii);
            nP = P1+P2;
            P1P2 = P1/nP.*normpdf(Ep,E0+hw,0.06)+P2/nP.*normpdf(Ep,E0+2*hw,0.06);
            Hr = kki.*P1P2;
            if sum(Hr) ~= 0
                Hr = Hr/sum(Hr);
            end
            gi = fi_temp(ii).*dE.*P12*Hr;
            pi = gi;
            gi(ii) = fi_temp(ii)-sum(pi)+pi(ii);
            Qi = gi+Qi;
        end
    
        wk = exp(-0.1.*(WF-Ep)./E12);
        wk(wk > 1) = 1;
        wk(Ep > WF) = 1;
        W = W+sum(Qi.*wk)*dE;
        WW(tt) = W;
        Qi = Qi.*(1-wk);
        
        % 考虑载流子的衰减
        delt_pi = Qi-fi0;
        delt_pi = exp(-dt/Tau_carrier)*delt_pi;
        Qi = fi0+delt_pi;
    
        % 最低价带处电子始终充足
        Qi(Ep < Evm) = fi0(Ep < Evm);
    
        Qi(Qi > dos) = dos(Qi > dos);
        Qi(Qi < 0) = 0;
    
        ff(:,tt+1) = Qi;
        ki_temp = abs(dos-Qi);
        kki = ki_temp/(sum(ki_temp)*dE);
    end
    
    fm = ff(:,length(t));
    Wk = W;
end