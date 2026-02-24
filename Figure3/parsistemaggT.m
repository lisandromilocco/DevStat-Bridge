function dX = parsistemaggT(t, X, theta1, theta2, IntegrationStep, devNoise_x1, devNoise_x2, u)
    i = fix(t / IntegrationStep) + 1;
    g1 = X(1);
    g2 = X(2);
    
    dg1 = (2 + theta1(i)) / (1 + (g2 / 2)^2) - 0.4 * g1 + devNoise_x1(i);
    dg2 = (2 + theta2(i)) / (1 + (g1 / (3 + u(i)))^2) - 0.4 * g2 + devNoise_x2(i);
    
    dX = [dg1; dg2];
end


%%%%%%%% exspresiones analiticas
% clear
% syms g1 g2 theta1 theta2 u
% dg1=(2+theta1)/(1+(g2/2)^2)-0.4*g1;
% dg2=(2+theta2)/(1+(g1/(3+u))^2)-0.4*g2;
% 
% a11=diff(dg1,g1);a12=diff(dg1,g2);
% a21=diff(dg2,g1);a22=diff(dg2,g2);
% 
% b11=diff(dg1,theta1);b12=diff(dg1,theta2);b13=diff(dg1,u);
% b21=diff(dg2,theta1);b22=diff(dg2,theta2);b23=diff(dg2,u);
% 
% 
% Ar=[a11 a12 ;a21 a22 ];
% Br=[b11 b12 b13;b21 b22 b23 ];
% save ABr Ar Br