function dX = sistemaggT(t,X)

global  theta1 theta2 IntegrationStep devNoise_x1 devNoise_x2 

i=fix(t/IntegrationStep)+1;
g1=X(1);g2=X(2);
dg1=(2+theta1(i))/(1+(g2/2)^2)-0.4*g1+devNoise_x1(i);
dg2=(2+theta2(i))/(1+(g1/(3))^2)-0.4*g2+devNoise_x2(i);
dX=[dg1;dg2];

end
