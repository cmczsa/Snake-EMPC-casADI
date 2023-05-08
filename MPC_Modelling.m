function [XNext]=MPC_Modelling(X,U)
global USR;
X=mat2cell(X,[USR.N-1 1 1 1 USR.N-1 1 1 1]);
PHI=X{1};%N-1
theta=X{2};%1
Px=X{3};%1
Py=X{4};%1
dPHI=X{5};%N-1 Vphi
dtheta=X{6};%1 Vtheta
Vt=X{7};%1
Vn=X{8};%1

%先算下个时刻的P坐标
dPx=Vt*cos(theta)-Vn*sin(theta);
dPy=Vt*sin(theta)+Vn*cos(theta);
Px=Px+dPx*USR.tStep;
Py=Py+dPy*USR.tStep;
%再用当前时刻PHI&dtheta算下个时刻其他参数
dVt=-USR.ct/USR.m*Vt+(2*USR.cp/(USR.N*USR.m))*Vn*USR.e_'*PHI-(USR.cp/(USR.N*USR.m))*...
    PHI'*USR.A*USR.D_*dPHI;
dVn=-USR.cn/USR.m*Vn+2*USR.cp/(USR.N*USR.m)*Vt*USR.e_'*PHI;
%之后用这个时刻dPHI算下个时刻PHI，用这个时刻dtheta算下个时刻theta
PHI=PHI+dPHI*USR.tStep;
theta=theta+dtheta*USR.tStep;
%最后更新新的dPHI和dtheta
dPHI=dPHI+U*USR.tStep;
ddtheta=-USR.lamda1*dtheta+USR.lamda2/(USR.N-1)*Vt*USR.e_'*PHI;
dtheta=dtheta+ddtheta*USR.tStep;

Vt=Vt+dVt*USR.tStep;
Vn=Vn+dVn*USR.tStep;
VPHI=dPHI;
Vtheta=dtheta;

XNext=[PHI;theta;Px;Py;VPHI;Vtheta;Vt;Vn];
end



