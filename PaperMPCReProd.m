%% ----☆----☆----Paper Replication of EMPC (58th CDC)  ----☆----☆------%%

%You need to download casADi toolbox to run the code.
% Press the 'Run' button and wait for about 5 minutes.

clear;
clc

  %% ------☆-----☆-----☆----Calculate Basic Matrix----☆-----☆-----☆-------%%
global USR;
USR.N=9;

USR.A=zeros(USR.N-1,USR.N);
    for ii=1:USR.N-1
        USR.A(ii,ii)=1;
        USR.A(ii,ii+1)=1;
    end
    
USR.D=zeros(USR.N-1,USR.N);
    for ii=1:USR.N-1
        USR.D(ii,ii)=1;
        USR.D(ii,ii+1)=-1;
    end
    
USR.e=zeros(USR.N,1);
    for ii=1:USR.N
        USR.e(ii)=1;
    end
    
USR.e_=zeros(USR.N-1,1);
    for ii=1:USR.N-1
        USR.e_(ii)=1;
    end 
    
USR.D_=USR.D'*inv(USR.D*USR.D');

  %% --------☆--------☆--------Basic Vars Definition---------☆-------☆---------%%
USR.L=0.14;
USR.m=1;
USR.ct=1;
USR.cn=3;
USR.cp=(USR.cn-USR.ct)/(2*USR.L);
USR.lamda1=0.5;
USR.lamda2=20; 

%Controller Param
Kp=20;
Kd=5;

%Gait pattern params
alpha=0.05;
omega=120*(pi/180);
delta=40*(pi/180);

%Sample Params
USR.tStep=0.05;
tEnd=15;
tNum=tEnd/USR.tStep+1;
USR.t=0:USR.tStep:tEnd;

 %% ----☆-----☆----Lateral Undulation Definition & Initialization----☆----☆----%%
%Definition 
PHIRef=zeros(USR.N-1,tNum);
VPHIRef=zeros(USR.N-1,tNum);
URef=zeros(USR.N-1,tNum);

PHI=zeros(USR.N-1,tNum);
theta=zeros(1,tNum);
Px=zeros(1,tNum);
Py=zeros(1,tNum);
VPHI=zeros(USR.N-1,tNum);
Vtheta=zeros(1,tNum);
Vt=zeros(1,tNum);
Vn=zeros(1,tNum);
X_LU=zeros(2*USR.N+4,tNum);
U_LU=zeros(USR.N-1,tNum);

%Initialization
PHI(:,1)=[0;0.01;-0.01;0.01;0;0;0.01;-0.01];
X_LU(:,1)=[PHI(:,1);theta(1);Px(1);Py(1);VPHI(:,1);Vtheta(1);Vt(1);Vn(1)];

    %% ---------☆------☆-------LU Calculation-------☆------☆--------%%
for jj=1:tNum %时间
        for ii=1:USR.N-1 %关节
            PHIRef(ii,jj)=alpha*sin(omega*USR.t(jj)+(ii-1)*delta);
            VPHIRef(ii,jj)=alpha*omega*cos(omega*USR.t(jj)+(ii-1)*delta);%dPHIRef
            URef(ii,jj)=-alpha*omega*omega*sin(omega*USR.t(jj)+(ii-1)*delta);%ddPHIRef
        end
end

for m=1:tNum-1
    fprintf('LU计算  当前进展：第%d时刻 \n',m);
    U_LU(:,m)=URef(:,m)+Kd*(VPHIRef(:,m)-VPHI(:,m))+Kp*(PHIRef(:,m)-PHI(:,m));
     [PHI(:,m+1),theta(m+1),Px(m+1),Py(m+1),VPHI(:,m+1),Vtheta(m+1),...
     Vt(m+1),Vn(m+1),X_LU(:,m+1)]=Modelling(X_LU(:,m),U_LU(:,m));
end
fprintf('----☆----☆-----☆----☆----☆-----☆----☆----☆\n')

  %% ---------☆------☆-----MPC Configuration ------☆------☆--------%%
global State;
State.Np=20; %预测区间
State.N=USR.N;

PHIMax=0.052;
VPHIMax=0.109;
UMax=0.2276;

%casADi基本变量定义
U_EMPC = SX.sym('U',State.N-1,State.Np);
P = SX.sym('P',2*State.N+4);
X_EMPC = SX.sym('X',2*State.N+4,State.Np+1);
X_EMPC(:,1) = P;

for ii=1:State.Np
    X_EMPC( : ,ii+1)=MPC_Modelling(X_EMPC( : ,ii),U_EMPC(:,ii));
end
clc

%% -----☆------☆------MPC Calculation (1)-(gamma=0) -----☆-----☆-----%%
%Objective Function
ObjFun=-sum( X_EMPC(2*State.N+3,: ) );

%Nonlinear Constraints
NonLnrCons=[];
for ii = 1:State.Np+1
    NonLnrCons = [NonLnrCons ; X_EMPC(1:State.N-1,ii)];%PHI
    NonLnrCons = [NonLnrCons ; X_EMPC(State.N+3:2*State.N+1,ii)];%VPHI
end
%NLP Struct
OptimU=reshape( U_EMPC , (State.N-1)*State.Np,1);
NlpMPC=struct('f',ObjFun,'x',OptimU,'g',NonLnrCons ,'p',P);

%NLP Solver
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;%Acceptable Convergence Tolerance
%opts.ipopt.acceptable_obj_change_tol = 1e-6;
MPCSolver=nlpsol('MPCSolver', 'ipopt', NlpMPC,opts);

%Upper and lower bound constraints
args=struct;
args.lbg=[];
args.ubg=[];

for ii = 1:State.Np+1
    args.lbg = [args.lbg ; ones(State.N-1,1)* -PHIMax];
    args.lbg = [args.lbg ; ones(State.N-1,1)*  -VPHIMax];
    args.ubg = [args.ubg ; ones(State.N-1,1)* PHIMax];
    args.ubg = [args.ubg ; ones(State.N-1,1)*  VPHIMax];
end

args.lbx = -UMax*ones((State.N-1)*State.Np,1);
args.ubx = UMax*ones((State.N-1)*State.Np,1);

%Vars Initialization
PHIInit=[0,0.01,-0.01,0.01,0,0,0.01,-0.01]';
thetaInit=0;
PxInit=0;
PyInit=0;
VPHIInit=zeros(State.N-1,1);
VthetaInit=0;
VtInit=0;
VnInit=0;
XInit=[PHIInit ; thetaInit ; PxInit ; PyInit ; VPHIInit ; VthetaInit ; VtInit ; VnInit ];
UEveryOptim=zeros((State.N-1)*State.Np,1);%OptAnswer Per Time             每一次求解送入的U
ALLUMPC1=zeros(State.N-1,tNum-1);%OptAnswer-U                                   所有时刻的最优控制序列U
ALLXMPC1=zeros(2*State.N+4,tNum);%OptAnswer-X                                  所有时刻的最优状态序列X
ALLXMPC1( : ,1)=XInit;

%Start MPC  (1) Optim
StartOptimTime = tic;
for m=1:tNum-1
    fprintf('MPC(1) 当前进展：第%d时刻 \n',m);
    args.p=ALLXMPC1( : ,m);
%     args.x0 = reshape(UEveryOptim,(State.N-1)*State.Np,1);
    args.x0=UEveryOptim;
    OptimAns=MPCSolver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg,...
                        'ubg', args.ubg,'p',args.p);
    UEveryOptim=full(OptimAns.x);%将casADi的DM类型转化为double
    
     ALLUMPC1(:,m)=UEveryOptim(1:State.N-1);
     ALLXMPC1( : ,m+1)=MPC_Modelling(ALLXMPC1( : ,m),ALLUMPC1(:,m));
     UEveryOptim(1:(end-State.N+1))=UEveryOptim(State.N:end);
end
MPCOptimTime=toc(StartOptimTime);
fprintf('MPC(1)优化结束，用时%f 秒 \n',MPCOptimTime);
clc
fprintf('----☆----☆-----☆----☆----☆-----☆----☆----☆\n')

 %% -----☆-----☆-----MPC Calculation (2)-(gamma=0.025) ------☆-----☆----%%

 gamma=0.025;
%Objective Function
ObjFun=0;
for ii=1:State.Np
    TempU=U_EMPC(:,ii);
    ObjFun=ObjFun-X_EMPC(2*State.N+3,ii+1)+gamma*(TempU'*TempU );
end

%Nonlinear Constraints
NonLnrCons=[];
for ii = 1:State.Np+1
    NonLnrCons = [NonLnrCons ; X_EMPC(1:State.N-1,ii)];%state phi
    NonLnrCons = [NonLnrCons ; X_EMPC(State.N+3:2*State.N+1,ii)];%state Vphi
end

%NLP Struct
OptimU=reshape( U_EMPC , (State.N-1)*State.Np,1);
NlpMPC=struct('f',ObjFun,'x',OptimU,'g',NonLnrCons ,'p',P);

%NLP Solver
opts = struct;
opts.ipopt.max_iter = 150;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-9;%Acceptable Convergence Tolerance
%opts.ipopt.acceptable_obj_change_tol = 1e-6;
MPCSolver=nlpsol('MPCSolver', 'ipopt', NlpMPC,opts);

%Upper and lower bound constraints
args=struct;
args.lbg=[];
args.ubg=[];
for ii = 1:State.Np+1
    args.lbg = [args.lbg ; ones(State.N-1,1)* -PHIMax];
    args.lbg = [args.lbg ; ones(State.N-1,1)*  -VPHIMax];
    args.ubg = [args.ubg ; ones(State.N-1,1)* PHIMax];
    args.ubg = [args.ubg ; ones(State.N-1,1)*  VPHIMax];
end
args.lbx = -UMax*ones((State.N-1)*State.Np,1);
args.ubx = UMax*ones((State.N-1)*State.Np,1);

%Vars Initialization
PHIInit=[0,0.01,-0.01,0.01,0,0,0.01,-0.01]';
thetaInit=0;
PxInit=0;
PyInit=0;
VPHIInit=zeros(State.N-1,1);
VthetaInit=0;
VtInit=0;
VnInit=0;
XInit2=[PHIInit ; thetaInit ; PxInit ; PyInit ; VPHIInit ; VthetaInit ; VtInit ; VnInit ];
UEveryOptim2=zeros((State.N-1)*State.Np,1);%每一次求解送入的U
ALLUMPC2=zeros(State.N-1,tNum-1);%所有时刻的最优控制序列U
ALLXMPC2=zeros(2*State.N+4,tNum);%所有时刻的最优状态序列X
ALLXMPC2( : ,1)=XInit2;

%Start MPC  (2) Optim
StartOptimTime2 = tic;
for m=1:tNum-1
    fprintf('MPC(2) 当前进展：第%d时刻 \n',m);
    args.p=ALLXMPC2( : ,m);
%     args.x0 = reshape(UEveryOptim,(State.N-1)*State.Np,1);
    args.x0=UEveryOptim2;
    OptimAns=MPCSolver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg,...
                        'ubg', args.ubg,'p',args.p);
    UEveryOptim2=full(OptimAns.x);%将casADi的DM类型转化为double
    
     ALLUMPC2(:,m)=UEveryOptim2(1:State.N-1);
     ALLXMPC2( : ,m+1)=MPC_Modelling(ALLXMPC2( : ,m),ALLUMPC2(:,m));
     UEveryOptim2(1:(end-State.N+1))=UEveryOptim2(State.N:end);
end
MPCOptimTime2=toc(StartOptimTime2);
fprintf('MPC(2)优化结束，用时%f 秒 \n',MPCOptimTime2);
clc
fprintf('----☆----☆-----☆----☆----☆-----☆----☆----☆\n')

   %% ------☆------☆------☆-------PLOT (1) & (2) ------☆------☆------☆-----%%

PlotTime=0:tNum-1;

%Vt
figure(1)
hold on
LU=plot(PlotTime,X_LU(2*USR.N+3,:),'color',[0.3,0.6,1],'LineWidth',1);
EMPC1=plot(PlotTime,ALLXMPC1(2*State.N+3,:),'r','linewidth',1);
EMPC2=plot(PlotTime,ALLXMPC2(2*State.N+3,:),'color',[1,0.5,0],'linewidth',1);
gcaTemp=gca;
gcaTemp.YAxis.Exponent = -2;
box on
set(gcf,'Position',[494   349   600   370])
axis([0 PlotTime(end) 0 0.07])
set(gca,'ytick',[0:.02:.06])
xlabel('Time Steps (1/20 s)')
ylabel('Forward Velocity \it{v_t (m/s)}')
legend([LU EMPC1 EMPC2],'LU','EMPC,\gamma=0','EMPC,\gamma=0.025');

axes('Position',[0.55,0.2,0.3,0.28]);
hold on
plot(PlotTime(160:220),X_LU(2*USR.N+3,160:220),'color',[0.3,0.6,1],'linewidth',1);
plot(PlotTime(160:220),ALLXMPC1(2*State.N+3,160:220),'r','linewidth',1);
plot(PlotTime(160:220),ALLXMPC2(2*State.N+3,160:220),'color',[1,0.5,0],'linewidth',1);
ax=gca;
ax.YAxis.Exponent = -2;
set(ax,'xtick',[160:20:200])
box on

%VPHI3
figure(2)
hold on
LU=plot(PlotTime,X_LU(State.N+5,:),'color',[0.3,0.6,1],'linewidth',1);
MPC=plot(PlotTime,ALLXMPC1(State.N+5,:),'r','linewidth',1);
plot(PlotTime,VPHIMax*ones(size(PlotTime)),'--','color',[1,0.5,0],'linewidth',1.5);
plot(PlotTime,-VPHIMax*ones(size(PlotTime)),'--','color',[1,0.5,0],'linewidth',1.5);
legend([LU MPC],'LU','EMPC');
axis([0 PlotTime(end) -0.15 0.15])
xlabel('Time Steps \it{(1/20 s)}')
ylabel('Joint Velocity \it{v_{\phi,3}} (m/s)')
box on
set(gcf,'Position',[494   349   560   370])

axes('Position',[0.65,0.2,0.2,0.25]);
hold on
plot(PlotTime(150:180),X_LU(State.N+5,150:180),'color',[0.3,0.6,1],'linewidth',1);
plot(PlotTime(150:180),ALLXMPC1(State.N+5,150:180),'r','linewidth',1);
axis([150 180 0 0.13])
box on

   %% ------☆-------☆-------☆------- Cal Table 1  ------☆-------☆-------☆-----%%
E_LU=0;
E_MPC=0;
E_MPC2=0;

for t=100:300
    E_LU=E_LU+U_LU(:,t)'*U_LU(:,t);
    E_MPC=E_MPC+ALLUMPC1(:,t)'*ALLUMPC1(:,t);
    E_MPC2=E_MPC2+ALLUMPC2(:,t)'*ALLUMPC2(:,t);
end

E_LU=E_LU/201;
E_MPC=E_MPC/201;
E_MPC2=E_MPC2/201;
E_absolute=[E_LU;E_MPC;E_MPC2];
E_relative=[E_LU/E_LU;E_MPC/E_LU;E_MPC2/E_LU];

Vav_LU=sum(X_LU(2*USR.N+3,100:300))/201;
Vav_MPC=sum(ALLXMPC1(2*State.N+3,100:300))/201;
Vav_MPC2=sum(ALLXMPC2(2*State.N+3,100:300))/201;
Vav_absolute=[Vav_LU;Vav_MPC;Vav_MPC2];
Vav_relative=[Vav_LU/Vav_LU;Vav_MPC/Vav_LU;Vav_MPC2/Vav_LU];

Table1=[E_absolute,E_relative,Vav_absolute,Vav_relative];
disp(Table1);
fprintf('----☆----☆-----☆----☆----☆-----☆----☆----☆\n')

 %% -----☆------☆------☆------ MPC(3) AF Unaware  -----☆------☆------☆-----%%MPC3_t=25;
MPC3_tStep=0.05;
MPC3_t=25;
MPC3_tNum=MPC3_t/MPC3_tStep+1;
PlotTime3=1:MPC3_tNum;

%Objective Function
ObjFun=-sum( X_EMPC(2*State.N+3,: ) );

%Nonlinear Constraints
NonLnrConsUnA=[];
for ii = 1:State.Np+1
    NonLnrConsUnA = [NonLnrConsUnA ; X_EMPC(1:State.N-1,ii)];%state phi
    NonLnrConsUnA = [NonLnrConsUnA ; X_EMPC(State.N+3:2*State.N+1,ii)];%state Vphi
end

%NLP Struct
OptimU=reshape( U_EMPC , (State.N-1)*State.Np,1);
NlpMPCUnA=struct('f',ObjFun,'x',OptimU,'g',NonLnrConsUnA ,'p',P);

%NLP Solver
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;%Acceptable Convergence Tolerance
%opts.ipopt.acceptable_obj_change_tol = 1e-6;
MPCSolverUnA=nlpsol('MPCSolverUnA', 'ipopt', NlpMPCUnA,opts);

%Upper and lower bound constraints
argsUnA=struct;
argsUnA.lbg=[];
argsUnA.ubg=[];
for ii = 1:State.Np+1
    argsUnA.lbg = [argsUnA.lbg ; ones(State.N-1,1)* -PHIMax];
    argsUnA.lbg = [argsUnA.lbg ; ones(State.N-1,1)*  -VPHIMax];
    argsUnA.ubg = [argsUnA.ubg ; ones(State.N-1,1)* PHIMax];
    argsUnA.ubg = [argsUnA.ubg ; ones(State.N-1,1)*  VPHIMax];
end
argsUnA.lbx = -UMax*ones((State.N-1)*State.Np,1);
argsUnA.ubx = UMax*ones((State.N-1)*State.Np,1);

%Vars Initialization
PHIInit=[0,0.01,-0.01,0.01,0,0,0.01,-0.01]';
thetaInit=0;
PxInit=0;
PyInit=0;
VPHIInit=zeros(State.N-1,1);
VthetaInit=0;
VtInit=0;
VnInit=0;
XInit=[PHIInit ; thetaInit ; PxInit ; PyInit ; VPHIInit ; VthetaInit ; VtInit ; VnInit ];
UEveryOptim=zeros((State.N-1)*State.Np,1);%每一次求解送入的U
UMPCUnA=zeros(State.N-1,MPC3_tNum-1);%所有时刻的最优控制序列U
XMPCUnA=zeros(2*State.N+4,MPC3_tNum);%所有时刻的最优状态序列X
XMPCUnA( : ,1)=XInit;

%Start MPC  (3) Optim
StartOptimTime = tic;
for m=1:MPC3_tNum-1
    fprintf('MPC(3) 当前进展：第%d时刻 \n',m);
        argsUnA.p=XMPCUnA( : ,m);
        argsUnA.x0=UEveryOptim;
        OptimAnsUnA=MPCSolverUnA('x0', argsUnA.x0, 'lbx', argsUnA.lbx, 'ubx', argsUnA.ubx,'lbg', argsUnA.lbg,...
                        'ubg', argsUnA.ubg,'p',argsUnA.p);
        UEveryOptim=full(OptimAnsUnA.x);%将casADi的DM类型转化为double
        UMPCUnA(:,m)=UEveryOptim(1:State.N-1);
    if(m<200)
         XMPCUnA( : ,m+1)=MPC_Modelling(XMPCUnA( : ,m),UMPCUnA(:,m));
    else
        XMPCUnA( : ,m+1)=MPC_Modelling(XMPCUnA( : ,m),UMPCUnA(:,m));
        XMPCUnA(State.N+6,m+1)=0;
    end
    UEveryOptim(1:(end-State.N+1))=UEveryOptim(State.N:end);
end
MPCOptimTime=toc(StartOptimTime);
fprintf('MPC(3)优化结束，用时%f 秒 \n',MPCOptimTime);
clc
fprintf('----☆----☆-----☆----☆----☆-----☆----☆----☆\n')

  %% -----☆------☆------☆------ MPC(4) AF Aware  -----☆------☆------☆-----%%
%Objective Function
ObjFun=-sum( X_EMPC(2*State.N+3,: ) );

%Nonlinear Constraints
NonLnrConsAw_Unaw=[];
for ii = 1:State.Np+1
    NonLnrConsAw_Unaw = [NonLnrConsAw_Unaw ; X_EMPC(1:State.N-1,ii)];
    NonLnrConsAw_Unaw = [NonLnrConsAw_Unaw ; X_EMPC(State.N+3:2*State.N+1,ii)];
end

%NLP Struct
OptimU=reshape( U_EMPC , (State.N-1)*State.Np,1);
NlpMPC=struct('f',ObjFun,'x',OptimU,'g',NonLnrConsAw_Unaw ,'p',P);

%NLP Solver
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;%Acceptable Convergence Tolerance
%opts.ipopt.acceptable_obj_change_tol = 1e-6;
MPCSolverFormer=nlpsol('MPCSolverUnA', 'ipopt', NlpMPC,opts);
MPCSolverLater=nlpsol('MPCSolverUnA', 'ipopt', NlpMPC,opts);

%Upper and lower bound constraints
argsFormer=struct;
argsFormer.lbg=[];
argsFormer.ubg=[];
for ii = 1:State.Np+1
    argsFormer.lbg = [argsFormer.lbg ; ones(State.N-1,1)* -PHIMax];
    argsFormer.lbg = [argsFormer.lbg ; ones(State.N-1,1)*  -VPHIMax];
    argsFormer.ubg = [argsFormer.ubg ; ones(State.N-1,1)* PHIMax];
    argsFormer.ubg = [argsFormer.ubg ; ones(State.N-1,1)*  VPHIMax];
end
argsFormer.lbx = -UMax*ones((State.N-1)*State.Np,1);
argsFormer.ubx = UMax*ones((State.N-1)*State.Np,1);

%Change Constraints After 200s 
argsLater=struct;
argsLater.lbg=[];
argsLater.ubg=[];
NewVPHIConsLb= ones(State.N-1,1)*(-VPHIMax);
NewVPHIConsLb(4)=0;
NewVPHIConsUb= ones(State.N-1,1)*(VPHIMax);
NewVPHIConsUb(4)=0;

for ii = 1:State.Np+1
    argsLater.lbg = [argsLater.lbg ; ones(State.N-1,1)* -PHIMax];
    argsLater.lbg = [argsLater.lbg ; NewVPHIConsLb];
    argsLater.ubg = [argsLater.ubg ; ones(State.N-1,1)* PHIMax];
    argsLater.ubg = [argsLater.ubg ; NewVPHIConsUb];
end

NewUConsLb= -UMax*ones((State.N-1)*State.Np,1);
NewUConsLb(4:State.N-1:end)=0;
NewUConsUp= UMax*ones((State.N-1)*State.Np,1);
NewUConsUp(4:State.N-1:end)=0;

argsLater.lbx =NewUConsLb;
argsLater.ubx =NewUConsUp;

%Vars Initialization
PHIInit=[0,0.01,-0.01,0.01,0,0,0.01,-0.01]';
thetaInit=0;
PxInit=0;
PyInit=0;
VPHIInit=zeros(State.N-1,1);
VthetaInit=0;
VtInit=0;
VnInit=0;
XInit=[PHIInit ; thetaInit ; PxInit ; PyInit ; VPHIInit ; VthetaInit ; VtInit ; VnInit ];
UEveryOptim=zeros((State.N-1)*State.Np,1);%每一次求解送入的U
UMPCAw=zeros(State.N-1,MPC3_tNum-1);%所有时刻的最优控制序列U
XMPCAw=zeros(2*State.N+4,MPC3_tNum);%所有时刻的最优状态序列X
XMPCAw( : ,1)=XInit;

%Start MPC  (4) Optim
StartOptimTime = tic;
for m=1:MPC3_tNum-1
    fprintf('MPC(4) 当前进展：第%d时刻 \n',m);
    if(m<200)
        argsFormer.p=XMPCAw( : ,m);
        argsFormer.x0=UEveryOptim;
        OptimAnsFormer=MPCSolverFormer('x0', argsFormer.x0, ...
                                      'lbx', argsFormer.lbx, 'ubx', argsFormer.ubx,'lbg', ...
                                      argsFormer.lbg,'ubg', argsFormer.ubg,'p',argsFormer.p);
        UEveryOptim=full(OptimAnsFormer.x);%将casADi的DM类型转化为double
        UMPCAw(:,m)=UEveryOptim(1:State.N-1);
        XMPCAw( : ,m+1)=MPC_Modelling(XMPCAw( : ,m),UMPCAw(:,m));
    else
        argsLater.p=XMPCAw( : ,m);
        argsLater.x0=UEveryOptim;
        OptimAnsLater=MPCSolverLater('x0', argsLater.x0, ...
                                      'lbx', argsLater.lbx, 'ubx', argsLater.ubx,'lbg', ...
                                      argsLater.lbg,'ubg', argsLater.ubg,'p',argsLater.p);
        UEveryOptim=full(OptimAnsLater.x);%将casADi的DM类型转化为double
        UMPCAw(:,m)=UEveryOptim(1:State.N-1);
        XMPCAw( : ,m+1)=MPC_Modelling(XMPCAw( : ,m),UMPCAw(:,m));
        XMPCAw(State.N+6,m+1)=0;
    end
    UEveryOptim(1:(end-State.N+1))=UEveryOptim(State.N:end);
end

MPCOptimTime=toc(StartOptimTime);
fprintf('MPC(4)优化结束，用时%f 秒 \n',MPCOptimTime);
fprintf('----☆----☆-----☆----☆----☆-----☆----☆----☆\n')
fprintf('----☆----☆-----☆----☆----☆-----☆----☆----☆\n')
  %% -------☆-------☆-------☆------ PLOT (3)(4)  -----☆-------☆------☆-------%%

%Vt
figure(3)
hold on
Aw=plot(PlotTime3,XMPCAw(2*State.N+3,:),'color',[0.3,0.6,1],'linewidth',1);
UnAw=plot(PlotTime3,XMPCUnA(2*State.N+3,:),'r','linewidth',1);
plot([200 200],[0 0.07],'--','color',[1,0.5,0],'Linewidth',1.5);
legend([Aw,UnAw],'fault-aware','fault-unaware');
gcaTemp=gca;
gcaTemp.YAxis.Exponent = -2;
box on

set(gcf,'Position',[494   349   600   380])
set(gca,'xtick',[0:100:500]);
set(gca,'ytick',[0:.02:.07]);
axis([0 PlotTime3(end) 0 0.07])
xlabel('Time Steps ($\frac{1}{20}$ s)','interpreter','latex')
ylabel('Forward Velocity \it{v_t}  (m/s)')

%PHI3
figure(4)
hold on
Aw=plot(PlotTime3,XMPCAw(3,:),'color',[0.3,0.6,1],'linewidth',1);
UnAw=plot(PlotTime3,XMPCUnA(3,:),'r','linewidth',1);
plot([200 200],[-0.06 0.06],'--','color',[1,0.5,0],'Linewidth',1.5);
legend([Aw,UnAw],'fault-aware','fault-unaware');
gcaTemp=gca;
gcaTemp.YAxis.Exponent = -2;
box on

set(gcf,'Position',[494   349   600   380])
set(gca,'xtick',[0:100:500]);
set(gca,'ytick',[-0.05,0,0.05]);
axis([0 PlotTime3(end) -0.06 0.06])
xlabel('Time Steps ($\frac{1}{20}$ s)','interpreter','latex')
ylabel('Joint Distance \it{\phi_3} (m)')

%% 画图参考函数
% 
% subplot(2,2,2)
% hold on
% plot(USR.t,dtheta,'b','LineWidth',2);
% plot(USR.t,dthetaRef,'r--','LineWidth',2);
% grid on;
% xlabel('Time(s)','FontSize',15,'FontWeight','Bold')
% ylabel('$\dot{\theta}$(deg)','FontSize',15,'FontWeight','Bold','interpreter','latex')
% legend('$\dot{\theta}$','$\dot{\theta}_{ref}$','FontSize',13,'interpreter','latex') 
% title('角速度')
% hold off
