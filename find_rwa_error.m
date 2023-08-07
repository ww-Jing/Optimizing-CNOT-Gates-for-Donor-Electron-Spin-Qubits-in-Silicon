%6.290055290098358e+02 . 8.339533255001352e+02

delete(gcp('nocreate'));
clear all
clc;
i = sqrt(-1);


%find the best sol of both rwa and no_rwa
%{
ep_array = linspace(0,2500, 500);
B_I_rwa_arr =[];
I_no_rwa_arr = [];
B_ep0 = [];
B_t = [];

for en = 1: length(ep_array )
    
ep0 = ep_array(en) ; 
t_f = 1/sqrt((10)^2+ ep0^2)/4;
Parameters = [t_f ep0]; 

%Lower bound
t_f_L = t_f/100;%1/28000/3;
ep0_L = ep_array(en)-100;
ParaLimL = [t_f_L ep0_L];

%Upper bound
t_f_U = t_f*100 ;
ep0_U = ep_array(en)+100;
ParaLimU = [t_f_U ep0_U];

%{
n1 = 1 ; 
ep0 =500 ; 
ep_osc = 300^-3*ones(1,n1-1);
t_f = 1/sqrt((10)^2+ ep0^2)/4;
Parameters = [t_f ep0 ep_osc]; 

%Lower bound
t_f_L = 10^-6;%1/28000/3;
ep0_L = -1000;
ep_osc_L = -1000*ones(1,n1-1);
ParaLimL = [t_f_L ep0_L ep_osc_L];

%Upper bound
t_f_U = 10^-1 ;
ep0_U = 1000;
ep_osc_U = 1000*ones(1,n1-1);
ParaLimU = [t_f_U ep0_U ep_osc_U];
%}
%options = optimset('Display','iter','MaxFunEvals',100,'MaxIter' , 100,'TolFun',1E-10,'TolX',1E-10);
options = optimset('MaxFunEvals',1000,'MaxIter' , 1000,'TolFun',1E-10,'TolX',1E-10);
I_rwa = infidelity_rwa(Parameters)  ;
[BestParameters, B_I_rwa] = fminsearchbnd(@infidelity_rwa ,Parameters,ParaLimL,ParaLimU,options);
%disp(BestInfidelity);
%disp(BestParameters);
I_no_rwa = infidelity(BestParameters) ;

B_ep0 = [B_ep0 BestParameters(2)];
B_t = [B_t BestParameters(1)];
B_I_rwa_arr =[B_I_rwa_arr B_I_rwa];
I_no_rwa_arr = [I_no_rwa_arr I_no_rwa];
end
I_diff = I_no_rwa_arr-B_I_rwa_arr ; 
filt_sol = find((B_I_rwa_arr<=10^-5) & (I_no_rwa_arr<=10^-5) & (I_diff<=10^-5));
A = B_ep0(filt_sol) ; 
B = B_t (filt_sol) ; 
C = B_I_rwa_arr(filt_sol) ; 
D = I_no_rwa_arr(filt_sol) ;
E = I_diff(filt_sol) ;
Table = [A;B;C;D;E];

%Table = [B_ep0(filt_sol); B_t (filt_sol); B_I_rwa_arr(filt_sol); I_no_rwa_arr(filt_sol);  I_diff(filt_sol)];
 disp(Table);

subplot(2,1,1);
semilogy(B_ep0, B_I_rwa_arr,'r-o');
hold on; 
semilogy(B_ep0, I_no_rwa_arr,'b-o');
hold off
title('infidelity-rwa')
xlabel('J(MHz)');ylabel('Infidelity');
legend('RWA','noRWA')

subplot(2,1,2);
semilogy(B_ep0, abs(B_I_rwa_arr-I_no_rwa_arr),'g-o');
title('infidelity-difference')
xlabel('J(MHz)');ylabel('Infidelity');
toc
%disp(2*pi/(BestParameters(1)*BestParameters(2)))
%}

%find the best sol of rwa
%%{
ep_array = linspace(0,1250, 250);
B_I_rwa_arr =[];
I_no_rwa_arr = [];
B_ep0 = [];
B_t = [];

for en = 1: length(ep_array )
    
ep0 = ep_array(en) ; 
t_f = 1/sqrt((10)^2+ ep0^2)/4;
Parameters = [t_f ep0]; 

%Lower bound
t_f_L = t_f/100;%1/28000/3;
ep0_L = ep_array(en)-100;
ParaLimL = [t_f_L ep0_L];

%Upper bound
t_f_U = t_f*100 ;
ep0_U = ep_array(en)+100;
ParaLimU = [t_f_U ep0_U];

%options = optimset('Display','iter','MaxFunEvals',100,'MaxIter' , 100,'TolFun',1E-10,'TolX',1E-10);
options = optimset('MaxFunEvals',1000,'MaxIter' , 1000,'TolFun',1E-10,'TolX',1E-10);
I_rwa = infidelity_rwa(Parameters)  ;
[BestParameters, B_I_rwa] = fminsearchbnd(@infidelity_rwa ,Parameters,ParaLimL,ParaLimU,options);
%disp(BestInfidelity);
%disp(BestParameters);
%I_no_rwa = infidelity(BestParameters) ;
I_no_rwa = infidelity(BestParameters) ;

B_ep0 = [B_ep0 BestParameters(2)];
B_t = [B_t BestParameters(1)];
B_I_rwa_arr =[B_I_rwa_arr B_I_rwa];
I_no_rwa_arr = [I_no_rwa_arr I_no_rwa];
end
%I_diff = I_no_rwa_arr-B_I_rwa_arr ; 
filt_sol = find((B_I_rwa_arr<=10^-5) );
A = B_ep0(filt_sol) ; 
B = B_t (filt_sol) ; 
C = B_I_rwa_arr(filt_sol) ; 
D = I_no_rwa_arr(filt_sol);
Table = [A;B;C;D];
 disp(Table);

semilogy(B_ep0, B_I_rwa_arr,'r-o');
hold on; 
semilogy(B_ep0, I_no_rwa_arr,'b-o');
hold off
title('infidelity-rwa')
xlabel('J(MHz)');ylabel('Infidelity');
legend('RWA','noRWA')

%%}

%find the best sol of both rwa and no_rwa
%{
n1 = 3;
xmin = 0; xmax = 0;
ep_ocsi = xmin + rand(1,(n1-1))*(xmax-xmin); 
ep_0 = 9.333333333319962e+02;
t_f = 2.678478266178126e-04;%1/sqrt((10)^2+ep_0^2)/4;
%Parameters = [ep_0]; 
Parameters = [t_f ep_0 ep_ocsi]; 
%Parameters = [0.00511420348690857,282.861598314170,264.630702657190,846.603403325863,988.624390520767,378.595421678796];
%Parameters = [0.00627789836467371,-98.2856605120174,905.984401133760,973.298326137333,972.145120678244,470.681933218292];
%[0.0942477800093119,564.645947894514,533.089047938932,431.836823480820,326.870397047504,957.080444794357]
%[0.0319791769159354,21.5050591164855,79.6878089741870,109.133064049507,152.582232240405,216.851176859626]
%Parameters = [5.35702327850724e-05,23317.0304433992,99.8755617999287,99.9427237885166,-95.5068742401096];

%Lower bound
ep_L = ep_0-100;
ep_ocsi_L = -10*ones(1,n1-1);
t_f_L = t_f/100;%1/28000/3;
%ParaLimL = [ep_L ];
ParaLimL = [t_f_L ep_L ep_ocsi_L];

%Upper bound
ep_U = ep_0+100;
ep_ocsi_U = 10*ones(1,n1-1);
t_f_U = t_f*100 ;
%ParaLimU = [ep_U];
ParaLimU = [t_f_U ep_U ep_ocsi_U]; 

options = optimset('Display','iter','MaxFunEvals',1000,'MaxIter' , 1000,'TolFun',1E-8,'TolX',1E-8);
%options = optimset('MaxFunEvals',1500,'MaxIter' , 1500,'TolFun',1E-10,'TolX',1E-10);
I_rwa = infidelity_rwa(Parameters) ;
[BestParameters, BestInfidelity] = fminsearchbnd(@ infidelity_rwa ,Parameters,ParaLimL,ParaLimU,options);
disp(BestInfidelity);
disp(BestParameters);

I = infidelity(BestParameters) ;
disp(I);
%}
%{
n1 = 3; n2 = 0;n3 = 0;
xmin = 0; xmax = 0;
ep_ocsi = xmin + rand(1,(n1-1))*(xmax-xmin); 
ep_0 = 9.333333333319962e+02;
t_f = 2.678478266178126e-04;%
%ep_0 = 900;
%t_f = 1/sqrt((10)^2+ep_0^2)/4;

B_x = zeros(1,n2);
B_y = zeros(1,n3);
%Parameters = [ep_0]; 
Parameters = [t_f ep_0 ep_ocsi]; 
%Parameters = [0.00511420348690857,282.861598314170,264.630702657190,846.603403325863,988.624390520767,378.595421678796];
%Parameters = [0.00627789836467371,-98.2856605120174,905.984401133760,973.298326137333,972.145120678244,470.681933218292];
%[0.0942477800093119,564.645947894514,533.089047938932,431.836823480820,326.870397047504,957.080444794357]
%[0.0319791769159354,21.5050591164855,79.6878089741870,109.133064049507,152.582232240405,216.851176859626]
%Parameters = [5.35702327850724e-05,23317.0304433992,99.8755617999287,99.9427237885166,-95.5068742401096];
%infidelity = 10^-6vv
%Parameters =[0.000267845184924782,935.628994663012,-9.95507308995696,-9.96339932997745];

%Parameters = [0.000248711769718507,1005.41618461699,-4.65170037316463,2.59310878530346,3.13197745543566,5.29162822494773,-8.63870623278909e-05,5.57971234881673e-05 0 ];
%Parameters =[0.000157135566979210, 1590.96385541871 , ep_ocsi, B_x B_y];

%Lower bound
ep_L = ep_0-1000;
ep_ocsi_L = -10*ones(1,n1-1);
t_f_L = t_f/100;
B_x_L = -10^-4*ones(1,n2);
B_y_L = -10^-4*ones(1,n3);
ParaLimL = [t_f_L ep_L ep_ocsi_L B_x_L B_y_L];

%Upper bound
ep_U = ep_0+1000;
ep_ocsi_U = 10*ones(1,n1-1);
t_f_U = t_f*100 ;
B_x_U = 10^-4*ones(1,n2);
B_y_U = 10^-4*ones(1,n3);
ParaLimU = [t_f_U ep_U ep_ocsi_U B_x_U B_y_U]; 

options = optimset('Display','iter','MaxFunEvals',1000,'MaxIter' , 1000,'TolFun',1E-8,'TolX',1E-8);
%options = optimset('MaxFunEvals',1500,'MaxIter' , 1500,'TolFun',1E-10,'TolX',1E-10);
I_rwa = infidelity_rwa(Parameters) ;
[BestParameters, BestInfidelity] = fminsearchbnd(@ infidelity_rwa ,Parameters,ParaLimL,ParaLimU,options);
disp(BestInfidelity);
disp(BestParameters);

I = infidelity(BestParameters) ;
disp(I);
%}
%%

function I_rwa = infidelity_rwa(Parameters) 
Uf_i = evo_rwa (Parameters);
U_t = [1 0 0 0; 0 1/2*(1+i) 1/2*(1-i) 0; 0 1/2*(1-i) 1/2*(1+i) 0; 0 0 0 1];
I_rwa = 1-(abs((trace(Uf_i*U_t'))))^2/16;
end

function Uf_i = evo_rwa (Parameters)
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
    if isempty(odeget(opts,'Stats'))
        odeset(opts,'Stats','on');
    end

t_i = 0;    t_f = Parameters(1); 
U_i = eye(4);   U_i_r = reshape(U_i , [] , 1);
[Tsol,Usol] = ode45(@(t,U)schrodinger_rwa(t,U,Parameters),[t_i t_f] , U_i_r , opts);
Uf_i = Usol(end,:);   
Uf_i = reshape(Uf_i,[],4);

%trans_unitary
gamma_e = 28000; B_0 = 2;
E_0 = gamma_e*B_0;
U_trans_E0 =eye(4) ;
U_trans_E0(1,1) = exp(-i*2*pi*E_0*t_f);
U_trans_E0(4,4) = exp(i*2*pi*E_0*t_f);

Ezr = 5;
sig_x = [0 1;1 0];
Ez = i*2*pi*Ezr*t_f * sig_x ;
U_trans_Ez_sing = kron(Ez, eye(2)) + kron(eye(2),-Ez);
U_trans_Ez = expm(U_trans_Ez_sing);


Uf_i = U_trans_Ez*U_trans_E0*Uf_i;
end

%%{
function dU_rwadt = schrodinger_rwa(t,U,Parameters)
i = sqrt(-1);
t_f = Parameters(1);
J_arr = Parameters(2:length(Parameters));

J = interaction(J_arr , t, t_f);

gamma_e = 28000; B_0 = 2;
E_0 = gamma_e*B_0;
Delta_E = 5;%MHz

H = zeros(4);

    H(1,1) = 0;    
    H(4,4) = 0;    
    H(2,2) = - J/2 ;% -J/2 ; %+ Delta_E/2 - J/2 ; %-J/2;    
    H(3,3) = - J/2 ;% -J/2; %- Delta_E/2  - J/2 ; %- J/2;
    H(2,3) =  J/2;    H(3,2) =  J/2;
    
H_un = kron(eye(4),H);    
dU_rwadt = -j*2*pi*H_un*U;
end
%%}
%{
function dU_rwadt = schrodinger_rwa(t,U,Parameters)
i = sqrt(-1);
n1 = 3; n2 = 0;n3 = 0;
t_f = Parameters(1);
J_arr = Parameters(2:(1+n1));
omega_x1_arr = Parameters((2+n1):(1+n1+n2));
omega_y1_arr = Parameters((2+n1+n2):(1+n1+n2+n3));
B_x = Control_MF_x ( omega_x1_arr , t , Parameters, t_f )    ;
B_y = Control_MF_y ( omega_y1_arr , t , Parameters, t_f )    ;

J = interaction(J_arr , t, t_f);

gamma_e = 28000; B_0 = 2;
E_0 = gamma_e*B_0;
Delta_E = 5;%MHz


H = zeros(4);

    H(1,1) =  0 ;%(Delta_E/2) ;    
    H(4,4) =  0 ;%- (Delta_E/2) ;    
    H(2,2) =  - J/2 ;% -J/2 ; %+ Delta_E/2 - J/2 ; %-J/2;    
    H(3,3) =  - J/2 ;% -J/2; %- Delta_E/2  - J/2 ; %- J/2;
    H(2,3) =  J/2;    H(3,2) =  J/2;
 
    H(1,2) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(1,3) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(2,1) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(2,4) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(3,1) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(3,4) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(4,2) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(4,3) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n; 
  
H_un = kron(eye(4),H);    
dU_rwadt = -j*2*pi*H_un*U;
end
%}



%%

function I = infidelity(BestParameters) 
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    if isempty(odeget(opts,'Stats'))
        odeset(opts,'Stats','on');
    end
t_i = 0;    t_f = BestParameters(1); 

U_i = eye(4);   U_i_r = reshape(U_i , [] , 1);
[Tsol,Usol] = ode45(@(t,U)schrodinger(t,U,BestParameters),[t_i t_f] , U_i_r , opts);

Uf_i = Usol(end,:);   Uf_i = reshape(Uf_i,[],4);
i = sqrt(-1);
U_t = [1 0 0 0; 0 1/2*(1+i) 1/2*(1-i) 0; 0 1/2*(1-i) 1/2*(1+i) 0; 0 0 0 1];
I = 1-(abs((trace(Uf_i*U_t'))))^2/16;
end
%{
function dUdt = schrodinger(t,U,BestParameters)
t_f = BestParameters(1);
J_arr = BestParameters(2:length(BestParameters));
J = interaction(J_arr , t, t_f);

gamma_e = 28000; B_0 = 2.5;
E_0 = gamma_e*B_0;
Delta_E = 5;%MHz

%H = [200 J/2;J/2 -200];
H = zeros(4);

    
    H(1,1) = E_0 + (Delta_E/2); %Delta_E/2;%E_0 + (Delta_E/2) ;    
    H(4,4) =-E_0 - (Delta_E/2) ;  %-Delta_E/2;%-E_0 - (Delta_E/2) ;    
    H(2,2) =  + Delta_E/2 - J/2 ;  %+ Delta_E/2 - J/2 ; %-J/2;    
    H(3,3) =   - Delta_E/2  - J/2 ; %- Delta_E/2  - J/2 ; %- J/2;
    H(2,3) =  J/2;    H(3,2) =  J/2;
    
    i = sqrt(-1);
 
   
H_un = kron(eye(4),H);    
dUdt = -j*2*pi*H_un*U;
end
%}

function dUdt = schrodinger(t,U,BestParameters)

t_f = BestParameters(1);
J_arr = BestParameters(2:length(BestParameters));
J = interaction(J_arr , t, t_f);

gamma_e = 28000; B_0 = 2;
E_0 = gamma_e*B_0;
Delta_E = 5;%MHz

%H = [200 J/2;J/2 -200];
H = zeros(4);

    
    H(1,1) = E_0 + (Delta_E/2); %Delta_E/2;%E_0 + (Delta_E/2) ;    
    H(4,4) =-E_0 - (Delta_E/2) ;  %-Delta_E/2;%-E_0 - (Delta_E/2) ;    
    H(2,2) =  + Delta_E/2 - J/2 ;  %+ Delta_E/2 - J/2 ; %-J/2;    
    H(3,3) =   - Delta_E/2  - J/2 ; %- Delta_E/2  - J/2 ; %- J/2;
    H(2,3) =  J/2;    H(3,2) =  J/2;
    
    i = sqrt(-1);
 
   
H_un = kron(eye(4),H);    
dUdt = -j*2*pi*H_un*U;
end


function J = interaction(J_arr , t, t_f)
Parameters_J = J_arr;
J_0 = Parameters_J(1);

J_1 = 0 ;
        for n = 2: length(Parameters_J)
            B_x_dr_1 = Parameters_J(n) * (sin((2*n-1)*pi*t/t_f))^3 ;
            %B_x_dr_1 = Parameters_x(n) * (sin((n)*pi*t/t_f))^3 ;
            J_1 = J_1 + B_x_dr_1;
        end
    J = J_0 + J_1;%*cos((J/2-A_1)*t);%* cos(w_dr*2*pi*t)  ; 
end


function B_x = Control_MF_x ( omega_x1_arr , t , Parameters, t_f )    
    B_x_s_dr_1 = 0 ;
    Parameters_x = omega_x1_arr;
        for n = 1: length(Parameters_x)
            B_x_dr_1 = Parameters_x(n) * (sin((2*n-1)*pi*t/t_f))^3 ;
            %B_x_dr_1 = Parameters_x(n) * (sin((n)*pi*t/t_f))^3 ;
            B_x_s_dr_1 = B_x_s_dr_1 + B_x_dr_1;
        end
    B_x = B_x_s_dr_1;%*cos((J/2-A_1)*t);%* cos(w_dr*2*pi*t)  ; 
end

function B_y = Control_MF_y ( omega_y1_arr , t , Parameters , t_f) 
    B_y_s_dr_1 = 0 ;
    Parameters_y = omega_y1_arr;
        for n = 1: length(Parameters_y)
            B_y_dr_1 = Parameters_y(n) * (sin((2*n)*pi*t/t_f))^3 ;
            %B_y_dr_1 = Parameters_y(n) * (sin((n)*pi*t/t_f))^3 ;
            B_y_s_dr_1 = B_y_s_dr_1 + B_y_dr_1;
        end
    B_y = B_y_s_dr_1;%*cos((J/2-A_1)*t);%*cos(w_dr*2*pi*t)  ;  
end


%%


