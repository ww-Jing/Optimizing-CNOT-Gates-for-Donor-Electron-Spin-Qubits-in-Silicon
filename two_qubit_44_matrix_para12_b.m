clc;
clear all;

delete(gcp('nocreate'));
%poolobj = parpool ;
%parpool local 2
parpool(12);
spmd

%Pauli_matrix
i = sqrt(-1);
Sigma_x = [0  1; 1  0]/2;
Sigma_y = [0 -i; i  0]/2;
Sigma_z = [1  0; 0 -1]/2;
Identity = eye(2);


%Control_Parameters
B_0 = 1; % 1T
gamma_e =28000;% 28000;%28024.95164;%B*10^6(Hz/T);
W_0 = B_0*gamma_e ; 

A_1 = 80;%80;%(MHZ)
A_2 = 120;%120;
J  = 10;%(80)
%Parameter_for_schrodinger
t_i = 0;
t_f = 0.1;%0.6;
U_i = eye(4);
U_t = [ 1 0 0 0 ;0 1 0 0 ;0 0 0 1 ; 0 0 1 0 ];
%Set_pulse_amplitude

n=8 ;
xmin = -0.0001;
xmax = 0.0001;
omega_x1_arr = xmin + rand(1,n)*(xmax-xmin);
omega_y1_arr = xmin + rand(1,n)*(xmax-xmin);
%omega_x2_arr = xmin + rand(1,n)*(xmax-xmin);
%omega_y2_arr = xmin + rand(1,n)*(xmax-xmin);
%{
Parameters = 1.0e-4 *[0.350861424310527   0.364875488779825 0.332028042688656   0.313898164080148...
    0.347642923309482   0.336469376443514  0.378990301931722   0.378148194573606...
    0.367095006718500   0.313315280486129 0.302131099614483   0.356136817962639...
    0.330060397066256   0.411437677388939 0.398016279731130   0.328496177025165...
    0.380197165617752   0.389693243943765 0.359612540805857   0.387568615643843];
%}

Parameters = [omega_x1_arr  omega_y1_arr ];%omega_x2_arr  omega_y2_arr ];

ParaLimL = -0.001*ones(1,16);%1 1 1 1 1 1 1 1 1 1 1 1];
ParaLimU =  0.001*ones(1,16);%1 1 1 1 1 1 1 1 1 1 1 1];
options = optimset('Display','iter','PlotFcns',@optimplotfval,...
    'MaxFunEvals',10000,'MaxIter',10000,'TolFun',1E-10,'TolX',1E-10);

[BestParameters, BestInfidelity] = fminsearchbnd(@Cal_Infidelity,Parameters,ParaLimL,ParaLimU,options);
%[x,fval,exitflag,output] = fminsearchbnd(@Cal_Infidelity,Parameters,ParaLimL,ParaLimU,options);
disp(BestInfidelity);
disp(BestParameters);
%disp(x);
%disp(fival);

end
%delete(gcp)


function infidelity = Cal_Infidelity(Parameters)

    t_i = 0;
    t_f = 0.1;%0.6;
    U_i = eye(4);
    U_i_r = reshape(U_i , [] , 1);
    U_t = [ 1 0 0 0 ;0 1 0 0 ;0 0 0 1 ; 0 0 1 0 ];
    
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    if isempty(odeget(opts,'Stats'))
        odeset(opts,'Stats','on');
    end
    
    [Tsol,Usol] = ode45(@(t,U) Schrodinger_H_p_rf(t,U,Parameters), [t_i  t_f] , U_i_r , opts);
    
    Uf = Usol(end,:);
    Uf = reshape(Uf,[],4);
    
    %Uf = reshape(Uf,[],round(sqrt(length(Uf))));
    %    finalState = U_0'*Uf;
    %Uf = Uf'*Uf;
    %U_t = U_t'*U_t;
    %F_ele = trace(sqrtm(sqrtm(U_t)*Uf*sqrtm(U_t)));
    %infidelity = 1-sqrt(conj(F_ele) *F_ele)/16;
    infidelity = 1-(abs(trace(Uf*U_t')))^2/16;
 
end


function dUdt = Schrodinger_H_p_rf(t,U,Parameters)
    i = sqrt(-1);
    gamma_e = 28000;
    A_1 = 80;
    A_2 = 120;
    delta_A = (A_2 - A_1)/2;
    ave_A = (A_2 + A_1)/2;
    J  = 10;
    t_f = 0.1;%0.6;
    
    
    omega_x1_arr = Parameters(1:8);
    omega_y1_arr = Parameters(9:16);
    %omega_x2_arr = Parameters(13:18);
    %omega_y2_arr = Parameters(19:24);
    B_x_s_1 = Control_MF_x_1(omega_x1_arr, t, t_f , A_1);
    B_y_s_1 = Control_MF_y_1(omega_y1_arr, t, t_f , A_1);
    B_x_s_2 = Control_MF_x_2(omega_x1_arr, t, t_f , A_2);
    B_y_s_2 = Control_MF_y_2(omega_y1_arr, t, t_f , A_2);
    
    H = zeros(4);
    
    H(1,1) =  J/4 ; %delta_A + J/4;
    H(4,4) =  J/4 ; %-delta_A + J/4;
    
    H(2,2) = ave_A/2 - J/4;%-J/2;
    H(3,3) = -ave_A/2 -J/4 ;%- J/2;
    
    H(2,3) =  J/2;
    H(3,2) =  J/2;
    %{
    H(1,2) = 1/4*gamma_e*(B_x_s_2/2 -  (B_y_s_2*i)/2)*RF_com_p_2;
    H(1,3) = 1/4*gamma_e*(B_x_s_1/2 -  (B_y_s_1*i)/2)*RF_com_n_1;
    H(2,1) = 1/4*gamma_e*(B_x_s_2/2 + (B_y_s_2*i)/2)*RF_com_n_2;
    H(2,4) = 1/4*gamma_e*(B_x_s_1/2 -  (B_y_s_1*i)/2)*RF_com_n_1;
    H(3,1) = 1/4*gamma_e*(B_x_s_1/2 + (B_y_s_1*i)/2)*RF_com_p_1;
    H(3,4) = 1/4*gamma_e*(B_x_s_2/2 -  (B_y_s_2*i)/2)*RF_com_p_2;
    H(4,2) = 1/4*gamma_e*(B_x_s_1/2 + (B_y_s_1*i)/2)*RF_com_p_1;
    H(4,3) = 1/4*gamma_e*(B_x_s_2/2 + (B_y_s_2*i)/2)*RF_com_n_2;
    %}
   
    %%{
    H(1,2) = 1/2*gamma_e*(B_x_s_2/2 - (B_y_s_2)*i/2);%*RF_com_p;
    H(1,3) = 1/2*gamma_e*(B_x_s_1/2 - (B_y_s_1)*i/2);%*RF_com_p;
    H(2,1) = 1/2*gamma_e*(B_x_s_2/2 + (B_y_s_2)*i/2);%*RF_com_n;
    H(2,4) = 1/2*gamma_e*(B_x_s_1/2 - (B_y_s_1)*i/2);%*RF_com_p;
    H(3,1) = 1/2*gamma_e*(B_x_s_1/2 + (B_y_s_1)*i/2);%*RF_com_n;
    H(3,4) = 1/2*gamma_e*(B_x_s_2/2 - (B_y_s_2)*i/2);%*RF_com_p;
    H(4,2) = 1/2*gamma_e*(B_x_s_1/2 + (B_y_s_1)*i/2);%*RF_com_n;
    H(4,3) = 1/2*gamma_e*(B_x_s_2/2 + (B_y_s_2)*i/2);%*RF_com_n; 
    %%}
    
    H_p_rf = kron(eye(4),H);
    
    %U_m = reshape(U,[],4);
    dUdt = -i * 2*pi*  H_p_rf * U ;
    %dUdt_m = -i * 2*pi*  H_p_rf * U_m ;
    
    %dUdt   = reshape(dUdt_m,[],1);
 
end


function B_x_s_1 = Control_MF_x_1 (omega_x1_arr, t, t_f , A_1)
    %B_0 = 1;
    %gamma_e = 28000;
    w_dr = - (A_1);%+W_0   
    B_x_s_dr_1 = 0 ;
        for n = 1:8
            B_x_dr_1 = omega_x1_arr(n) * sin((2*n-1)*pi*t/t_f) ;
            B_x_s_dr_1 = B_x_s_dr_1 + B_x_dr_1;
        end
    B_x_s_1 = B_x_s_dr_1;%* cos(w_dr*2*pi*t)  ; 
end


function B_y_s_1 = Control_MF_y_1 (omega_y1_arr, t, t_f , A_1)
    %B_0 = 1;
    %gamma_e = 28000;     
    %A_1 = 1440;
    w_dr = - (A_1);%+W_0    
    B_y_s_dr_1 = 0 ;
        for n = 1:8
            B_y_dr_1 = omega_y1_arr(n) * sin((2*n)*pi*t/t_f) ;
            B_y_s_dr_1 = B_y_s_dr_1 + B_y_dr_1;
        end
    B_y_s_1 = B_y_s_dr_1;%*cos(w_dr*2*pi*t)  ; 
end

function B_x_s_2 = Control_MF_x_2 (omega_x1_arr, t, t_f , A_2)   
    %B_0 = 1;
    %gamma_e = 28000;    
    %A_2 = 2440;
    w_ur = (A_2);%+W_0
    B_x_s_dr_2 = 0 ;
        for n = 1:8
            B_x_dr_2 = omega_x1_arr(n) * sin((2*n-1)*pi*t/t_f);
            B_x_s_dr_2 = B_x_s_dr_2 + B_x_dr_2;
        end
    B_x_s_2 = B_x_s_dr_2;%*cos(w_ur*2*pi*t)   ; 
end

function B_y_s_2 = Control_MF_y_2 (omega_y1_arr, t, t_f , A_2)
    %B_0 = 1;
    %gamma_e = 28000;
    %A_2 = 2440;
    w_ur = (A_2);%+W_0     
    B_y_s_dr_2 = 0 ;
        for n = 1:8
            B_y_dr_2 = omega_y1_arr(n) * sin((2*n)*pi*t/t_f) ;
            B_y_s_dr_2 = B_y_s_dr_2 + B_y_dr_2;
        end
    B_y_s_2 = B_y_s_dr_2;%*cos(w_ur*2*pi*t)   ; 
end