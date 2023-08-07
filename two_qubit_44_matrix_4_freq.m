clc;
clear all;

%delete(gcp('nocreate'));
%poolobj = parpool ;
%parpool local 2
%parpool(2);
%spmd

%Control_Parameters
B_0 = 1; % 1T
gamma_e = 28000;%28024.95164;%B*10^6(Hz/T);
gamma_n = 107;
A_1 =   1.0e+04 *1.552787182574073;%1440;%1440;
A_2 =1.0e+04 * 3.167945452306466;%2440;%2420;
J  = 1.0e+04 *-2.103550365443601;%;35;
%Parameter_for_schrodinger
t_i = 0;
t_f =1.0e+04 * 0.000000003739751;%3.4*10^-5  ;
U_i = eye(16);
U_t = [ 1 0 0 0 ;0 1 0 0 ;0 0 0 1 ; 0 0 1 0 ];
%Set_pulse_amplitude
%{
n=5 ;
xmin = -0.01;
xmax = 0.01;
omega_x1_1_arr = xmin + rand(1,n)*(xmax-xmin);
omega_x1_2_arr = xmin + rand(1,n)*(xmax-xmin);
omega_x1_3_arr = xmin + rand(1,n)*(xmax-xmin);
omega_x1_4_arr = xmin + rand(1,n)*(xmax-xmin);

omega_y1_1_arr = xmin + rand(1,n)*(xmax-xmin);
omega_y1_2_arr = xmin + rand(1,n)*(xmax-xmin);
omega_y1_3_arr = xmin + rand(1,n)*(xmax-xmin);
omega_y1_4_arr = xmin + rand(1,n)*(xmax-xmin);

omega_x2_1_arr = xmin + rand(1,n)*(xmax-xmin);
omega_x2_2_arr = xmin + rand(1,n)*(xmax-xmin);
omega_x2_3_arr = xmin + rand(1,n)*(xmax-xmin);
omega_x2_4_arr = xmin + rand(1,n)*(xmax-xmin);

omega_y2_1_arr = xmin + rand(1,n)*(xmax-xmin);
omega_y2_2_arr = xmin + rand(1,n)*(xmax-xmin);
omega_y2_3_arr = xmin + rand(1,n)*(xmax-xmin);
omega_y2_4_arr = xmin + rand(1,n)*(xmax-xmin);
%}
omega_x1_1_arr  = 1.0e+04 *[-0.000077878510665   0.000049950814252 -0.000050769111175   0.000036440158777  -0.000010211462892 ];
omega_x1_2_arr = 1.0e+04 *[0.000099146013447   0.000056369995338   0.000023732673947  -0.000023640207772   0.000052671669253 ];
omega_x1_3_arr = 1.0e+04 *[-0.000037922791201 0.000068852121593   0.000052113335639  -0.000007147750924  0.000080803770661];
omega_x1_4_arr = 1.0e+04 *[0.000003060307923   0.000086030423706 0.000018772686025   0.000025231726330   0.000005764716681];
omega_y1_1_arr = 1.0e+04 *[0.000086188285213   0.000009869680530   0.000051382356726  -0.000046867666781   0.000036440207227];
omega_y1_2_arr = 1.0e+04 *[0.000055478289316  -0.000034651584797   0.000014795179055   0.000004647293110  -0.000005191938845];
omega_y1_3_arr = 1.0e+04 *[0.000035814390324   0.000083313101585  0.000021805897120  -0.000067864555600   0.000034441401559];
omega_y1_4_arr = 1.0e+04 *[0.000066114975488  -0.000016479384696  -0.000019598319236 0.000082261664682   0.000048327562811] ;
omega_x2_1_arr = 1.0e+04 *[0.000086050380418   0.000027669555654  -0.000080688141319   0.000068034929578   -0.000045391746158];
omega_x2_2_arr = 1.0e+04 *[0.000026701517890   0.000002464417305  0.000063464313304   0.000065669577346   0.000020068815983];
omega_x2_3_arr = 1.0e+04 *[0.000096366977372   0.000013870841977   0.000077743344871    0.000028228476079   0.000017926660873] ;
omega_x2_4_arr = 1.0e+04 *[0.000062695194405   0.000001363756508  -0.000047868748359   0.000005184874579   0.000056821944894];
omega_y2_1_arr = 1.0e+04 *[0.000033834891997   0.000091682810666   0.000066137428783  -0.000001465925800   0.000079635770751];
omega_y2_2_arr = 1.0e+04 *[-0.000052808163887  -0.000056877190178  -0.000051424682206    -0.000096825244715  -0.000048012026294];

omega_y2_3_arr = 1.0e+04 *[0.000050534056383     -0.000017484779239   0.000030267782380   0.000070384675464   -0.000078551550493];
omega_y2_4_arr = 1.0e+04 *[ 0.000011984670479   0.000098672116843  -0.000071085014696   0.000039084232052   0.000026667705673];
Parameters = [A_1 A_2 J t_f...
    omega_x1_1_arr omega_x1_2_arr omega_x1_3_arr omega_x1_4_arr...
    omega_y1_1_arr omega_y1_2_arr omega_y1_3_arr omega_y1_4_arr...
    omega_x2_1_arr omega_x2_2_arr omega_x2_3_arr omega_x2_4_arr ...
    omega_y2_1_arr omega_y2_2_arr omega_y2_3_arr omega_y2_4_arr];

ParaLimL = -10*Parameters;
ParaLimU = 10*Parameters;        
options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',50000,'MaxIter',50000,'TolFun',1E-10,'TolX',1E-10);

[BestParameters, BestInfidelity] = fminsearchbnd(@Cal_Infidelity,Parameters,ParaLimL,ParaLimU,options);
%[x,fval,exitflag,output] = fminsearchbnd(@Cal_Infidelity,Parameters,ParaLimL,ParaLimU,options);
disp(BestInfidelity);
disp(BestParameters);
%disp(x);
%disp(fival);

%end
%delete(gcp)


function infidelity = Cal_Infidelity(Parameters)

    t_i = 0;
    t_f = Parameters(4);%3.4*10^-3 ;
    U_i = eye(4);
    U_t = [ 1 0 0 0 ;0 1 0 0 ;0 0 0 1 ; 0 0 1 0 ];
    
    U_i_r = reshape(U_i , [] , 1);
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    if isempty(odeget(opts,'Stats'))
        odeset(opts,'Stats','on');
    end
    
    [Tsol,Usol] = ode45(@(t,U) Schrodinger_H_p_rf(t,U,Parameters), [t_i t_f] , U_i_r , opts);
    
    Uf = Usol(end,:);
    Uf = reshape(Uf,[],4);
    
    %Uf = reshape(Uf,[],round(sqrt(length(Uf))));
    %    finalState = U_0'*Uf;
    F_ele = trace(sqrtm(sqrtm(U_t)*Uf*sqrtm(U_t)));
    infidelity = 1-(F_ele *F_ele')/16;
    %infidelity = 1-(abs(trace(Uf*U_t')))^2/16;
 
end


function dUdt = Schrodinger_H_p_rf(t,U,Parameters)


% t    : timing [arbitrary unit: multiple of unit time T0]
% U    : input state [16x1 flatten matrix]
% dUdt : output dU/dt [16x1 flatten matrix]
    B_0 =1;
    A_1 = Parameters(1);
    A_2 = Parameters(2);
    J = Parameters(3);
    
    omega_x1_1_arr = Parameters(5:9);
    omega_x1_2_arr = Parameters(10:14);
    omega_x1_3_arr = Parameters(15:19);
    omega_x1_4_arr = Parameters(20:24);

    omega_y1_1_arr = Parameters(25:29);
    omega_y1_2_arr = Parameters(30:34);
    omega_y1_3_arr = Parameters(35:39);
    omega_y1_4_arr = Parameters(40:44);

    omega_x2_1_arr = Parameters(45:49);
    omega_x2_2_arr = Parameters(50:54);
    omega_x2_3_arr = Parameters(55:59);
    omega_x2_4_arr = Parameters(60:64);

    omega_y2_1_arr = Parameters(65:69);
    omega_y2_2_arr = Parameters(70:74);
    omega_y2_3_arr = Parameters(75:79);
    omega_y2_4_arr = Parameters(80:84);
    
    
    gamma_e = 28000;
    
    B_x_s_1 = Control_MF_x_1(B_0, gamma_e, A_1, omega_x1_1_arr,omega_x1_2_arr,omega_x1_3_arr,omega_x1_4_arr, t);
    B_y_s_1 = Control_MF_y_1(B_0, gamma_e, A_1, omega_y1_1_arr,omega_y1_2_arr,omega_y1_3_arr,omega_y1_4_arr, t);
    B_x_s_2 = Control_MF_x_2(B_0, gamma_e, A_2, omega_x2_1_arr,omega_x2_2_arr,omega_x2_3_arr,omega_x2_4_arr, t);
    B_y_s_2 = Control_MF_y_2(B_0, gamma_e, A_2, omega_y2_1_arr,omega_y2_2_arr,omega_y2_3_arr,omega_y2_4_arr, t);

    i = sqrt(-1);
    
    U_m = reshape(U,[],4);

    %B1 = ControlField.SecondOrderField_atom1(t);
    %B2 = ControlField.SecondOrderField_atom2(t);

    W_0 = gamma_e*B_0;
    delta_A = (A_1 - A_2)/2;
    ave_A = (A_1 + A_2)/2;
    
    E_uu =  W_0 + J/4 + delta_A/2;
    E_ud = -J/4 + ave_A/2;
    E_du = -J/4 - ave_A/2;
    E_dd = -W_0 + J/4 - delta_A/2;
    H = [ E_uu 0    0    0    ;
              0    E_ud J/2  0    ;
              0    J/2  E_du 0    ;
              0    0    0    E_dd ];

    %Sigma_x = [0  1; 1  0]/2;
    %Sigma_y = [0 -i; i  0]/2;
    %Sigma_z = [1  0; 0 -1]/2;
    %Identity = eye(2);
    
    Pulse_1x =[0         0    0.5000         0;...
          0         0         0    0.5000;...
        0.5000         0         0         0;...
          0    0.5000         0         0];

    Pulse_1y =[0   0 -0.5i 0;...
                        0   0 0 -0.5i;...
                       0.5i   0 0 0;...
                      0   0.5i 0 0];
    
    Pulse_2x =[0    0.5000         0         0;...
                      0.5000         0         0         0;...
                       0         0         0    0.5000;...
                         0         0    0.5000         0];

    Pulse_2y=[0  -0.5i 0 0;...
                       0.5i 0 0 0;...
                       0 0 0 -0.5i;...
                       0 0 0.5i 0]; 
    
    
    P_1 = gamma_e*B_x_s_1*(Pulse_1x) +...
              gamma_e*B_y_s_1*(Pulse_1y)+...
              gamma_e*B_x_s_2*(Pulse_2x ) +...
              gamma_e*B_y_s_2*(Pulse_2y );
    
    
    
          
    H_p = H  + P_1 ;

%    E_w = i*(gamma_e*B0+delta_A/2)*2*pi*t;
%    U_0 = expm(E_w*(kron(Sigma_z,eye(2))+kron(eye(2),Sigma_z)));
%    dU_0dt = i*(gamma_e*B0+delta_A/2)*2*pi/T0*(kron(Sigma_z,eye(2))+kron(eye(2),Sigma_z))*U_0;
%    H_p_rf = U_0'*H_p*U_0 - i*U_0'*dU_0dt;
    H_p_rf = H_p;
    dUdt_m = -i * H_p_rf * U_m ;
    
    dUdt   = reshape(dUdt_m,[],1);
 
end


function B_x_s_1 = Control_MF_x_1 (B_0, gamma_e, A_1, omega_x1_1_arr,omega_x1_2_arr,omega_x1_3_arr,omega_x1_4_arr, t)
    %B_0 = 1;
    %gamma_e = 28000;
    w_dr_1 = gamma_e * B_0 + (A_1);
    w_dr_2 = gamma_e * B_0 + (-A_1);
    w_dr_3 = -gamma_e * B_0 + (A_1);
    w_dr_4 = -gamma_e * B_0 + (-A_1);
            
    B_x_s_dr_1 = 0 ;
        for n = 1:5
            B_x_dr_1 = omega_x1_1_arr(n) * cos(n*w_dr_1*2*pi*t) ;
            B_x_s_dr_1 = B_x_s_dr_1 + B_x_dr_1;
        end
        
        B_x_s_dr_2 = 0 ;
        for n = 1:5
            B_x_dr_2 = omega_x1_2_arr(n) * cos(n*w_dr_2*2*pi*t) ;
            B_x_s_dr_2 = B_x_s_dr_2 + B_x_dr_2;
        end
        
        B_x_s_dr_3 = 0 ;
        for n = 1:5
            B_x_dr_3 = omega_x1_3_arr(n) * cos(n*w_dr_3*2*pi*t) ;
            B_x_s_dr_3 = B_x_s_dr_3 + B_x_dr_3;
        end
        
        B_x_s_dr_4 = 0 ;
        for n = 1:5
            B_x_dr_4 = omega_x1_4_arr(n) * cos(n*w_dr_4*2*pi*t) ;
            B_x_s_dr_4 = B_x_s_dr_4 + B_x_dr_4;
        end
    B_x_s_1 = B_x_s_dr_1 + B_x_s_dr_2 + B_x_s_dr_3 + B_x_s_dr_4  ; 
end


function B_y_s_1 = Control_MF_y_1 (B_0, gamma_e, A_1, omega_y1_1_arr,omega_y1_2_arr,omega_y1_3_arr,omega_y1_4_arr, t)
    %B_0 = 1;
    %gamma_e = 28000;     
    %A_1 = 1440;
    w_dr_1 = gamma_e * B_0 + (A_1);
    w_dr_2 = gamma_e * B_0 + (-A_1);
    w_dr_3 = -gamma_e * B_0 + (A_1);
    w_dr_4 = -gamma_e * B_0 + (-A_1);
            
    B_y_s_dr_1 = 0 ;
        for n = 1:5
            B_y_dr_1 = omega_y1_1_arr(n) * cos(n*w_dr_1*2*pi*t) ;
            B_y_s_dr_1 = B_y_s_dr_1 + B_y_dr_1;
        end
        
        B_y_s_dr_2 = 0 ;
        for n = 1:5
            B_y_dr_2 = omega_y1_2_arr(n) * cos(n*w_dr_2*2*pi*t) ;
            B_y_s_dr_2 = B_y_s_dr_2 + B_y_dr_2;
        end
        
        B_y_s_dr_3 = 0 ;
        for n = 1:5
            B_y_dr_3 = omega_y1_3_arr(n) * cos(n*w_dr_3*2*pi*t) ;
            B_y_s_dr_3 = B_y_s_dr_3 + B_y_dr_3;
        end
        
        B_y_s_dr_4 = 0 ;
        for n = 1:5
            B_y_dr_4 = omega_y1_4_arr(n) * cos(n*w_dr_4*2*pi*t) ;
            B_y_s_dr_4 = B_y_s_dr_4 + B_y_dr_4;
        end
    B_y_s_1 = B_y_s_dr_1 + B_y_s_dr_2 + B_y_s_dr_3 + B_y_s_dr_4  ; 
end

function B_x_s_2 = Control_MF_x_2 (B_0, gamma_e, A_2, omega_x2_1_arr, omega_x2_2_arr, omega_x2_3_arr, omega_x2_4_arr, t)   
    %B_0 = 1;
    %gamma_e = 28000;    
    %A_2 = 2440;
    w_ur_1 = gamma_e * B_0 + (A_2);
    w_ur_2 = gamma_e * B_0 + (-A_2);
    w_ur_3 = -gamma_e * B_0 + (A_2);
    w_ur_4 = -gamma_e * B_0 + (-A_2);
    
    B_x_s_ur_1 = 0 ;
        for n = 1:5
            B_x_ur_1 = omega_x2_1_arr(n) * cos(n*w_ur_1*2*pi*t) ;
            B_x_s_ur_1 = B_x_s_ur_1 + B_x_ur_1;
        end
        
        B_x_s_ur_2 = 0 ;
        for n = 1:5
            B_x_ur_2 = omega_x2_2_arr(n) * cos(n*w_ur_2*2*pi*t) ;
            B_x_s_ur_2 = B_x_s_ur_2 + B_x_ur_2;
        end
        
        B_x_s_ur_3 = 0 ;
        for n = 1:5
            B_x_ur_3 = omega_x2_3_arr(n) * cos(n*w_ur_3*2*pi*t) ;
            B_x_s_ur_3 = B_x_s_ur_3 + B_x_ur_3;
        end
        
        B_x_s_ur_4 = 0 ;
        for n = 1:5
            B_x_ur_4 = omega_x2_4_arr(n) * cos(n*w_ur_4*2*pi*t) ;
            B_x_s_ur_4 = B_x_s_ur_4 + B_x_ur_4;
        end
    B_x_s_2 = B_x_s_ur_1 +  B_x_s_ur_2 +B_x_s_ur_3 + B_x_s_ur_4 ; 
end

function B_y_s_2 = Control_MF_y_2 (B_0, gamma_e, A_2, omega_y2_1_arr,omega_y2_2_arr,omega_y2_3_arr,omega_y2_4_arr, t)
    %B_0 = 1;
    %gamma_e = 28000;
    %A_2 = 2440;
    w_ur_1 = gamma_e * B_0 + (A_2);
    w_ur_2 = gamma_e * B_0 + (-A_2);
    w_ur_3 = -gamma_e * B_0 + (A_2);
    w_ur_4 = -gamma_e * B_0 + (-A_2);
            
    B_y_s_ur_1 = 0 ;
        for n = 1:5
            B_y_ur_1 = omega_y2_1_arr(n) * cos(n*w_ur_1*2*pi*t) ;
            B_y_s_ur_1 = B_y_s_ur_1 + B_y_ur_1;
        end
        
        B_y_s_ur_2 = 0 ;
        for n = 1:5
            B_y_ur_2 = omega_y2_2_arr(n) * cos(n*w_ur_2*2*pi*t) ;
            B_y_s_ur_2 = B_y_s_ur_2 + B_y_ur_2;
        end
        
        B_y_s_ur_3 = 0 ;
        for n = 1:5
            B_y_ur_3 = omega_y2_3_arr(n) * cos(n*w_ur_3*2*pi*t) ;
            B_y_s_ur_3 = B_y_s_ur_3 + B_y_ur_3;
        end
        
        B_y_s_ur_4 = 0 ;
        for n = 1:5
            B_y_ur_4 = omega_y2_4_arr(n) * cos(n*w_ur_4*2*pi*t) ;
            B_y_s_ur_4 = B_y_s_ur_4 + B_y_ur_4;
        end
    B_y_s_2 = B_y_s_ur_1 + B_y_s_ur_2 + B_y_s_ur_3 + B_y_s_ur_4   ; 
end

   
       