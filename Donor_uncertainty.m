clc;
clear all;

delete(gcp('nocreate'));
parpool(12);
spmd
tic
%Pauli_matrix
i = sqrt(-1);
Sigma_x = [0  1; 1  0]/2; Sigma_y = [0 -i; i  0]/2; Sigma_z = [1  0; 0 -1]/2; Identity = eye(2);
%Parameter_for_schrodinger
t_i = 0; t_f = 0.05;%0.6;

U_i = eye(4); U_t = [ 1 0 0 0 ;0 1 0 0 ;0 0 0 1 ; 0 0 1 0 ];
%Control_Parameters
B_0 = 1; gamma_e =28000; W_0 = B_0*gamma_e ; 
%Set_pulse_amplitude
%uncertainty = 0.0;
%omega_x1_arr = [-0.000226614825291575,0.000420734747623292,0.000900524478435870,9.09025557221950e-05,-4.29369920738970e-05];
%omega_y1_arr = [-0.000745284681823740,7.44470814293546e-06,0.000579036998733818,-0.000372980203188311,0.000937058631651698]; 
%uncertainty = 0.1;
%omega_x1_arr = [-0.000229562765284195,0.000418371744979030,0.000908109022758516,7.28140635265502e-05,-3.50380374589541e-05];
%omega_y1_arr = [-0.000735701681501532,1.52434836316730e-05,0.000572261730212563,-0.000375414896325241,0.000929747508235442];
%uncertainty = 0.2;
%omega_x1_arr = [-0.000233621456314741,0.000412755036554090,0.000917194037233896,3.96590199576849e-05,-1.93196269656166e-05];
%omega_y1_arr = [-0.000718322460473626,2.62367223499680e-05,0.000564304104244053,-0.000384380517641673,0.000916876323521469];



n=5 ;
xmin = -0.001;  xmax = 0.001;
omega_x1_arr = xmin + rand(1,n)*(xmax-xmin);
omega_y1_arr = xmin + rand(1,n)*(xmax-xmin);

%omega_x1_arr = [-0.000534612742183390,0.000782574185951429,0.000732626710638530,-0.000145155679564514,-6.35214404420113e-05];
%omega_y1_arr = [-0.000534882236966587,-0.000622172704841780,0.000456625310235316,0.000112906409946470,8.90789950780157e-05];

%{
%n=5
%omega_x1_arr =  1.0e-03 *[-0.132620660507666   0.374877394234259   0.715839256065953 0.131110512677927  -0.112322763419107 ];
%omega_y1_arr =  1.0e-03 *[-0.907915432237853 -0.047151372199587   0.758836609668051  -0.488523025925990 0.906888570870210];    

%uncertainty = 0.0
omega_x1_arr =  1.0e-03 *[-0.132620660507666   0.374877394234259   0.715839256065953 0.131110512677927  -0.112322763419107 ];
omega_y1_arr =  1.0e-03 *[-0.907915432237853 -0.047151372199587   0.758836609668051  -0.488523025925990 0.906888570870210];    
omega_x1_arr = [-0.000226614825291575,0.000420734747623292,0.000900524478435870,9.09025557221950e-05,-4.29369920738970e-05];
omega_y1_arr = [-0.000745284681823740,7.44470814293546e-06,0.000579036998733818,-0.000372980203188311,0.000937058631651698]; 

%uncertainty = 0.01
%omega_x1_arr = [-0.000230278968465414,0.000418443536596899,0.000859265222968262,6.06378749974294e-05,-2.12713961369501e-05];
%omega_y1_arr = [-0.000715388351447065,-2.96776350985364e-06,0.000622405877360154,-0.000425027457949524,0.000922885275836387];

%uncertainty = 0.02
omega_x1_arr = [-0.000245017022917049,0.000419358453824606,0.000874408689293908,-1.08386166524071e-07,3.16402288520951e-06];
omega_y1_arr = [-0.000669950708825652,2.57191299156154e-05,0.000610322617663502,-0.000435929150969362,0.000887027212529738];
%}

%omega_x1_arr = [0.000749866989939498,0.000418616071050761,-0.000110664760753723,0.000309841280510374,-8.09614602878843e-07,0.000507436558545293,-0.000145011758438182,0.000387922243958466];
%omega_y1_arr = [0.000268185653058170,-0.000625132635920432,-0.000806113436819898,0.000402526244896166,0.000490319000558091,-0.000371230238514204,0.000257683147462110,0.000313978534634828];
Parameters = [omega_x1_arr  omega_y1_arr ];%omega_x2_arr  omega_y2_arr ];
ParaLimL = -0.01*ones(1,length(Parameters));%1 1 1 1 1 1 1 1 1 1 1 1];
ParaLimU =  0.01*ones(1,length(Parameters));%1 1 1 1 1 1 1 1 1 1 1 1];

options = optimset('Display','iter','PlotFcns',@optimplotfval,...
    'MaxFunEvals',10000,'MaxIter',10000,'TolFun',1E-10,'TolX',1E-10);

%options = optimset('MaxFunEvals',4000,'MaxIter',4000,'TolFun',1E-10,'TolX',1E-10);
[BestParameters, BestInfidelity] = fminsearchbnd(@Cal_Infidelity,Parameters,ParaLimL,ParaLimU,options);
%[x,fval,exitflag,output] = fminsearchbnd(@Cal_Infidelity,Parameters,ParaLimL,ParaLimU,options);

disp(BestInfidelity);
disp(BestParameters);
%disp(x);
%disp(fival);
toc
end

function infidelity = Cal_Infidelity(Parameters)  
    %J1
    t_i = 0;    t_f = 0.05;%0.6;
    A_1 = 40; A_2 = 80; 
    J_0 = 60;
    
    U_i = eye(4);   U_i_r = reshape(U_i , [] , 1);
    U_t = [ 1 0 0 0 ;0 1 0 0 ;0 0 0 1 ; 0 0 1 0 ];
    
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    if isempty(odeget(opts,'Stats'))
        odeset(opts,'Stats','on');
    end
    
    [Tsol_i,Usol_i] = ode45(@(t,U) Schrodinger_H_i_rf(t,U,Parameters,t_f,A_1, A_2,J_0), [t_i t_f] , U_i_r , opts);
    uncertainty = 0.1;
    [Tsol_un,Usol_pun] = ode45(@(t,U) Schrodinger_H_pun_rf(t,U,Parameters,t_f,A_1, A_2,J_0,uncertainty), [t_i t_f] , U_i_r , opts);
    uncertainty = -0.1;
    [Tsol_un,Usol_nun] = ode45(@(t,U) Schrodinger_H_nun_rf(t,U,Parameters,t_f,A_1, A_2,J_0,uncertainty), [t_i t_f] , U_i_r , opts);
    Uf_i = Usol_i(end,:);
    Uf_i = reshape(Uf_i,[],4);
    
    Uf_pu = Usol_pun(end,:);
    Uf_pu = reshape(Uf_pu,[],4);
    
    Uf_nu = Usol_nun(end,:);
    Uf_nu = reshape(Uf_nu,[],4);
    
    infidelity_i = (1-(abs(trace(Uf_i*U_t')))^2/16 );
    infidelity_pu = (1-(abs(trace(Uf_pu*U_t')))^2/16 );
    infidelity_nu = (1-(abs(trace(Uf_nu*U_t')))^2/16 );
    
    infidelity = 10*infidelity_i + infidelity_pu + infidelity_nu;
end


function dUdt = Schrodinger_H_i_rf(t,U,Parameters,t_f,A_1,A_2,J_0)
    i = sqrt(-1); gamma_e = 28000;
    %Dt = 0;
    %global dt ;   
    %if (dt)
    %    Dt = dt;
    %end
 
    %t_f = 0.05;%0.6;
    % A_1 = hyperfine_A1(t,t_f);     A_2 = hyperfine_A2(t,t_f);
    %A_1 = 40; A_2 = 80; 
    delta_A = (A_2 - A_1)/2;    ave_A = (A_2 + A_1)/2;
    J = interaction(t,t_f,J_0);
    %J  = 40;
    
    omega_x1_arr = Parameters(1:length(Parameters)/2);
    omega_y1_arr = Parameters(length(Parameters)/2+1:length(Parameters));
    
    B_x = Control_MF_x (omega_x1_arr, t, t_f);
    B_y = Control_MF_y (omega_y1_arr, t, t_f);
    
    H = zeros(4);
    
    H(1,1) =  J/4  + delta_A ;    
    H(4,4) =  J/4  - delta_A;
    H(2,2) = ave_A/2 - J/4;%-J/2;    
    H(3,3) = -ave_A/2 -J/4 ;%- J/2;
    H(2,3) =  J/2;    H(3,2) =  J/2;
    
    H(1,2) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(1,3) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(2,1) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(2,4) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(3,1) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(3,4) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(4,2) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(4,3) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n; 
    
    H_i = H ; 
    
    
    H_i_rf = kron(eye(4),H_i);    
    %U_m = reshape(U,[],4);
    dUdt = -i * 2*pi*  H_i_rf * U ;
    %dUdt_m = -i * 2*pi*  H_p_rf * U_m ;
    %dUdt   = reshape(dUdt_m,[],1);
end



function dUdt = Schrodinger_H_pun_rf(t,U,Parameters,t_f,A_1,A_2,J_0,uncertainty)
    i = sqrt(-1); gamma_e = 28000;
    %Dt = 0;
    %global dt ;   
    %if (dt)
    %    Dt = dt;
    %end
 
    %t_f = 0.05;%0.6;
   % A_1 = hyperfine_A1(t,t_f);     A_2 = hyperfine_A2(t,t_f);
    %A_1 = 80; A_2 = 160; 
    delta_A = (A_2 - A_1)/2;    ave_A = (A_2 + A_1)/2;
    J = interaction(t,t_f,J_0);
    %uncertainty = -0.1 ; 
    %J  = 40;
    
    omega_x1_arr = Parameters(1:length(Parameters)/2);
    omega_y1_arr = Parameters(length(Parameters)/2+1:length(Parameters));
    
    B_x = Control_MF_x (omega_x1_arr, t, t_f);
    B_y = Control_MF_y (omega_y1_arr, t, t_f);
    
    H = zeros(4);
    
    H(1,1) =  J/4  + delta_A ;    
    H(4,4) =  J/4  - delta_A;
    H(2,2) = ave_A/2 - J/4;%-J/2;    
    H(3,3) = -ave_A/2 -J/4 ;%- J/2;
    H(2,3) =  J/2;    H(3,2) =  J/2;
    
    H(1,2) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(1,3) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(2,1) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(2,4) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(3,1) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(3,4) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(4,2) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(4,3) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n; 
    %%%
    J_u = J*uncertainty;
    H_u = zeros(4);
    H_u(1,1) =  J_u/4;  H_u(4,4) =  J_u/4;  H_u(2,2) = - J_u/4;
    H_u(3,3) = -J_u/4 ; H_u(2,3) =  J_u/2;  H_u(3,2) =  J_u/2;
    
    H_i = H + H_u; 
    
    
    H_un_rf = kron(eye(4),H_i);    
    %U_m = reshape(U,[],4);
    dUdt = -i * 2*pi*  H_un_rf * U ;
    %dUdt_m = -i * 2*pi*  H_p_rf * U_m ;
    %dUdt   = reshape(dUdt_m,[],1);
end


function dUdt = Schrodinger_H_nun_rf(t,U,Parameters,t_f,A_1,A_2,J_0,uncertainty)
    i = sqrt(-1); gamma_e = 28000;
    %Dt = 0;
    %global dt ;   
    %if (dt)
    %    Dt = dt;
    %end
 
    %t_f = 0.05;%0.6;
   % A_1 = hyperfine_A1(t,t_f);     A_2 = hyperfine_A2(t,t_f);
    %A_1 = 80; A_2 = 160; 
    delta_A = (A_2 - A_1)/2;    ave_A = (A_2 + A_1)/2;
    J = interaction(t,t_f,J_0);
    %uncertainty = 0.1 ; 
    %J  = 40;
    
    omega_x1_arr = Parameters(1:length(Parameters)/2);
    omega_y1_arr = Parameters(length(Parameters)/2+1:length(Parameters));
    
    B_x = Control_MF_x (omega_x1_arr, t, t_f);
    B_y = Control_MF_y (omega_y1_arr, t, t_f);
    
    H = zeros(4);
    
    H(1,1) =  J/4  + delta_A ;    
    H(4,4) =  J/4  - delta_A;
    H(2,2) = ave_A/2 - J/4;%-J/2;    
    H(3,3) = -ave_A/2 -J/4 ;%- J/2;
    H(2,3) =  J/2;    H(3,2) =  J/2;
    
    H(1,2) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(1,3) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(2,1) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(2,4) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(3,1) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(3,4) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(4,2) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(4,3) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n; 
    %%%
    J_u = J*uncertainty;
    H_u = zeros(4);
    H_u(1,1) =  J_u/4;  H_u(4,4) =  J_u/4;  H_u(2,2) = - J_u/4;
    H_u(3,3) = -J_u/4 ; H_u(2,3) =  J_u/2;  H_u(3,2) =  J_u/2;
    
    H_i = H + H_u; 
    
    
    H_un_rf = kron(eye(4),H_i);    
    %U_m = reshape(U,[],4);
    dUdt = -i * 2*pi*  H_un_rf * U ;
    %dUdt_m = -i * 2*pi*  H_p_rf * U_m ;
    %dUdt   = reshape(dUdt_m,[],1);
end


function J = interaction(t,t_f,J_0)
   J = J_0*pi/2*sin(pi*t/t_f) ;         %40*pi/2*sin(pi*t/t_f) 
end
%{
function A_1 = hyperfine_A1(t,t_f)
   A_1 = 80*pi/2*sin(pi*t/t_f) ;         
end

function A_2 = hyperfine_A2(t,t_f)
   A_2 = 160*pi/2*sin(pi*t/t_f) ;         
end
%}
function B_x = Control_MF_x (omega_x1_arr, t, t_f)    
    B_x_s_dr_1 = 0 ;
        for n = 1:length(omega_x1_arr)
            B_x_dr_1 = omega_x1_arr(n) * sin((2*n-1)*pi*t/t_f) ;
            B_x_s_dr_1 = B_x_s_dr_1 + B_x_dr_1;
        end
    B_x = B_x_s_dr_1;%* cos(w_dr*2*pi*t)  ; 
end


function B_y = Control_MF_y (omega_y1_arr, t, t_f ) 
    B_y_s_dr_1 = 0 ;
        for n = 1:length(omega_y1_arr)
            B_y_dr_1 = omega_y1_arr(n) * sin((2*n)*pi*t/t_f) ;
            B_y_s_dr_1 = B_y_s_dr_1 + B_y_dr_1;
        end
    B_y = B_y_s_dr_1;%*cos(w_dr*2*pi*t)  ; 
end






   