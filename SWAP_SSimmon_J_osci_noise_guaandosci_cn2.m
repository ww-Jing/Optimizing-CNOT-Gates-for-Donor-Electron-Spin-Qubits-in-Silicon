%6.290055290098358e+02 . 8.339533255001352e+02

delete(gcp('nocreate'));
clc;
i = sqrt(-1);

%{
n1 = 5;
xmin = 100; xmax = 1000;
ep_0 = xmin + rand(1,n1)*(xmax-xmin); 

tmin = 5*10^-4;%1/28000/3;
tmax = 5*10^-3; 
t_f = tmin + rand(1,1)*(tmax-tmin) ;
%}
n1 = 4;
xmin = 10^-6; xmax = 10^-5;
ep_ocsi = xmin + rand(1,(n1-1))*(xmax-xmin); 

t_f = 1/sqrt((10)^2+80^2)/4;
ep_0 = 80;
%Parameters = [ep_0]; 
Parameters = [t_f ep_0 ep_ocsi]; 
Parameters = [0.0005004887146 , 617.8737072964088, -513.8299237935033, -533.1315058477579 0 ];
Parameters = [0.000500125952479926,747.340612293612,-942.413094240474,-934.691165513210,-565.686658306661];
Parameters = [5.35702327850724e-05,23317.0304433992,99.8755617999287,99.9427237885166,-95.5068742401096];

%Parameters = [0.00511420348690857,282.861598314170,264.630702657190,846.603403325863,988.624390520767,378.595421678796];
%Parameters = [0.00627789836467371,-98.2856605120174,905.984401133760,973.298326137333,972.145120678244,470.681933218292];
%[0.0942477800093119,564.645947894514,533.089047938932,431.836823480820,326.870397047504,957.080444794357]
%[0.0319791769159354,21.5050591164855,79.6878089741870,109.133064049507,152.582232240405,216.851176859626]

%Lower bound
%Lower bound
ep_L = -100000*ones(1,1);
ep_ocsi_L = -100*ones(1,n1-1);
t_f_L = 10^-10;%1/28000/3;
%ParaLimL = [ep_L ];
ParaLimL = [t_f_L ep_L ep_ocsi_L];

%Upper bound
ep_U = 1000000*ones(1,1);
ep_ocsi_U = 100*ones(1,n1-1);
t_f_U = 10^-1 ;
%ParaLimU = [ep_U];
ParaLimU = [t_f_U ep_U ep_ocsi_U];

options = optimset('Display','iter','MaxFunEvals',150,'MaxIter' , 150,'TolFun',1E-10,'TolX',1E-10);
%options = optimset('MaxFunEvals',1500,'MaxIter' , 1500,'TolFun',1E-10,'TolX',1E-10);
I = infidelity(Parameters) ;
[BestParameters, BestInfidelity] = fminsearchbnd(@infidelity ,Parameters,ParaLimL,ParaLimU,options);
disp(BestInfidelity);
disp(BestParameters);



%%
%{
function I = infidelity(Parameters) 
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    if isempty(odeget(opts,'Stats'))
        odeset(opts,'Stats','on');
    end
%n_t = Parameters(length(Parameters));   
t_i = 0;    %t_f = Parameters(length(Parameters)); 
t_f =Parameters(1);%360/28000/3;%0.0628;%pi/50;

U_i = eye(4);   U_i_r = reshape(U_i , [] , 1);
[Tsol,Usol] = ode45(@(t,U)schrodinger(t,U,Parameters),[t_i t_f] , U_i_r , opts);

Uf_i = Usol(end,:);   Uf_i = reshape(Uf_i,[],4);
i = sqrt(-1);
U_t = [1 0 0 0; 0 1/2*(1+i) 1/2*(1-i) 0; 0 1/2*(1-i) 1/2*(1+i) 0; 0 0 0 1];
I = 1-(abs(trace(Uf_i*U_t')))^2/16;
%test forbieus norm
end
%}


function I = infidelity(Parameters)
t_i = 0;    %t_f = Parameters(length(Parameters)); 
t_f =Parameters(1);%360/28000/3;%0.0628;%pi/50;

t_slice = 1000;
t = linspace(t_i,t_f,t_slice);     
Usol = zeros(16,length(t)); 
h=(t_f - t_i)/t_slice;  
Uud_ud_arr_num = zeros(1,length(t));

k = 100;%1000
J_arr = [];
for tol = 1: k
U_i = eye(4);   
U_i_r = reshape(U_i , [] , 1);    Usol(:,1) = U_i_r; 

%gamma_ou_arr = 10^-1;%linspace(10^-8,1,25);%[10^-3/2/pi]; %linspace(10^-8,1,25);
%gamma_ou_arr = linspace(10^-5,10^0,20);
gamma_ou_arr_n = linspace(-4.4559,-3.5528 ,20);%[10^-3/2/pi]; %linspace(10^-8,1,25);
gamma_ou_arr = (10.^gamma_ou_arr_n)/10^6;
for ga = 1: length(gamma_ou_arr)
    
    gamma_ou_t = gamma_ou_arr(ga);
    beta_s = beta_arr(t, gamma_ou_t);
    %U_n+1 = f ( t_n, U_n )
    dU = @(t,U, beta_s) schrodinger(t,U, beta_s, Parameters);
for m=1:(length(t)-1)      
    k_1 = dU(t(m),Usol(:,m),beta_s(m));
    k_2 = dU((t(m)+0.5*h),(Usol(:,m)+0.5*h*k_1),beta_s(m));
    k_3 = dU((t(m)+0.5*h),(Usol(:,m)+0.5*h*k_2),beta_s(m));
    k_4 = dU((t(m)+h),(Usol(:,m)+h*k_3 ),beta_s(m));
    Usol(:,m+1) = Usol(:,m) + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*h;  % main equation
end
%{
for m=1:(length(t)-1)     
    beta_s = beta_s(m);
    dU = @(t,U) schrodinger(t,U, beta_s, Parameters);

    k_1 = dU(t(m),Usol(:,m));
    k_2 = dU((t(m)+0.5*h),(Usol(:,m)+0.5*h*k_1 ));
    k_3 = dU((t(m)+0.5*h),(Usol(:,m)+0.5*h*k_2 ));
    k_4 = dU((t(m)+h),(Usol(:,m)+h*k_3 ));
    
    Usol(:,m+1) = Usol(:,m) + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*h;  % main equation
end
%}
end
    Uf = Usol(:,end);   
    Uf_i = reshape(Uf,[],4);
    U_t = [1 0 0 0; 0 1/2*(1+i) 1/2*(1-i) 0; 0 1/2*(1-i) 1/2*(1+i) 0; 0 0 0 1];
    J_1 = (1-(((abs(trace(Uf_i*U_t')))^2)/16));
    J_arr = [J_arr J_1] ;   
end
I = sum(J_arr)/k/length(gamma_ou_arr);
end


function dU = schrodinger(t,U, beta_s, Parameters )
t_f = Parameters(1);%360/28000/3;%0.0628;%pi/50;
%B_x = 0; B_y = 0;
J_arr = Parameters(2:length(Parameters));
J_0 = interaction(J_arr ,t, t_f);
J = J_0 *(1); %J_0  *( 1 + beta_s ) ;
gamma_e = 28000; B_0 = 8; E_0 = gamma_e*B_0;
Delta_E = 10;%MHz

H = zeros(4);
i = sqrt(-1);
    H(1,1) = E_0 +Delta_E; %Delta_E/2;%E_0 + (Delta_E/2) ;    
    H(4,4) = -E_0 -Delta_E; %-Delta_E/2;%-E_0 - (Delta_E/2) ;    
    H(2,2) =  Delta_E/2 - J/2 ; %+ Delta_E/2 - J/2 ; %-J/2;    
    H(3,3) =  - Delta_E/2 - J/2; %- Delta_E/2  - J/2 ; %- J/2;
    H(2,3) =  J/2;    
    H(3,2) =  J/2;
    %{
    H(1,2) = B_x/2;%Delta_E/2;
    H(2,1) = B_x/2;%Delta_E/2;
    H(3,4) = B_x/2;%Delta_E/2;
    H(4,3) = B_x/2;%Delta_E/2;
     %}
H_un = kron(eye(4),H);    
dU = -i*2*pi*H_un*U;
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

function beta_s = beta_arr(t, gamma_ou_t )
%beta_s = 1*ones(1,length(t));; %zeros(1,length(t));% 
%beta_s = randn()*0.1*ones(1,length(t));
%0.008;%0.087518605511222*(0.056^2)/(0.3^2);%1.496882793017457e-24;%0.087518605511222*0.01663;;
S_1 = 0.0088*(10^-6*2.418*10^14*0.0011)^2;%0.551110629813258/(0.25)^2;
kappa = 0.2051;%0.316;
gamma_ou_arr_n = linspace(-6.3,1,20);%[10^-3/2/pi]; %linspace(10^-8,1,25);
gamma_ou_arr = 10.^gamma_ou_arr_n;
dominator_arr = [];
for i = 1: length(gamma_ou_arr)
    dominator = (((kappa^(log10(gamma_ou_arr(i)))*gamma_ou_arr(i)))/(gamma_ou_arr(i)^2 +1));
    dominator_arr = [dominator_arr dominator];
end
dominator = sum(dominator_arr);
%sigma_ou_t = 0.015;%0.145 
sigma_ou_t = sqrt(pi*S_1/dominator) / 10^6;   %sqrt(pi*S_1/dominator) ;       

beta_ss_t = [randn()*sigma_ou_t];
for m = 1 : length(t)-1
    dW = randn()*sqrt(t(m+1)-t(m));
    beta_ss_t(m+1) = (1-gamma_ou_t*(t(m+1)-t(m))) * beta_ss_t(m) ...
        + sigma_ou_t*sqrt(2*gamma_ou_t)*dW;
end

beta_s  =  beta_ss_t ;
end

