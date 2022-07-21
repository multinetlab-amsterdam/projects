function timeseries = Jansen_network_RK2(h,T,noise,P,coupling,Nsize,SC,delay_matrix)
% generate timecourses (timeseries) based on the neural mass model of
% Jansen and Rit. In this case we model a system of N (Nsize) coupled neural
% masses. Equations are solved using a stochastic fourth order RK scheme
% Tewarie 2019

tic
% parameters JR (as in Grimbert and Faugeras 2006)
beta_E  = 100;                      % timescale for excitatory populations     ~100ms      %
beta_I  = 50;                       % timescale for inhibitory population      ~50ms       %
A       = 3.25;                     % average excitatory synaptic gain         ~3.25mV     %
B       = 22.0;                     % average inhibitory synaptic gain         ~22.0mV     %
C       = 135*[1,0.8,0.25,0.25];    % 4 dimensional vector, average number of synapses     %
                                    % between populations     ~135mV*[1,0.8,0.25,0.25]     %
nu      =  5;                       % threshold of sigmoid                     ~5s^-1      %
r       =  0.56;                    % slope of sigmoid                         ~0.56mV^-1  %
theta   =  6;                       % amplitude of sigmoid                     ~6mV        %
% P       =  200;                   % external input                           ~200mV
k       =  coupling*10;       % internodal coupling 

% simulation parameters
Tsize = ceil(T/h);
q_std = noise; % noise parameter
max_delay = ceil(max(max(delay_matrix))/h);

% initialization of variables
y = zeros(6,Nsize,Tsize);
index = 5; % noise index
suminput = zeros(Nsize,Tsize);
phe_delay = zeros(1,Nsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main loop
for n=(1+max_delay):Tsize-1

        % noise in the model
        noise = zeros(6,Nsize);
        noise(index,:) = beta_E*beta_I*q_std*randn(1,Nsize); % noise input

        % external input from other cortical populations
        suminput(:,n) = f(y(2,:,n) - y(3,:,n),nu,r,theta);
        delays = (n - ceil(delay_matrix/h));
        for m=1:Nsize
            phe_delay(m,:) = suminput(m,delays(m,:)).*SC(m,:); % external input based on sc_matrix and delays (every connection is based on distance and speed)
        end
        suminput_delay = sum(phe_delay,1);  
        
        % Stochastic Fourth order RK scheme 
        x1  = y(:,:,n);
        Fx1 = jr_equations(x1,beta_E,beta_I,A,B,C,nu,r,theta,P,Nsize,suminput_delay,k);
        w1  = noise .* sqrt(h);
        k1  = h.*Fx1 + w1;

        x2  = x1 + h/2.*k1;
        Fx2 = jr_equations(x2,beta_E,beta_I,A,B,C,nu,r,theta,P,Nsize,suminput_delay,k);
        w2  = noise .* sqrt(h);
        k2  = h.*Fx2 + w2;
        
        x3  =  x1 + h/2.*k2;
        Fx3 = jr_equations(x3,beta_E,beta_I,A,B,C,nu,r,theta,P,Nsize,suminput_delay,k);
        w3  = noise .* sqrt(h);
        k3  = h.*Fx3 + w3;
        
        x4  = x1 + h.*k3; 
        Fx4 = jr_equations(x4,beta_E,beta_I,A,B,C,nu,r,theta,P,Nsize,suminput_delay,k);
        w4  = noise .* sqrt(h);
        k4  = h.*Fx4 + w4;
        
        y(:,:,n+1) = x1 + (k1+k2.*2+k3.*2+k4).*(1/6);

%    Reference: Hansen en Penland 2006, Efficient approximate techniques for 
%    integrating stochastic differential equations
             
end
timeseries = squeeze(y(2,:,:)-y(3,:,:));
timeseries(:,1:2/h)=[];
time = toc;
fprintf('simulated timeseries Jansen Rit in %d seconds \n',time)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Jansen and Rit equations
function dy = jr_equations(y,beta_E,beta_I,A,B,C,nu,r,theta,P,Nsize,suminput,k)
dy=zeros(6,Nsize);
dy(1,:) = y(4,:); 
dy(2,:) = y(5,:); 
dy(3,:) = y(6,:); 

dy(4,:) = A.*beta_E.*(f(y(2,:) - y(3,:),nu,r,theta)) - 2.*beta_E.*y(4,:)- beta_E^2.*y(1,:);
dy(5,:) = A.*beta_E.*(P + k.*suminput + C(2).*f(C(1).*y(1,:),nu,r,theta)) - 2*beta_E.*y(5,:) - beta_E^2.*y(2,:);
dy(6,:) = B.*beta_I.*C(4).*f(C(3).*y(1,:),nu,r,theta) - 2*beta_I.*y(6,:) - beta_I^2.*y(3,:);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sigmoid function to convert potential to firing                                          
function y = f(x,nu,r,theta)
y = nu ./ (1 + exp(-r .* (x - theta))) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
