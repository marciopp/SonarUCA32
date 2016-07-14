%{
***************************************************************************
* Conformação de Feixes com Frequência Invariante 
* para Arranjo Uniforme Circular
* de 32 Hidrofones Omnidirecionais
* Marcio Pinto Pereira - junho de 2016
* Programado em Matlab R2016a com Phased Array Toolbox
* Baseado no trabalho de Fernando Pereira Gonçalves de Sá
* para Arranjos Uniformes Lineares
***************************************************************************
%}
clear all;
close all;
% Parâmetros e constantes
v=1500; % velocidade do som
L=32;   % número de sensores
J=32;   % número de taps
fs=3125; %Frequência de sincronismo decimada
fl=625./fs; %Frequencia minima da faixa analisada normalizada
fu=1250./fs; %Frequencia máxima da faixa analisada normalizada
Ts=1; % atraso entre taps
d=0.5*v/fu; % raio do conjunto - distância entre sensores normalizada 
% Geometria do arranjo
array_geometry=[zeros(1,L);zeros(1,L);zeros(1,L)];
for i=1:L
    array_geometry(1,i)=d*cosd((i-1)*360/L); %distancias x dos sensores
    array_geometry(2,i)=d*sind((i-1)*360/L); %distancias y dos sensores
end
PesoW=zeros(1024,361);
for look_direction=-180:180; % direção a conformar
    %região de frequencia invariante (manutenção da largura de feixe na banda de frequencia - tunel)
    sidelobe_theta=[-180:5:look_direction-15 look_direction+15:5:180]; 
    threshold=[1.0e-06;] %limiar gamma (flutuação da resposta do arranjo - restrição SRV)
%***************************************************************************
% Preparação para otimização
    b=sparse([-1;zeros(L*J,1)]);
    freq1=[linspace(fl,fu,21)]; % resolução em frequência
    c1=[];
    A1t=[];
    K.f=0;
    for i=1:length(freq1)
        a_look=kron(exp(-j*2*pi*freq1(i)*[0:J-1]'*Ts),spv(array_geometry,[90 look_direction],freq1(i),v));
        c1=[c1;1];
        A1t=[A1t;0 a_look'];
        K.f=K.f+1;
    end
    c2=[];
    A2t=[];
    K.q=[];
    for i=1:length(freq1)
        for m=1:length(sidelobe_theta) 
            Steering_vector=kron(exp(-j*2*pi*freq1(i)*[0:J-1]'*Ts),spv(array_geometry,[90 sidelobe_theta(m)],freq1(i),v));               
            c2=[c2;0;0];
            A2t=[A2t;-1 zeros(1,L*J);0 -Steering_vector'];
            K.q=[K.q;2];
        end
    end
    f0=fl; %frequencia de referencia  - adota a frequencia minima fl
    Q0=zeros(L*J,L*J);
    freq1=fl:0.01:fu;
    theta1=[-180:5:-5 5:5:180]; % 
    for m=1:length(theta1) 	
        for n=1:length(freq1)          
            S_Tau0=spv(array_geometry,[90 theta1(m)],f0,v);
            S_T0=exp(-j*2*pi*f0*Ts*[0:J-1]');
            Steering_vector0=kron(S_T0,S_Tau0);        
            S_Tau=spv(array_geometry,[90 theta1(m)],freq1(n),v);
            S_T=exp(-j*2*pi*freq1(n)*Ts*[0:J-1]');
            Steering_vector=kron(S_T,S_Tau);        
            Q0=Q0+real((Steering_vector-Steering_vector0)*(Steering_vector-Steering_vector0)');
        end
    end
    Q0=Q0/length(theta1)/length(freq1);
    [U0 D0]=eig(Q0);
    M0=rank(Q0);
    [sort_D0 sort_index0]=sort(diag(D0));
    D_r0=D0(sort_index0(L*J-M0+1:L*J),sort_index0(L*J-M0+1:L*J));
    U_r0=U0(:,sort_index0(L*J-M0+1:L*J));
    R0=sqrt(D_r0)*U_r0';
    c2=[c2;sqrt(threshold);zeros(size(R0,1),1)];
    A2t=[A2t;0 zeros(1,size(R0,2));zeros(size(R0,1),1) -R0];
    K.q=[K.q;size(R0,1)+1];
    epsilon=10000;
    c2=[c2;sqrt(epsilon); zeros(L*J,1)];
    A2t=[A2t; 0 zeros(1,L*J); zeros(L*J,1) eye(L*J)];
    K.q=[K.q;L*J+1];

% Chamando a otimização
    c=[c1;c2]; %Define c
    At=[A1t;A2t]; %Define At
    K.xcomplex=[1:length(c)]; %Define K
    [x,y,info]=sedumi(At,b,c,K); % Aqui acontece a otimização
    w=y(2:L*J+1); % !!!! <- Pesos ótimos dos taps do conformador

    PesoW(:,look_direction+181)=w; 
    % 1...32 = primeiro hidrofone
    % ... 992...1024 = 32º hidrofone
    % Dimensao 2 é o angulo de conformação de -180 a 180
end

