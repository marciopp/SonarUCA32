%{
***************************************************************************
* Conformação de Feixes por Atraso e Soma 
* para Arranjo Uniforme Circular
* de 32 Hidrofones Omnidirecionais
* Marcio Pinto Pereira - junho de 2016
* Programado em Matlab R2016a com Phased Array Toolbox
* Licenciado sob CC-BY-SA
***************************************************************************
%}
clc;
close all;
clear all;

%%Ler dados
num_canais      = 32;     % Número de hidrofones
num_amostras    = 4050000;  % Número de amostras para ler a cada iteração
Fs              = 31250;  % Frequência de amostragem (Hz)

[filename pathname] = uigetfile('*.bin','Selecione o arq .bin');
arquivo             = [pathname char(filename)];
fid                 = fopen(arquivo,'rb');
arquivo             = dir(arquivo);
tam_arq             = arquivo.bytes;          %tamanho do arquivo em bytes
rot_max = tam_arq/num_canais/num_amostras/4;  %limite da iteração

for rot=1:rot_max;
    sinal = fread(fid,num_amostras*num_canais,'float32');
    sinal_mx = reshape(sinal,num_canais,[])';
    sinal_mx = sinal_mx*10./2^31;                   %passando de inteiro com sinal para +-V

end
%% Conformador
M=32; %número de sensores
d=1; r=d/2; %raio
c=1500;%velocidade do som
fs=31250;
dir=150; dir360=180+dir; % direção a conformar entre -180 e 180.
intervalo=.1; %intervalo=25/31250; % segundos entre amostras
amostra_central=intervalo*31250;
tamanho_da_amostra=22; % diametro/velocidade*fs=20.8 amostras
for i=1:M
    x(i)=r*sind((i-1)*360/M); %distancias x dos sensores
    y(i)=r*cosd((i-1)*360/M); %distancias y dos sensores
end
for s=1:M
    delay(s) = -(x(s)*sind(dir360)+y(s)*cosd(dir360))/c; 
end
SinalEmFase = delayseq(sinal_mx,delay,fs);
y=sum(SinalEmFase,2); % <- Este é o sinal conformado: y !!!
yy=y; % preserva y em yy, pois rotinas LOFAR e DEMON utilizam y também
%%
% Análise LOFAR
y=yy;

inicio = 1;
fs = 3125;
fmax = fs/2;
npts = 1024;
nfft = npts*2;
novr = 0;
R = 1;

B=[];
Bm=[];
cnt=[];
Et=[];

y=y-mean(y);
if R>1
   y=decimate(y,R);
   fs=fs/3;
   fmax=fs/2;
end
B=spectrogram(y,hanning(nfft),novr,nfft,fs);
B=abs(B);
Et=var(B);
Bm=mean(B',1);
%B=log10(B);
%B=B-tpsw(B,npts+1,10,1,1.3);
B=B./tpsw(B);
B=log10(B);
B(find(B<-.2))=0;
f=(0:npts)*fmax/npts;
t=inicio/(fs*R)+(0:size(B,2)-1)*(npts-novr/2)/fmax;

   imagesc(f,t,B');
   xlabel('Frequência (Hz)')
   ylabel ('Tempo (s)')

%%
% Analise DEMON
y=yy;

nfft2=2048/2;					% Numero de pontos para calcular a FFT
Tesp=0.5;						% Tempo entre espectros apresentados
fs1=fs;							% Frequencia de amostragem do sinal original

% Demodulacao em amplitude do sinal
y=abs(y);       			% Demodula sinal

R1=25;					% Primeira decimacao
y=decimate(y,R1);
fs=fs/R1;
Fmax=fs/2;

% decimacao pelo segundo fator (dependente da faixa a ser analisada)
R2=25;
y=decimate(y,R2);
fs=fs/R2;
Fmax=fs/2;

% Calculo dos espectros DEMON, com overlap novr3
novr3=floor(nfft2-2*Fmax*Tesp);	% Calcula overlap para calculo da FFT
[Y,f,t] = spectrogram(y-mean(y),hanning(nfft2,'periodic'),novr3,nfft2,fs);
Y=abs(Y);						% Modulo da FFT
ind=(1:8);Y(ind,:)=repmat(Y(length(ind),:),[length(ind) 1]); % Descarta 8 primeiros bins
Y=Y./tpsw(Y);					% Normaliza usando TPSW

% Apresenta espectros DEMON
figure						% Abre nova janela		
subplot(2,1,1)				% Seleciona parte superior da figura
imagesc(f*60,t,Y')			% Desenha DEMONgrama
xlabel('Rotação (rpm)','fontsize',14)		% Legenda do Eixo x
ylabel('Tempo (s)','fontsize',14)         % Legenda do Eixo y
colormap(1-gray)			% Escala de cores

subplot(2,1,2)				% Seleciona parte inferior da figura
plot(f*60,10*log10(normaliza(mean(Y'),0))),grid			% Plota espectro DEMON medio

axis tight					% Ajusta eixos
xlabel('Rotação (rpm)','fontsize',14)		% Legenda do Eixo x
ylabel('Amplitude','fontsize',14)
