%{
***************************************************************************
* Gráfico Marcação x Tempo 
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
%%
%sinal_processado=sinal_mx./tpsw(sinal_mx); % opção para TPSW
%sinal_processado=tpsw(sinal_mx,num_amostras,10,1,1.3); % opção para TPSW
sinal_processado=sinal_mx; % opção sem TPSW
%% Conformador
Phi_max=360;
M=32; %número de sensores
d=1; r=d/2; %raio
c=1500;%velocidade do som
fs=31250;
for i=1:M
    % 
    x(i)=r*sind((i-1)*360/M); %distancias x dos sensores
    y(i)=r*cosd((i-1)*360/M); %distancias y dos sensores
end
intervalo=.1; %intervalo=25/31250; % segundos entre amostras
amostra_central=round(intervalo*31250);
tamanho_da_amostra=1000; % diametro/velocidade*fs=20.8 amostras
for tempo=1:(num_amostras/intervalo/31250)-1 %número de linhas
    sinal_agora=sinal_processado([amostra_central*(tempo)-tamanho_da_amostra/2:amostra_central*(tempo)+tamanho_da_amostra/2],:);
    disp(tempo);
  for i =1:Phi_max
    for s=1:M
       delay(i,s) = -(x(s)*sind(i)+y(s)*cosd(i))/c; 
    end
    SinalEmFase = delayseq(sinal_agora,delay(i,:),fs);
    SinalSomado(i,:) = sum(SinalEmFase');
    PressaoAcustica(i) = var(SinalSomado(i,:));
  end
  M_T(tempo,:)=PressaoAcustica;
end
%%
aa=M_T;
dB=zeros(1295,360);
dB=20*log10(aa./max(aa(:)));
figure1 = figure('Colormap',...
    [0 0 0;0 0.000708616804331541 0;0 0.00141723360866308 0;0 0.00212585041299462 0;0 0.00283446721732616 0;0 0.00354308402165771 0;0 0.00425170082598925 0;0 0.00496031763032079 0;0 0.00566893443465233 0;0 0.00637755123898387 0;0 0.00708616804331541 0;0 0.00779478438198566 0;0 0.00850340165197849 0;0 0.00921201799064875 0;0 0.00992063526064157 0;0 0.0106292515993118 0;0 0.0113378688693047 0;0 0.0120464852079749 0;0 0.0127551024779677 0;0 0.013463718816638 0;0 0.0141723360866308 0;0 0.0148809524253011 0;0 0.0308274533599615 0;0 0.0467739552259445 0;0 0.0627204552292824 0;0 0.0786669626832008 0;0 0.0946134626865387 0;0 0.110559962689877 0;0 0.126506462693214 0;0 0.142452970147133 0;0 0.15839946269989 0;0 0.174345970153809 0;0 0.190292477607727 0;0 0.206238970160484 0;0 0.222185477614403 0;0 0.23813197016716 0;0 0.254078477621078 0;0 0.285808771848679 0;0 0.317539036273956 0;0 0.349269330501556 0;0 0.380999624729156 0;0 0.412729918956757 0;0 0.444460183382034 0;0 0.476190477609634 0;0 0.509863972663879 0;0 0.543537437915802 0;0 0.577210903167725 0;0 0.610884368419647 0;0 0.64455783367157 0;0 0.678231298923492 0;0 0.711904764175415 0;0 0.745578229427338 0;0 0.77925169467926 0;0 0.812925159931183 0;0 0.831632614135742 0;0 0.850340127944946 0;0 0.86904764175415 0;0 0.88775509595871 0;0 0.906462550163269 0;0 0.925170063972473 0;0 0.943877577781677 0;0 0.962585031986237 0;0 0.981292486190796 0;0 1 0]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
image(dB,'Parent',axes1,'CDataMapping','scaled');
xlabel('Marcação (º)');
ylabel('Tempo (s)');
xlim(axes1,[0.5 360.5]);
ylim(axes1,[0.5 1295.5]);
box(axes1,'on');
set(axes1,'CLim',[-76.0748 0],'Layer','top','XTick',...
    [0 60 120 180 240 300 360],'XTickLabel',...
    {'-180','-120','-60','0','60','120','180'},'YTickLabel',...
    {'20','40','60','80','100','120'});
colorbar('peer',axes1);
