%{
***************************************************************************
* Padrão de Diretividade para Arranjo Uniforme Circular 
* de 32 Hidrofones Omnidirecionais
* Broadside a 0º
* Marcio Pinto Pereira - junho de 2016
* Programado em Matlab R2016a com Phased Array Toolbox
* Licenciado sob CC-BY
***************************************************************************
%}
h = phased.UCA; % Arranjo Cicular Uniforme (UCA)
h.NumElements = 32; % Quantidade de Hidrofones
h.Radius = 0.5; % Raio do Arranjo
h.ArrayNormal = 'z'; % Definição do eixo de elevação
% Calculando Taper
wind = ones(1,32);
h.Taper = wind;
% Definindo Hidrofone Omnidirecional
el = phased.OmnidirectionalMicrophoneElement; 
h.Element = el;
% Direção a conformar: [Azimute;Elevação]
SA = [0;0]; 
% Quantização de bits de deslocamento de fase
PSB = 0; 
% Velocidade de Propagação
PS = 1500;
cont=0;
% Frequências a calcular
freq=[1 100:100:10000]; 
for f=1:length(freq) 
    F=freq(f);
    disp(F)
    NumCurves = length(F);
    %Cálculo dos pesos
    w = zeros(getDOF(h), NumCurves);
    for idx = 1:length(F)
        SV = phased.SteeringVector('SensorArray',h, 'PropagationSpeed', PS, ...
        'NumPhaseShifterBits', PSB(idx));
        w(:, idx) = step(SV, F(idx), SA(:, idx));
    end
    fmt = 'rectangular';
    cutAngle = 0;
    cont=cont+1;
    pat(:,cont)=pattern(h, F, -180:180, cutAngle, 'PropagationSpeed', PS, 'Type', ...
    'directivity', 'CoordinateSystem', fmt ,'weights', w);
end

%%
% Plota padrão
patt=pat-max(pat(:));
figure1 = figure('Colormap',...
    [0 0 0.5625;0 0 0.57267439365387;0 0 0.582848846912384;0 0 0.593023240566254;0 0 0.603197693824768;0 0 0.613372087478638;0 0 0.623546540737152;0 0 0.633720934391022;0 0 0.643895328044891;0 0 0.654069781303406;0 0 0.664244174957275;0 0 0.67441862821579;0 0 0.684593021869659;0 0 0.694767415523529;0 0 0.704941868782043;0 0 0.715116262435913;0 0 0.725290715694427;0 0 0.735465109348297;0 0 0.745639562606812;0 0 0.755813956260681;0 0 0.765988349914551;0 0 0.776162803173065;0 0 0.786337196826935;0 0 0.796511650085449;0 0 0.806686043739319;0 0 0.816860437393188;0 0 0.827034890651703;0 0 0.837209284305573;0 0 0.847383737564087;0 0 0.857558131217957;0 0 0.867732584476471;0 0 0.877906978130341;0 0 0.88808137178421;0 0 0.898255825042725;0 0 0.908430218696594;0 0 0.918604671955109;0 0 0.928779065608978;0 0 0.938953459262848;0 0 0.949127912521362;0 0 0.959302306175232;0 0 0.969476759433746;0 0 0.979651153087616;0 0 0.98982560634613;0 0 1;0 0.200000002980232 1;0 0.400000005960464 1;0 0.600000023841858 1;0 0.800000011920929 1;0 1 1;0.200000002980232 1 0.800000011920929;0.400000005960464 1 0.600000023841858;0.600000023841858 1 0.400000005960464;0.800000011920929 1 0.200000002980232;1 1 0;1 0.800000011920929 0;1 0.600000023841858 0;1 0.400000005960464 0;1 0.200000002980232 0;1 0 0;0.899999976158142 0 0;0.800000011920929 0 0;0.699999988079071 0 0;0.600000023841858 0 0;0.5 0 0]);
axes1 = axes('Parent',figure1,...
    'Position',[0.258874388254487 0.126318788305677 0.646125611745515 0.798681211694324]);
hold(axes1,'on');
mesh(patt,'Parent',axes1);
xlabel('Frequência (Hz)');
ylabel('Direção de Chegada (graus)');
xlim(axes1,[0 100]);
ylim(axes1,[0 360]);
view(axes1,[-37.5 30]);
grid(axes1,'on');
set(axes1,'CLim',[-60 0],'XAxisLocation','origin','XTick',...
    [0 20 40 60 80 100],'XTickLabel',{'0','2000','4000','6000','8000','10000'},...
    'YAxisLocation','origin','YTick',[0 90 180 270 360],'YTickLabel',...
    {'180','90','0','-90','-180'});
colorbar('peer',axes1,'Position',...
    [0.109298531810768 0.401709401709404 0.0603588907014684 0.369658119658119]);
