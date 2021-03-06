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
clear all;
M=32; % número de sensores
d=1; r=d/2; % diametro em m
c=1500;% velocidade do som em m/s
dir=90; % direção a olhar em º 0 a 360 - broadside é 0
npts=360;
angulo=[-180:360/npts:180];
cont=0;
somareal=0;somaimag=0;
freq=[100:100:10000];
for i=1:M
    x(i)=r*sind((i-1)*360/M); %distancias x dos sensores
    y(i)=r*cosd((i-1)*360/M); %distancias y dos sensores
    delay_dir(i) = -((x(i)*sind(dir))+(y(i)*cosd(dir)))/c;
end
for f=1:length(freq)
    disp(freq(f));
    cont=cont+1;
    for i =1:length(angulo)
        somareal=0;somaimag=0;
        for s=1:M
            % Delay da onda plana
            delay_angulo = -((x(s)*sind(angulo(i)))+(y(s)*cosd(angulo(i))))/c;
            % O delay que maximiza a direção a observar está travado no sensor.
            % A fase será a diferença entre o delay da direção e o
            %   delay real de chegada da onda plana.
            fase =(freq(f))*(delay_dir(s)-delay_angulo); 
            real(i,s)= cosd(360*fase); imag(i,s)= sind(360*fase);
            somareal= somareal + real(i,s);
            somaimag= somaimag + imag(i,s);
        end
    %disp(somareal); disp(somaimag);
    g(i)=sqrt(somareal^2+somaimag^2)/32;
    dB(i)=20*log10(g(i));
    end
    diag(:,cont)=dB;
    diag2(:,cont)=g;
end
%%
figure;
xlim([1 360]);
%axes1 = axes('Parent',Parent1);
%hold(axes,'on');
title('Diretividade');
%box(axes,'on');
%set(axes,'XGrid','on','XMinorTick','on','XTick',...
%    [0 360 720 1080 1440 1800 2160 2520 2880 3240 3600],'YGrid','on');
plot(dB);

figure;
theta = [-pi:2*pi/npts:pi];
polarplot(theta,dB-min(dB),'--r');

figure;
theta = [-pi:2*pi/npts:pi];
polarplot(theta,g,'--r');

figure1 = figure('Colormap',...
    [0 0 0.5625;0 0 0.57267439365387;0 0 0.582848846912384;0 0 0.593023240566254;0 0 0.603197693824768;0 0 0.613372087478638;0 0 0.623546540737152;0 0 0.633720934391022;0 0 0.643895328044891;0 0 0.654069781303406;0 0 0.664244174957275;0 0 0.67441862821579;0 0 0.684593021869659;0 0 0.694767415523529;0 0 0.704941868782043;0 0 0.715116262435913;0 0 0.725290715694427;0 0 0.735465109348297;0 0 0.745639562606812;0 0 0.755813956260681;0 0 0.765988349914551;0 0 0.776162803173065;0 0 0.786337196826935;0 0 0.796511650085449;0 0 0.806686043739319;0 0 0.816860437393188;0 0 0.827034890651703;0 0 0.837209284305573;0 0 0.847383737564087;0 0 0.857558131217957;0 0 0.867732584476471;0 0 0.877906978130341;0 0 0.88808137178421;0 0 0.898255825042725;0 0 0.908430218696594;0 0 0.918604671955109;0 0 0.928779065608978;0 0 0.938953459262848;0 0 0.949127912521362;0 0 0.959302306175232;0 0 0.969476759433746;0 0 0.979651153087616;0 0 0.98982560634613;0 0 1;0 0.200000002980232 1;0 0.400000005960464 1;0 0.600000023841858 1;0 0.800000011920929 1;0 1 1;0.200000002980232 1 0.800000011920929;0.400000005960464 1 0.600000023841858;0.600000023841858 1 0.400000005960464;0.800000011920929 1 0.200000002980232;1 1 0;1 0.800000011920929 0;1 0.600000023841858 0;1 0.400000005960464 0;1 0.200000002980232 0;1 0 0;0.899999976158142 0 0;0.800000011920929 0 0;0.699999988079071 0 0;0.600000023841858 0 0;0.5 0 0]);
axes1 = axes('Parent',figure1,...
    'Position',[0.258874388254487 0.126318788305677 0.646125611745515 0.798681211694324]);
hold(axes1,'on');
mesh(diag,'Parent',axes1);
xlabel('Frequência (Hz)');
ylabel('Direção de Chegada (graus)');
xlim(axes1,[0 100]);
ylim(axes1,[0 360]);
view(axes1,[-37.5 30]);
grid(axes1,'on');
set(axes1,'CLim',[-60 0],'XAxisLocation','origin','XTick',...
    [0 20 40 60 80 100],'XTickLabel',{'0','2000','4000','6000','8000','10000'},...
    'YAxisLocation','origin','YTick',[0 90 180 270 360],'YTickLabel',...
    {'-180','-90','0','90','180'});
colorbar('peer',axes1,'Position',...
    [0.109298531810768 0.401709401709404 0.0603588907014684 0.369658119658119]);
