function render_figures(dataFile,out)
%RENDER_FIGURES 重绘32张同星表实测图、2张共用前端曲线和原理示意。
root=fileparts(mfilename('fullpath'));
if nargin<1,dataFile=fullfile(root,'spoof_receiver','outputs','report','figure_data.mat');end
assert(isfile(dataFile),'Run run_ds23 to generate figure_data.mat or pass an existing MAT file.');
if nargin<2,out=fullfile(fileparts(dataFile),'figures');end
if ~isfolder(out),mkdir(out);end
S=load(dataFile,'D','B','IQ','cl','Z');
assert(numel(S.D)==32&&nnz(startsWith({S.D.id},'sun_'))==8&&numel(S.B)==2);
for j=1:numel(S.D)
 E=S.D(j); ix=find(strcmp({S.Z.id},E.id));
 if isempty(ix)
  draw_case(E,S.cl,out,[95 240],[]);
 else
  Z=S.Z(ix);draw_case(E,S.cl,out,[95 240],Z.upper);
  E.id=['detail_' E.id]; E.title=[E.title ' | Full range detail'];E.arrows=false;
  t=E.time>=Z.start&E.time<=Z.stop;E.time=E.time(t);E.value=E.value(t);
  assert(max(E.value,[],'omitnan')==Z.peak);
  draw_case(E,S.cl,out,[Z.start Z.stop],[]);
 end
end
for j=1:numel(S.B),draw_frontend(S.B(j),S.cl,out);end
draw_iq(S.IQ,out);
draw_principles(out);
fprintf('32 same-case results, 2 front-end comparisons, 1 IQ overview and 2 schematics rendered.\n');
end

function draw_case(E,cl,out,tRange,upper)
%DRAW_CASE 上方是一种方法的检测量及门限，下方是同时间轴的I、Q。
C=E.capture; left=.12; w=.755; y=[.687 .385 .083]; h=.243;
f=figure('Visible','off','Color','w','Position',[35 50 1130 870]); ax=gobjects(1,3);
for j=1:3
    p=[left y(j) w h]; a=axes(f,'Position',p,'PositionConstraint','innerposition'); ax(j)=a;
    if j==1
        draw_score(a,E,upper);
    else
        branch=j-1; if branch==1,Z=C.I;tag='I';else,Z=C.Q;tag='Q';end
        imagesc(a,C.t,C.d,Z.','AlphaData',isfinite(Z.')); axis(a,'xy');hold(a,'on');
        colormap(a,parula(256)); clim(a,cl(branch,:)); cb=colorbar(a);
        cb.Position=[left+w+.022 y(j) .015 h]; cb.Label.String='dB';
        ylabel(a,'Delay from Prompt (chip)'); ylim(a,[-2 2]); yticks(a,-2:2);
        title(a,sprintf('%s branch | PRN %d | 81 taps',tag,C.prn)); draw_prompt(a,C,branch,E.arrows);
    end
    a.Position=p; set(a,'FontSize',11,'Layer','top','Box','on'); xlim(a,tRange);
    step=5;if diff(tRange)>60,step=20;end
    xticks(a,ceil(tRange(1)/step)*step:step:tRange(2));
end
xlabel(ax(end),'Original recording time (s)'); linkaxes(ax,'x'); drawnow;
p=vertcat(ax.Position); assert(max(abs(p(:,1)-left))<1e-10&&max(abs(p(:,3)-w))<1e-10);
exportgraphics(f,fullfile(out,[E.id '.png']),'Resolution',165); close(f);
end

function draw_score(a,E,upper)
%DRAW_SCORE 标出本图的门限来源及实际数值。
red=plot(a,E.time,E.value,'Color',[.75 .14 .11],'LineWidth',1.25); hold(a,'on');
g=yline(a,E.threshold,'--','Color',[.24 .24 .24],'LineWidth',1.1); grid(a,'on');
where='northwest'; txt=E.thresholdLabel;
if ismember(E.id,{'sun_ds2_prn32','sun_ds3_prn07','detail_sun_ds3_prn07','feng_ds3_prn07','detail_feng_ds3_prn07'}),where='northeast';end
if startsWith(E.id,'feng'),txt='Nominal chi-square (40)';end
if startsWith(E.id,'bhuiyan'),where='east';end
keys=[red g]; names={'Detector statistic',sprintf('%s = %.6g',txt,E.threshold)};
if ~isempty(upper)
  over=E.value>upper;
  o=plot(a,E.time(over),repmat(.975*upper,nnz(over),1),'.','Color',[.80 .45 .10],'MarkerSize',6);
  keys=[keys o];names=[names {'Above scale (see detail)'}]; ylim(a,[0 upper]);
end
legend(a,keys,names, ...
    'Location',where,'Box','off');
ylabel(a,E.label); title(a,E.title);
if ~any(isfinite(E.value))
 text(a,.5,.55,'No calculable pseudorange window','Units','normalized','HorizontalAlignment','center','FontSize',11);
end
if isempty(upper),lim=ylim(a); ylim(a,[lim(1),lim(2)+.06*diff(lim)]);end
end

function draw_frontend(B,cl,out)
%DRAW_FRONTEND 上方仅有一条录制级曲线，下方四颗星提供同期相关峰参照。
f=figure('Visible','off','Color','w','Position',[35 40 1240 1100]);
t=tiledlayout(f,5,2,'Padding','compact','TileSpacing','compact');
E=struct('time',B.time,'value',B.value,'threshold',B.threshold,'thresholdLabel','Baseline 99th percentile', ...
 'id',['bhuiyan_' B.record],'title',[upper(B.record) ' | Shared receiver input statistic'],'label','Pearson statistic (dB)');
a=nexttile(t,1,[1 2]);draw_score(a,E,[]);set(a,'FontSize',11,'Box','on');
xlim(a,[95 240]);xticks(a,100:20:240);
for j=1:4
 C=B.capture(j);
 for branch=1:2
  a=nexttile(t);if branch==1,Z=C.I;name='I';else,Z=C.Q;name='Q';end
  imagesc(a,C.t,C.d,Z.','AlphaData',isfinite(Z.'));axis(a,'xy');hold(a,'on');yline(a,0,'k:');
  colormap(a,parula(256));clim(a,cl(branch,:));cb=colorbar(a);cb.Label.String='dB';
  ylim(a,[-2 2]);yticks(a,-2:2);xlim(a,[95 240]);xticks(a,100:20:240);
  title(a,sprintf('PRN %d | %s',C.prn,name));set(a,'FontSize',10,'Layer','top','Box','on');
  if branch==1,ylabel(a,'Delay (chip)');end
  if j==4,xlabel(a,'Original recording time (s)');end
 end
end
exportgraphics(f,fullfile(out,[E.id '.png']),'Resolution',150);close(f);
end

function draw_prompt(a,C,b,annotate)
%DRAW_PROMPT 黑点线是Prompt；仅在代表性图中标注已分离的两条峰迹。
if ~annotate
    yline(a,0,'k:','LineWidth',1.0); return
end
t=C.t; d=C.clean;
t0=205;t1=235;ix=t>=t0&t<=t1;
assert(all(isfinite(d(ix))&abs(d(ix))<=2)&&d(t==t1)<d(t==t0));
if b==1
    plot(a,[95 t0],[0 0],'k:',[t1 240],[0 0],'k:','LineWidth',1.0);
    c=[.79 .12 .12];plot(a,[t0 t1],[0 0],'-','Color',c,'LineWidth',2.05);tip(a,t1,0,0,c);
    text(a,194,.53,'Spoof signal','Color',c,'FontWeight','bold','BackgroundColor','w','Margin',3);
else
    yline(a,0,'k:','LineWidth',1.0);
    plot(a,t(ix),d(ix),'w-','LineWidth',2.05); tip(a,t1,d(t==t1),(d(t==t1)-d(t==t1-5))/5,'w');
    text(a,194,-1.62,'Authentic signal','Color','w','FontWeight','bold','BackgroundColor',[.14 .20 .29],'Margin',3);
end
end

function tip(a,t,d,slope,c)
%TIP 箭头沿时间正方向，尖端斜率跟随真实峰的局部方向。
b=3.2;d0=d-b*slope;plot(a,[t-b t t-b],[d0-.14 d d0+.14],'-','Color',c,'LineWidth',2.05);
end

function draw_iq(IQ,out)
%DRAW_IQ 各录制原始ADC的每秒首20ms RMS，使用共用坐标范围。
f=figure('Visible','off','Color','w','Position',[35 40 1140 365]);
tl=tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
for k=1:numel(IQ)
    a=nexttile(tl);t=IQ(k).time;v=IQ(k).rms;
    plot(a,t,v(:,1),'Color',[.08 .35 .68],'LineWidth',1.05);hold(a,'on');
    plot(a,t,v(:,2),'Color',[.85 .32 .14],'LineWidth',1.05);grid(a,'on');xlim(a,[t(1),t(end)]);
    ylim(a,[550 1200]); ylabel(a,'RMS (ADC count)');title(a,upper(IQ(k).name));xlabel(a,'Recording time (s)');
    if k==1,legend(a,{'I','Q'},'Location','southeast','Box','off');end
end
title(tl,'TEXBAT | Raw I/Q amplitude');
exportgraphics(f,fullfile(out,'texbat_iq.png'),'Resolution',165);close(f);
end

function draw_principles(out)
%DRAW_PRINCIPLES 以解析三角峰示意检测量所用的特征，与实测图分开。
c=[.08 .38 .67]; r=[.75 .23 .16]; g=[.46 .49 .53]; font='Microsoft YaHei';
x=(-2:.005:2); t=(-40:40)/20; tri=@(v) max(0,1-abs(v)); taps=[-.5 0 .5];
f=figure('Visible','off','Color','w','Position',[35 40 1140 390]);
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
a=nexttile(tl); hold(a,'on'); set(a,'FontName',font,'FontSize',11);
h1=plot(a,x,tri(x),'--','Color',g,'LineWidth',1.5);
h2=plot(a,x,tri(x)+.65*tri(x+1.5),'-','Color',c,'LineWidth',2);
h3=scatter(a,taps,tri(taps),38,r,'filled');
legend(a,[h1 h2 h3],{'Prompt 缩放的单峰模板','叠加后曲线','常规早、准、迟三点'},'Location','northeast','Box','off','FontSize',10);
grid(a,'on'); xlim(a,[-2 2]);ylim(a,[0 1.3]);xlabel(a,'相对 Prompt 的码延迟（chip）');ylabel(a,'相关幅度（归一化）');
title(a,'三点可以重合，峰旁仍有差异');
a=nexttile(tl); hold(a,'on');set(a,'FontName',font,'FontSize',11);
plot(a,x,.65*tri(x+1.5),'-','Color',c,'LineWidth',1.5);
stem(a,t,.65*tri(t+1.5),'Marker','.','MarkerSize',10,'Color',c,'LineWidth',.7);
text(a,-.44,.66,'81 点，间隔 0.05 chip','FontName',font,'Color',c);
grid(a,'on');xlim(a,[-2 2]);ylim(a,[0 .8]);xlabel(a,'相对 Prompt 的码延迟（chip）');ylabel(a,'相减后的残差（归一化）');
title(a,'多相关器把差异留在残差向量中');
exportgraphics(f,fullfile(out,'principle_shang.png'),'Resolution',165);close(f);

f=figure('Visible','off','Color','w','Position',[35 40 1140 390]);
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
a=nexttile(tl);hold(a,'on');set(a,'FontName',font,'FontSize',11);
q=sind(60); iv=cosd(60);
quiver(a,0,0,.97,0,0,'Color',r,'LineWidth',2,'MaxHeadSize',.16);
quiver(a,0,0,iv,q,0,'Color',c,'LineWidth',2,'MaxHeadSize',.18);
plot(a,[iv 0],[q q],'--','Color',c); plot(a,[0 0],[0 q],'-','Color',c,'LineWidth',2);
text(a,.57,.86,'真实信号相量','Color',c,'FontName',font);
text(a,.49,.11,'欺骗信号相量','Color',r,'FontName',font);
text(a,.035,.46,'Q 投影','Color',c,'FontName',font);
axis(a,'equal');xlim(a,[-.15 1.23]);ylim(a,[-.18 1.05]);grid(a,'on');
xlabel(a,'I');ylabel(a,'Q');title(a,'相位差使真实信号投到 Q 轴');
a=nexttile(tl);hold(a,'on');set(a,'FontName',font,'FontSize',11);
h1=plot(a,x,q*tri(x+.65),'-','Color',c,'LineWidth',2);
h2=plot(a,x,q*tri(x+1.6),'--','Color',g,'LineWidth',1.6);
for z=[-.4 .4],xline(a,z,':','Color',r,'LineWidth',1.15);end
scatter(a,[-.4 .4],q*tri([-.4 .4]+.65),36,r,'filled');
text(a,-.46,.69,'E','FontWeight','bold','Color',r,'HorizontalAlignment','right');
text(a,.46,.69,'L','FontWeight','bold','Color',r);
legend(a,[h1 h2],{'真实峰靠近早点','真实峰继续退开'},'Location','northeast','Box','off');
grid(a,'on');xlim(a,[-2 2]);ylim(a,[0 1.15]);xlabel(a,'相对 Prompt 的码延迟（chip）');ylabel(a,'Q 相关幅度（归一化）');
title(a,'固定的早迟点只读到局部峰形');
exportgraphics(f,fullfile(out,'principle_sun.png'),'Resolution',165);close(f);
end
