 clear;
 close all;




% d=readtable('LCCa_PKAp.csv');
% d=readtable('RyRp.csv');


% filename = 'cAMPtot.csv';


x=dir( '*.csv');

for f = x'
	plot_single_SA(f.name)
end



function plot_single_SA(filename)

disp(filename);
d=readtable(filename);
ID= 1:1000;%d.ID+1;
para = load('para.log');

X = log(para);


% X(:,12) = [];
% Y = d.sstate_no_ISO;
% Y = d.pTime;
Y = d.peak_value;
% Y = d.minvalue;
% Y = (AP_log_dat(:, 2) + AP_log_dat(:, 2+15))/2;
% Y=log(Y);


sel=ID;

XX= X(sel,:);
YY= Y(sel,:);



[T,P,W,Wstar,U,b,C,B_pls,...
Bpls_star,Xori_rec,Yori_rec,...
R2_X,R2_Y]=PLS_nipals(XX(1:end-400,:),YY(1:end-400),rank(X));
% bar(Bpls_star)
% bar(Bpls)


% SSYT = sum((YY-ones(length(YY),1)*mean(YY)).^2);
% SSYR = sum((Yori_rec-ones(length(Yori_rec),1)*mean(Yori_rec)).^2);
% R2each = SSYR./SSYT
B_pls_all = B_pls;



XX = zscore(XX);
YY = zscore(YY);

[nn, mm ] = size(XX);

New_X = [ones(nn, 1), XX];

[b,bint,r,rint,stats] = regress(YY,New_X);
stats


yfit = New_X *b;
y=YY;
TSS = sum((y-mean(y)).^2);
RSS = sum((y-yfit).^2);
Rsquared = 1 - RSS/TSS




fig1= figure('units','inch','position',[0,0,13,5]);
bar([B_pls_all, b(2:end)]);



labels = {'LCCtot'       ,
'RyRtot'       ,
'PLBtot'       ,
'TnItot'       ,
'PLMtot'       ,
'b1ARtot'      ,
'Gstot'        ,
'ACtot'        ,
'ATP'          ,
'PDE3tot'      ,
'PDE4tot'      ,
'PKItot'       ,
'PKAIItot'     ,
'I1tot'        ,
'PP1tot'       ,
'PKACII_LCCtot',
'PP1_LCC'      ,
'PP2A_LCC'     ,
'PKAIIryrtot'  ,
'PP1ryr'       ,
'PP2Aryr'      ,
'PP1_ikstot'   ,
'PKAII_ikstot' ,
'Km_AC_basal'  ,
'Kd_AC_Gsa'    ,
'Km_PDE3_cAMP' ,
'Km_PDE4_cAMP' ,
'Km_PKA_I1'    ,
'Km_PP2A_I1'   ,
'Km_PKA_PLB'   ,
'Km_PP1_PLB'   ,
'Km_PKA_LCC'   ,
'Km_PP1_LCC'   ,
'Km_PP2A_LCC'  ,
'Km_pka_ryr'   ,
'Km_pp1_ryr'   ,
'Km_pp2a_ryr'  ,
'Km_PKA_TnI'   ,
'Km_PP2A_TnI'  ,
'Km_pka_iks'   ,
'Km_pp1_iks'   ,
'Kinetic_scale'   


};

set(gca,'XTick',1:42)
set(gca,'XTickLabel',labels)
ylabel('Sensitivity Coefficient');
str=[];
str= sprintf('in %s, R^2 = %f',filename, Rsquared);
% text(17, 0.85, str)
title(str)
xtickangle(90)


box off

set(gca,'FontSize',13, 'FontName', 'Arial') % Creates an axes and sets its FontSize to 18

% str= sprintf('%s',filename, Rsquared);
% str = strrep(filename, '.csv', [])
str = erase(filename, '.csv');
str= sprintf('Figs/%s',str);

% saveas(gcf, str);
print('-dpng','-r300',str)
close all;
% saveas(gcf,'Synergy_Ik2pIKur_B3.eps')
end