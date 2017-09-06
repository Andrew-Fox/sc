function param = fit_Sc(ID)



Mod = Get_Model(1);
Mod.t_array = linspace(0,100,100);
Mod.Cpt(1).Initial = 1;
Mod.Cpt(2).Initial = 1;
%Mod.Par = [100;2;1;3;.3;20;.1;2;1;2;1];

Mod.U.HCO3 = @(t) 15;



Mod.Data.tit_x = [0.155000000000000;0.355000000000000;0.655000000000000;1.05500000000000;1.75500000000000;2.45500000000000;0.00140000000000000;0.000120000000000000;2.70000000000000e-05;8.40000000000000e-06;2.40000000000000e-06];
Mod.Data.tit_y = [1.40000000000000;2.10000000000000;4.40000000000000;5.50000000000000;6.20000000000000;6.50000000000000;1.10000000000000;1.10000000000000;3.30000000000000;5.20000000000000;4.60000000000000];
Mod.Data.tit_sem = [0.100000000000000;0.100000000000000;0.200000000000000;0.600000000000000;0.300000000000000;0.300000000000000;0.100000000000000;0.100000000000000;0.100000000000000;0.200000000000000;0.600000000000000];


j = num2str(ID);
P = 10*rand(11,100);
Par_Guess = P(:,j);




Ig = [3 3];
func = @(t,X,Par) Mod.fun(t,X,Par,Mod,Mod.Cpt,Mod.U);
OBJ =  @(Par) Mod.obj_ss(Par,Mod,Mod.U,2,'Ca',Mod.Data,Ig,func);
options = optimset('Display','off');
param = fminsearch(OBJ,Par_Guess,options);
save(['out/par_',ID],'param')


end

