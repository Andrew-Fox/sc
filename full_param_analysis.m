function m =  full_param_analysis(ID)
load('Data.mat')

N_Data = 22
N_Chain = 1000;

%% for each data set:
Mod = Get_Hog_Models(str2num(ID));
fchain = N_Chain
grad_params = zeros(N_Data,Mod.Num_Par);
mcmc_params = zeros(N_Data,Mod.Num_Par);

for i=1:22
    
    % Gradient search
    [tmp] = Hog1p_Fitting(Mod,Data(i),Mod.Par_Guess,'fit');
    grad_params(i,:) = tmp ;
    % MCMC search, starting with the best parameters from the gradient search 
   % [Param,ParChain,FunChain] = Hog1p_Fitting(Mod,Data_i,tmp,'mcmc');
   % mcmc_params(i,:) = Param;
    Mod.Par_Guess = tmp;
    f = num2str(i);
    save(['out/out_',ID,f])
end
save(['out/fits_model_',ID], 'grad_params')

end 
