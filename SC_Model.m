classdef SC_Model
    
    properties
        Par = struct();            
        Cpt = struct();
        t_array = linspace(0,20,1000);
        Data = struct();
        Param_guess = struct();
        fs_guess = struct();
        U = struct();
        sigma_D = struct();
        Index = struct();
        N = struct();
        I = struct();
        func = struct();
    end
    
    
    properties (Constant)
    end
    
    properties (Dependent)
        Xo
        Numerical
        Parameters_TE
        Parameters_SS
        
    end
    
    methods
        
        function Xo = get.Xo(obj)    % Define Initial Condition for dy/dt
             for i = 1:length(obj.Cpt)
                Xo(i) = obj.Cpt(i).Initial;    % Initial condition
            end
        end
        
        
        function Trajectory = get.Numerical(obj) % Solve dy/dt for Par property
            
            Trajectory = obj.solveODE(obj,obj.Par);
            
        end  
        

        function fit = get.Parameters_TE(obj) % Fit Parameters for Temporal Evolution Data
            
            options = optimset('Display','off');
            %func = @(t,X,Par) obj.fun(t,X,Par,obj,obj.Cpt,obj.U);
            OBJ = @(Param) obj.objective_TE(obj.Data,obj.func,obj.t_array,obj.Xo,Param,obj.Index);
            fit  = fminsearch(OBJ,obj.Param_guess,options);
            fit = fit';
        end
        
        
       
        

            

    end
    
    methods (Static)
        
        function out = fun(t,X,K,obj,Cpt,U) % Define dy/dt as y(i) = Input - Output
             
            for i = 1:length(Cpt)
                 f(i) = (Cpt(i).Input(obj,X,t,K,U) - Cpt(i).Output(obj,X,t,K,U)); % right hand side of ODE function;
            end
            out = f';
        end
        
        function Trajectory = solveODE(obj,Par) % Solve dy/dt for any parameter set
            
            func = @(t,X) obj.fun(t,X,Par,obj,obj.Cpt,obj.U); % Right hand side of ODE function
            
            [~,Trajectory] = ode23s(func,obj.t_array,(obj.Xo)); % Solve ode for times in t_array
        end
        
      
        
        function [J,x_Mod] = objective_TE(Data,func,t_array,Xo,Par,Index) % Objective function for temporal evolution fitting
            
            func = @(t,X) func(t,X,Par);
            [~,Trajectory] = ode23s(func,t_array,(Xo));
            x_Mod = Trajectory(:,Index);
            x_Mod = x_Mod./x_Mod(1);
            x_Mod(1) = [];
            J = (sum((((x_Mod - Data.TE_x)./Data.TE_sem).^2)));

        end
        
        function [Steady_State,SS] = ode_ss(Par,func,Ig,Index)
            
            options = optimset('Display','off');
            func = @(X) func(0,X,Par);
            SS = fsolve(func,Ig,options);
            Steady_State = SS(Index);
            
            
        end
        
       
        
        
        function [J,SS] = obj_ss(Par,obj,U,Index,Ut,Data,Ig,func)
           
           tf = strcmp(Ut,'Ca');
           
           if tf == 1
               
                for i = 1:length(Data.tit_x)
                    U.Ca = @(t) Data.tit_x(i);
                    [SS(i,:)] = obj.ode_ss(Par,func,Ig,Index);
                end
            
           end
           
           if tf == 0
                
                for i = 1:length(Data.tit_x)
                    U.HCO3 = @(t) Data.tit_x(i);
                    [SS(i,:)] = obj.ode_ss(Par,func,Ig,Index);
                end
           end
            
             J = sum(sum(((SS - Data.tit_y)./Data.tit_sem).^2));
        end
        
        
        
        
        function x = MN(K,U,t)
            
            Ca = U.Ca(t);
            HCO3 = U.HCO3(t);
            
            A = K(1)*K(2)*(HCO3*(1 + K(3)*Ca))/(K(4) + HCO3*(1 + K(3)*Ca));
            B = K(2);
            C = K(5)*Ca/(K(6) + Ca);
            D = 2*K(7);
            E = K(8);
            F = K(9);
            G = K(10)/2;
            H = K(11) - K(8);
            
           
            L = [ -1,  0,           log(A) - log(B*C)
                  4, -3,     log(D*E) - log((F*G)/2)
                 -1,  0,           log(A) - log(B*C)
                  4, -2,       log(D*E) - log(F*G*H)
                 -1,  0,           log(A) - log(B*C)
                  0, -3,     log(D*E) - log((E*G)/2)
                 -1,  0,           log(A) - log(B*C)
                  0, -2,       log(D*E) - log(E*G*H)
                 -1,  0,           log(A) - log(B*C)
                  0, -1,         log(D*E) - log(D/2)
                 -1,  0,           log(A) - log(B*C)
                  4,  1,     log(G/4) - log((F*G)/2)
                 -1,  0,           log(A) - log(B*C)
                  4,  2,       log(G/4) - log(F*G*H)
                 -1,  0,           log(A) - log(B*C)
                  0,  1,     log(G/4) - log((E*G)/2)
                 -1,  0,           log(A) - log(B*C)
                  0,  2,       log(G/4) - log(E*G*H)
                 -1,  0,           log(A) - log(B*C)
                  0,  3,         log(G/4) - log(D/2)
                 -1,  0,           log(A) - log(B*C)
                  4,  0, log((G*H)/2) - log((F*G)/2)
                 -1,  0,           log(A) - log(B*C)
                  4,  1,   log((G*H)/2) - log(F*G*H)
                 -1,  0,           log(A) - log(B*C)
                  0,  0, log((G*H)/2) - log((E*G)/2)
                 -1,  0,           log(A) - log(B*C)
                  0,  1,   log((G*H)/2) - log(E*G*H)
                 -1,  0,           log(A) - log(B*C)
                  0,  2,     log((G*H)/2) - log(D/2)
                 -1, -1,             log(A) - log(C)
                  4, -3,     log(D*E) - log((F*G)/2)
                 -1, -1,             log(A) - log(C)
                  4, -2,       log(D*E) - log(F*G*H)
                 -1, -1,             log(A) - log(C)
                  0, -3,     log(D*E) - log((E*G)/2)
                 -1, -1,             log(A) - log(C)
                  0, -2,       log(D*E) - log(E*G*H)
                 -1, -1,             log(A) - log(C)
                  0, -1,         log(D*E) - log(D/2)
                 -1, -1,             log(A) - log(C)
                  4,  1,     log(G/4) - log((F*G)/2)
                 -1, -1,             log(A) - log(C)
                  4,  2,       log(G/4) - log(F*G*H)
                 -1, -1,             log(A) - log(C)
                  0,  1,     log(G/4) - log((E*G)/2)
                 -1, -1,             log(A) - log(C)
                  0,  2,       log(G/4) - log(E*G*H)
                 -1, -1,             log(A) - log(C)
                  0,  3,         log(G/4) - log(D/2)
                 -1, -1,             log(A) - log(C)
                  4,  0, log((G*H)/2) - log((F*G)/2)
                 -1, -1,             log(A) - log(C)
                  4,  1,   log((G*H)/2) - log(F*G*H)
                 -1, -1,             log(A) - log(C)
                  0,  0, log((G*H)/2) - log((E*G)/2)
                 -1, -1,             log(A) - log(C)
                  0,  1,   log((G*H)/2) - log(E*G*H)
                 -1, -1,             log(A) - log(C)
                  0,  2,     log((G*H)/2) - log(D/2)];

                k = 1;
                p = 1;
                for i = 1:30
                    M = L([p;p+1],[1 2]);
                    V = L([p;p+1],[3]);
                    x(i,:) = M\V;
                    k = k + 1;
                    p = k*2 - 1;
                end

                x = exp(x);  
            
        end
        
        
        
% % % % % % % % % %         function Steady_State = ode_SS(Par,func,obj,Cpt,U,Xo,I,f) % Solve for steady state of system dy/dt
% % % % % % % % % %             
% % % % % % % % % %             e = @(t,X) obj.diff_event(t,X,Par,obj,Cpt,U,f);
% % % % % % % % % %              opts = odeset('Events', e);
% % % % % % % % % %          
% % % % % % % % % %             func = @(t,X) func(t,X,Par);
% % % % % % % % % %              [t,X] = ode45(func,[0 1000],Xo,opts);
% % % % % % % % % %            % [t,X] = ode45(func,[0 10000],Xo);
% % % % % % % % % %             Steady_State = X(end,I);
% % % % % % % % % %         end
% % % % % % % % %        
% % % % % % % % % %         function J = obj_ss(Par,obj,Cpt,U,Xo,I,Ut,Data,f,func)
% % % % % % % % % %            
% % % % % % % % % %            tf = strcmp(Ut,'Ca');
% % % % % % % % % %            
% % % % % % % % % %            if tf == 1
% % % % % % % % % %                
% % % % % % % % % %                 for i = 1:length(Data.tit_x)
% % % % % % % % % %                     U.Ca = @(t) Data.tit_x(i);
% % % % % % % % % %                     SS(i,:) = obj.ode_SS(Par,func,obj,Cpt,U,Xo,I,f);
% % % % % % % % % %                 end
% % % % % % % % % %             
% % % % % % % % % %            end
% % % % % % % % % %            
% % % % % % % % % %            if tf == 0
% % % % % % % % % %                 
% % % % % % % % % %                 for i = 1:length(Data.tit_x)
% % % % % % % % % %                     U.HCO3 = @(t) Data.tit_x(i);
% % % % % % % % % %                     SS(i,:) = obj.ode_SS(Par,func,obj,Cpt,U,Xo,I,f);
% % % % % % % % % %                 end
% % % % % % % % % %            end
% % % % % % % % % %             
% % % % % % % % % %              J = sum(sum(((SS - Data.tit_y)./Data.tit_sem).^2));
% % % % % % % % % %         end
        
        
        function [x,isterm,dir] = diff_event(t,X,K,obj,Cpt,U,f) % Define termination event at dy/dt = 0
            
            %dy = obj.fun(t,X,K,obj,Cpt,U);
            dy = f(t,X,K,obj,Cpt,U);
            x = norm(dy) - 1e-9 ;
            isterm = 1;
            dir = 0;
        
        end

    end
        
        

        
            
end
            

    

    
