function Mod = Get_Model(i)

switch i
    case 1
        Mod = SC_Model;
        
        % K(1) = k1;
        % K(2) = m3;
        % K(3) = c;
        % K(4) = m1;
        % K(5) = k_1;
        % K(6) = kc;
        % K(7) = k2;
        % K(8) = p1;
        % K(9) = m2;
        % K(10) = k_2;
        % K(11) = p2;
        
        % X(1) = cAMP;
        % X(2) = PKAc;
        % X(5) = pY;
        k1 = @(obj,X,t,K,U) K(1)*K(2)/(K(2) + X(2)); 
        
        Mod.Cpt(1).Input = @(obj,X,t,K,U) k1(obj,X,t,K,U)*U.HCO3(t)*(1 + K(3)*U.Ca(t))/(K(4) + U.HCO3(t)*(1 + K(3)*U.Ca(t)));
        Mod.Cpt(1).Output = @(obj,X,t,K,U) K(5)*(U.Ca(t))*X(1)/(K(6) + U.Ca(t));


        Mod.Cpt(2).Input = @(obj,X,t,K,U) 2*K(7)*X(1)^4*(K(8) - X(2)/2)/(K(9) + X(1)^4*(K(8) - X(2)/2));
        Mod.Cpt(2).Output = @(obj,X,t,K,U) (K(10)*X(2)^2/2) * ((K(11) - K(8)) + X(2)/2);
        
       
        
        Mod.Cpt(1).Initial = 1;
        Mod.Cpt(2).Initial = 1;
        
  
     
        Mod.Par = [1;1;1;1;1;1;1;1;1;1;1];
        Mod.U.HCO3 = @(t) 15;
        Mod.U.Ca = @(t) 1.7;
        
    
      
        
end
end

