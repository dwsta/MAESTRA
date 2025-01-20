function [LinESol, NHSol] = runOptimNHLinEWPades(UCF,VCF,UBF,VBF,X,Y,runLinE,runNH)
    
% Setup and run linear Elasticity (E,nu) model and/or Neo-Hookean model
    % (G,K) for monolayer
    % INPUTS
    % UCF,VCF - Filtered Cell displacements
    % UBF, VBF - Filtered Bead displacements
    % dV - Divergence field
    % X,Y - Positions for UF/UBF

    LinESol = struct;  NHSol = struct;
    
    for i = 1 : size(UCF,3)
        
        Ut = UCF(:,:,i); Vt = VCF(:,:,i);
        UBt = UBF(:,:,i); VBt = VBF(:,:,i);

        dx = X(1,2) - X(1,1);
        dy = Y(2,1) - Y(1,1);
        
        % 2D strain matrix (epsilon) 
        dUdX1 = dfdx_pade_2D_rik(Ut,dx); % Exx
        dUdX2 = dfdy_pade_2D_rik(Ut,dy);
        dVdX1 = dfdx_pade_2D_rik(Vt,dx);
        dVdX2 = dfdy_pade_2D_rik(Vt,dy); % Eyy
        Exy = (dUdX2 + dVdX1)./2;

        if runNH
            % F - deformation gradient tensor
            % J - Jacobi matrix (|F|)
            % shTerma and bulkTerm are multiples of G and K in th energy term
            % for NH model
            [~,~,shTerm,bulkTerm] = calcDeformationGradient(dUdX1,Exy,[],Exy,dVdX2,[],[],[],[]);
        end
        
        if runLinE
            % Derivtive of Strain matrix
            dE = struct();
            Exx = dUdX1; Eyy = dVdX2; Exy = (dUdX2+dVdX1)./2;
            
            dExxdX1 = dfdx_pade_2D_rik(Exx,dx);
            dExxdX2 = dfdy_pade_2D_rik(Exx,dy);
            
            dEyydX1 = dfdx_pade_2D_rik(Eyy,dx);
            dEyydX2 = dfdy_pade_2D_rik(Eyy,dy);
            
            dExydX1 = dfdx_pade_2D_rik(Exy,dx);
            dExydX2 = dfdy_pade_2D_rik(Exy,dy);
            
            dE.dExxdX1 = dExxdX1;
            dE.dExxdX2 = dExxdX2;
            dE.dEyydX1 = dEyydX1;
            dE.dEyydX2 = dEyydX2;
            dE.dExydX1 = dExydX1;
            dE.dExydX2 = dExydX2;
        end

        uBeads = [reshape(UBt,[],1); reshape(VBt,[],1)];
        
        if runLinE 
            [solE,exitFlagE,fvalE,uPredictE] = runoptimsolveElasto_Eandnu([],...
                                        dx,dE,uBeads,{X,Y},@objFuncDefElasto_Eandnu_reform1_FourierTractions);
            uuXE = reshape(uPredictE(1:end/2),size(X)); uuYE = reshape(uPredictE(1+end/2:end),size(X));
            LinESol(i).sol = solE;
            LinESol(i).fval = fvalE;
            LinESol(i).uuX = uuXE;
            LinESol(i).uuY = uuYE;
            LinESol(i).exitflag = exitFlagE;

            % Calculate the stress tensor
            
            sxx = (solE.E./(1-solE.nu.^2)) .* (Exx + solE.nu .* Eyy);
            syy = (solE.E./(1-solE.nu.^2)) .* (Eyy + solE.nu .* Exx);
            sxy = (solE.E./(1+solE.nu)) .* (Exy);
            
            LinESol(i).sxx = sxx;
            LinESol(i).syy = syy;
            LinESol(i).sxy = sxy;
        end 
        
        if runNH
            [sol,exitflag,fval,uPredict] = runoptimsolveNH_Pades([],...
                                    dx,dx,shTerm,bulkTerm,uBeads,{X,Y});

            uuX = reshape(uPredict(1:end/2),size(X)); uuY = reshape(uPredict(1+end/2:end),size(X));
        
            NHSol(i).sol = sol;
            NHSol(i).fval = fval;
            NHSol(i).uuX = uuX;
            NHSol(i).uuY = uuY;
            NHSol(i).exitflag = exitflag;
            
            % Calculate the stress tensor (First Piola-Kirchoff tensor)       
            Pxx = sol.G * squeeze(shTerm(:,:,1,1)) + sol.K * squeeze( shTerm(:,:,1,1) ) ;
            Pyy = sol.G * squeeze(shTerm(:,:,2,2)) + sol.K * squeeze( shTerm(:,:,2,2) ) ;
            Pxy = sol.G * squeeze(shTerm(:,:,1,2)) + sol.K * squeeze( shTerm(:,:,1,2) ) ;
            Pyx = sol.G * squeeze(shTerm(:,:,2,1)) + sol.K * squeeze( shTerm(:,:,2,1) ) ;
            NHSol(i).Pxx = Pxx;
            NHSol(i).Pyy = Pyy;
            NHSol(i).Pxy = 0.5*(Pxy + Pyx);
        end
        
    end   
end