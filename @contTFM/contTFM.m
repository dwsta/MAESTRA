classdef contTFM

methods (Static = true)
    
    TFM_step_nPass(rootdir,alias,cfg_data,nPass);

    [tx,ty,nxx,nxy,nyy,ut0,vt0,X0,Y0,Lx,Ly]...
    = mytrac_msm_2D_rik_from_dp(E,h,h0,nx0,ny0,x,y,u,v,fcalx,Nyq)

    [tauxz,tauyz,tauzz,u,v,w] = forces3D(alpha,beta,h,h0,sigma,uph0,vph0,wph0);
end





end