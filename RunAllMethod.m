clear all
close all
addpath ./Quality_Assess
load Noisy_indian_CASE4
%subspace = 4; %Salinas
 subspace = 8; %simu_indian
k_subspace = subspace;
X = Img;
Y = Noisy_Img;
 bd = 110; % 显示第bd光谱段图像 
 figure
 Singleband = X(:,:,bd);
 imshow(Singleband);
 figure
 Singleband = Y(:,:,bd);
 imshow(Singleband);

[H,W,B]  = size(Noisy_Img);
clear Img Noisy_Img;
[PSNR(:,1), SSIM(:,1), ERGAS(1,1), SAM(1,1)] = MSIQA(X*255, Y*255);
 fprintf( 'Noisy Image: nSig = %2.2f, PSNR = %2.2f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,1)),mean(SSIM(:,1)),ERGAS(1,1),SAM(1,1));
%% ========== METHOD SWITCH==========
STATUS_LRMR  = 0; 
STATUS_BM4D  = 0;
STATUS_LRTDTV =0;
STATUS_Enhanced3DTV = 0;
STATUS_NMoG = 0;
 STATUS_GLF =0;
STATUS_SLRL4D = 0;
STATUS_L1HyMixDe = 0;
STATUS_FGSLR = 0;
STATUS_RCTV = 0;
STATUS_WMLRTR = 0;
STATUS_ETPTV = 0;
STATUS_My =1;


%% LRMR
if STATUS_LRMR
    addpath (genpath('./EH_LRMR'))
    r = 7;
    slide =20;
    s = 0.15;
    stepsize = 8;
%     fprintf('\n');
%     disp(['performing ',METHODSET{index}, ' ... ']);
    tic;
    E_Img    = LRMR_HSI_denoise( Y,r,slide,s,stepsize );
   toc;
[PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
fprintf( 'Denoised Image(LRMR): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
     rmpath('./EH_LRMR');
end

%% BM4D
if STATUS_BM4D
    addpath(genpath('./EH_BM4D'))
%     sigma=0.1;
%     fprintf('\n');
%     disp(['performing ',METHODSET{index}, ' ... ']);
    tic
%     if gausssigma~=0  %有Gaussian
         if true
        [~, E_Img,sigma_est] = bm4d(X, Y, sigma ,0);
    end
 toc;
    [PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
fprintf( 'Denoised Image(BM4D): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
    rmpath('./EH_BM4D');
end


%% LRTDTV
if STATUS_LRTDTV
    addpath(genpath('./EH_LRTDTV'))
%     disp(['performing ',METHODSET{index}, ' ... ']);
gausssigma = sigma;

        tau    = 1;
%        beta=100;
        lambda =10;
        Rank   = [round(H*0.8),round(W*0.8),10];
        fprintf('\n');
        tic
%        E_Img               = LRTDTV_G(Y, tau,lambda,beta,Rank);
     E_Img = LRTDTV(Y, tau,lambda,Rank);
         toc;
           [PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
 fprintf( 'Denoised Image(LRTDTV): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));


 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
    rmpath('./EH_LRTDTV');
 end


%% Enhanced3DTV
if STATUS_Enhanced3DTV
    addpath(genpath('./EH_Enhanced3DTV'))
    rank   = [13,13,13];
    tau    = 0.004 *sqrt(H*W);
%     fprintf('\n');
%     disp(['performing ',METHODSET{index}, ' ... ']);
    tic
    E_Img = EnhancedTV(Y,tau,rank);
 toc;
    [PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
fprintf( 'Denoised Image(Enhanced3DTV): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
   rmpath('./EH_Enhanced3DTV');
end
%% STATUS_NMoG
if STATUS_NMoG
    addpath(genpath('./NG_NMoG_RPCA'))
     muOn = 0;                     % muOn = 0: set mu as 0 without updating;
Rank = 5;                     % objective rank of low rank component
param.initial_rank = 30;      % initial rank of low rank component
param.rankDeRate = 7;         % the number of rank reduced in each iteration
param.mog_k = 3;              % the number of component reduced in each band
param.lr_init = 'SVD';
param.maxiter = 30;
param.tol = 1e-4;
param.display = 0;
     tic
    [prior, model] = InitialPara(param,muOn,B);  % set hyperparameters and initialize model parameters
    N_I = reshape(Y,H*W,B);
    [Model,Lr_model] = NMoG_RPCA( N_I,Rank,param,model,prior);
    U = Lr_model.U;
    V = Lr_model.V;
    E_Img = reshape(U*V',size(Y));
     toc;
    [PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
fprintf( 'Denoised Image(NMoG): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
  rmpath('./NG_NMoG_RPCA');
end
%% GLF
if STATUS_GLF
 
    addpath(genpath('./HSI-denoiser-GLF-master'))
    noise_type = 'additive';
    tic
    E_Img = GLF_denoiser(Y, subspace,  noise_type) ;
    toc
    [PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
fprintf( 'Denoised Image(GLF): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
  rmpath('./HSI-denoiser-GLF-master');
end
%% STATUS_SLRL4D 
if STATUS_SLRL4D
    addpath(genpath('./NG_SLRL4D'))
    k_num=4;
    lambda=0.03;
iterations = 100;
tic
[E_Img,~,~,~,~] =SLRL4D_fast(Y,'LAMBDA_S',lambda,'SUBSPACE_DIM',k_num,'AL_ITERS',iterations,'TRUE_X',X);
toc
    [PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
fprintf( 'Denoised Image(SLRL4D): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
  rmpath('./NG_SLRL4D');
end
%% STATUS_L1HyMixDe
if STATUS_L1HyMixDe
     addpath(genpath('./L1HyMixDe'))
      addpath(genpath('./BM3D'))
    p = 0.1;
     tic
         E_Img = L1HyMixDe(Y, k_subspace, p);
     toc
    [PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
fprintf( 'Denoised Image(L1HyMixDe): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
  rmpath('./L1HyMixDe');
  rmpath('./BM3D');
end
%% FGSLR 
if STATUS_FGSLR 
        addpath(genpath('./FGSLR'))
              opt.alpha = 1;
      opt.beta = 0.1;          % TV regularization:set as 0.1 or 0.5  
      opt.lambda = 0.01;       % sparse noise
      opt.mu = 5;              % penalty parameter:set as 5 or 10
      opt.rho = 0.1;           % proximal parameter
      opt.delta = 0.5;         % noise parameter:set as 0.5 or 5
      opt.r = 20;
      opt.regul_B = 'L21';
      [Q, ~, ~] = FGSLR_TV_PAM(Y, X, opt);
     E_Img = reshape(Q',H,W,B);
          toc
    [PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
fprintf( 'Denoised Image(FGSLR): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
  rmpath('./FGSLR');
end
%% RCTV
if STATUS_RCTV
    addpath(genpath('./RCTV')) 
r=k_subspace;
beta = 50;
lambda = 1;% 5,0.5
tau = [0.8,0.8];% need to fine tune
tic;
E_Img = RCTV(Y, beta,lambda, tau, r);
 toc;
    [PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
fprintf( 'Denoised Image(RCTV): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
  rmpath('./RCTV');
end
%% SMLRTR 
if STATUS_WMLRTR
  addpath(genpath('./WMLRTR')) 
 
  Par = ParSet(sigma, subspace)
 tic
 [E_Img] = HSI_Denoising(Y, Par);
 toc;
    [PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
fprintf( 'Denoised Image(SMLRTR): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
  rmpath('./WMLRTR');
end
%% ETPTV
if STATUS_ETPTV
     addpath(genpath('./ETPTV')) 
    j=10;    
    param.Rank   = [7,7,5];
    param.initial_rank = 2;
     param.maxIter = 50;
        param.lambda    = 4e-3*sqrt(H*W);  %0.1;   %mu2      = Alpha*mu1
        tic
        [  E_Img,U_x,V_x,E] = WETV(Y,X, param);
       toc;
    [PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
fprintf( 'Denoised Image(ETPTV): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband);
 rmpath('./ETPTV');
end
%% My_Method
if STATUS_My
    addpath ./nonconvex
    addpath ./Utilize

    Par = ParSet(sigma, subspace);
tic
[E_Img] = My_denoise(Y, Par);
toc

E_Img(E_Img>1)=1;
E_Img(E_Img<0)=0;

[PSNR(:,2), SSIM(:,2), ERGAS(1,2), SAM(1,2)] = MSIQA(X*255, E_Img*255);
 fprintf( 'Denoised Image(Our): nSig = %2.2f, PSNR = %2.4f, SSIM = %2.4f, ERGAS = %2.4f,SAM = %2.4f\n', sigma, mean(PSNR(:,2)),mean(SSIM(:,2)),ERGAS(1,2),SAM(1,2));
 figure
 Singleband = E_Img(:,:,bd);
 imshow(Singleband); 
 rmpath('./nonconvex');
 rmpath('./Utilize');
end