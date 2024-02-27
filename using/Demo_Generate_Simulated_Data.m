%% Demo on generate simulate HSI data
clc; 
clear all;
%% Please select dataset and case
dataset = 0; % 0 for Salinas, 1 for Indian
CaseNum = 4; % case number, 1 to 4
%% Load HSI
if dataset == 0
  load('Salinas.mat');
  Img = salinas;                    % not normalize
elseif dataset == 1
  load('simu_indian.mat');          % not normalize
   Img = simu_indian;
end

%% nomalize the image
[nr,nc,L] = size(Img);
%nomorlize the observed X
Img = Img./repmat(max(max(Img,[],1),[],2),nr,nc);

%% Add noise in different ways
Noisy_Img = Img;
CaseType  = 4 * dataset + CaseNum;
switch CaseType
    case 1
        % zero-mean Gaussian noise with the different standard deviationRL_sp
        % randomly sampled from [0.1,0.2] was added to the different band
        sig =  0.1+(0.2-0.1)*rand(L,1);
        for i = 1 : L
            Noisy_Img(:,:,i) = Img(:,:,i) + sig(i)*randn(nr,nc);
        end
        sigma = mean(sig);
        save Noisy_Salinas_CASE1 Img Noisy_Img sigma;
    case 2
        % Gaussian noise the same as case 1, and impulse noise with a
        % density percentage of 20% was added 20 bands randomly
        sig =  0.1+(0.2-0.1)*rand(L,1);
        for i = 1 : L
            Noisy_Img(:,:,i) = Img(:,:,i) + sig(i)*randn(nr,nc);
        end
        sigma = mean(sig);
        ratio = 0.2;
        RL_sp = randperm(L,20);
        for i = 1 : 20
            Noisy_Img(:,:,RL_sp(i)) = imnoise(Noisy_Img(:,:,RL_sp(i)), 'salt & pepper', ratio);
        end
        save Noisy_Salinas_CASE2 Img Noisy_Img sigma ratio;
    case 3
        %  Gaussian noise and impulse noise the same as case 2
        %  dead lines were added to 20% band randomly, and the width
        %  of each dead line was randomly generated from 1 to 3
        sig =  0.1+(0.2-0.1)*rand(L,1);
        for i = 1 : L
            Noisy_Img(:,:,i) = Img(:,:,i) + sig(i)*randn(nr,nc);
        end
        sigma = mean(sig);
        ratio = 0.2;
        RL_sp = randperm(L,20);
        for i = 1:20
            Noisy_Img(:,:,RL_sp(i)) = imnoise(Noisy_Img(:,:,RL_sp(i)), 'salt & pepper', ratio);
        end
        % of deadlines being randomly selected and the width
        % of each dead line was randomly generated from 1 to 3
        RL_dl = randperm(L,16);
        for i=1:16
            indp=randperm(8,1)+2;
            ind=randperm(nc-1,indp);
            an = zeros(length(ind),1);
            for k=1:length(ind)
                an(k)=randperm(3,1);
            end
            % searching the location of an which value is 1,2,3
            loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
            Noisy_Img(:,ind(loc1),RL_dl(i))=1; 
            Noisy_Img(:,ind(loc2):ind(loc2)+1,RL_dl(i))=1;
            Noisy_Img(:,ind(loc3)-1:ind(loc3)+1,RL_dl(i))=1;
        end
        % dead lines by another way
%         for i=126:145
%             loc_dl = ceil(rand(1,20)*nc);
%             loc_dl = [loc_dl, 20:22];
%             loc_dl = [loc_dl, 142:143];
%             Noisy_Img(:,loc_dl,i) = ones(nr,size(loc_dl,2));
%         end  
        save Noisy_Salinas_CASE3 Img Noisy_Img RL_sp RL_dl sigma;
    case 4
        %  same as case 3
        %  dead lines were added to 20% band randomly, and the width
        %  of each dead line was randomly generated from 1 to 3
        sig =  0.1+(0.2-0.1)*rand(L,1);
        for i = 1 : L
            Noisy_Img(:,:,i) = Img(:,:,i) + sig(i)*randn(nr,nc);
        end
        sigma = mean(sig);
        ratio = 0.2;
        RL_sp = randperm(L,20);
        for i = 1:20
            Noisy_Img(:,:,RL_sp(i)) = imnoise(Noisy_Img(:,:,RL_sp(i)), 'salt & pepper', ratio);
        end
        % of deadlines being randomly selected and the width
        % of each dead line was randomly generated from 1 to 3
        RL_dl = randperm(L,16);
        for i=1:16
            indp=randperm(8,1)+2;
            ind=randperm(nc-1,indp);
            an = zeros(length(ind),1);
            for k=1:length(ind)
                an(k)=randperm(3,1);
            end
            % searching the location of an which value is 1,2,3
            loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
            Noisy_Img(:,ind(loc1),RL_dl(i))=1; 
            Noisy_Img(:,ind(loc2):ind(loc2)+1,RL_dl(i))=1;
            Noisy_Img(:,ind(loc3)-1:ind(loc3)+1,RL_dl(i))=1;
        end
        % dead lines by another way
%         for i=126:145
%             loc_dl = ceil(rand(1,20)*nc);
%             loc_dl = [loc_dl, 20:22];
%             loc_dl = [loc_dl, 142:143];
%             Noisy_Img(:,loc_dl,i) = ones(nr,size(loc_dl,2));
%         end  
        % stripes: 40% band and each band was randomly generated from 6 to 15
        RL_s = randperm(L, 32);
        for i = 1:length(RL_s)
            stripenum = randperm(10,1)+5;
            locolumn    = randperm(nc,stripenum);
            Noisy_Img(:,locolumn,RL_s(i))=0.2*rand(1)+0.6;
        end
        save Noisy_Salinas_CASE4 Img Noisy_Img RL_sp RL_dl RL_s sigma;

    case 5
        % zero-mean Gaussian noise with the different standard deviation
        % randomly sampled from [0.1,0.2] was added to the different band
        sig =  0.1+(0.2-0.1)*rand(L,1);
        for i = 1 : L
            Noisy_Img(:,:,i) = Img(:,:,i) + sig(i)*randn(nr,nc);
        end
        sigma = mean(sig);
        save Noisy_indian_CASE1 Img Noisy_Img sigma;
    case 6
        % Gaussian noise the same as case 5, and impulse noise with a
        % density percentage of 20% was added 20 bands randomly
        sig =  0.1+(0.2-0.1)*rand(L,1);
        for i = 1 : L
            Noisy_Img(:,:,i) = Img(:,:,i) + sig(i)*randn(nr,nc);
        end
        sigma = mean(sig);
        ratio = 0.2;
        RL_sp = randperm(L,20);
        for i = 1 : 20
            Noisy_Img(:,:,RL_sp(i)) = imnoise(Noisy_Img(:,:,RL_sp(i)), 'salt & pepper', ratio);
        end
        save Noisy_indian_CASE2 Img Noisy_Img sigma ratio;
    case 7
        %  Gaussian noise and impulse noise the same as case 3
        %  dead lines were added to 20% band, and the width
        %  of each dead line was randomly generated from 1 to 3
        sig =  0.1+(0.2-0.1)*rand(L,1);
        for i = 1 : L
            Noisy_Img(:,:,i) = Img(:,:,i) + sig(i)*randn(nr,nc);
        end
        sigma = mean(sig);
        ratio = 0.2;
        RL_sp = randperm(L,20);
        for i = 1:20
            Noisy_Img(:,:,RL_sp(i)) = imnoise(Noisy_Img(:,:,RL_sp(i)), 'salt & pepper', ratio);
        end
        % dead lines were added from band 126 to band 145 with the number
        % of deadlines being randomly selected from 3 to 10, and the width
        % of each dead line was randomly generated from 1 to 3
        RL_dl = randperm(L,38);
%         RL_dl = sort(RL_dl);
        for i=1:38
            indp=randperm(8,1)+2;
            ind=randperm(nc-1,indp);
            an = zeros(length(ind),1);
            for k=1:length(ind)
                an(k)=randperm(3,1);
            end
            % searching the location of an which value is 1,2,3
            loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
            Noisy_Img(:,ind(loc1),RL_dl(i))=1; 
            Noisy_Img(:,ind(loc2):ind(loc2)+1,RL_dl(i))=1;
            Noisy_Img(:,ind(loc3)-1:ind(loc3)+1,RL_dl(i))=1;
        end
        % dead lines by another way
%         for i=126:145
%             loc_dl = ceil(rand(1,20)*nc);
%             loc_dl = [loc_dl, 20:22];
%             loc_dl = [loc_dl, 142:143];
%             Noisy_Img(:,loc_dl,i) = ones(nr,size(loc_dl,2));
%         end  
        save Noisy_indian_CASE3 Img Noisy_Img RL_sp RL_dl sigma; 
        case 8
        %  same as case 3
        %  dead lines were added to 20% band randomly, and the width
        %  of each dead line was randomly generated from 1 to 3
        sig =  0.1+(0.2-0.1)*rand(L,1);
        for i = 1 : L
            Noisy_Img(:,:,i) = Img(:,:,i) + sig(i)*randn(nr,nc);
        end
        sigma = mean(sig);
        ratio = 0.2;
        RL_sp = randperm(L,20);
        for i = 1:20
            Noisy_Img(:,:,RL_sp(i)) = imnoise(Noisy_Img(:,:,RL_sp(i)), 'salt & pepper', ratio);
        end
        % of deadlines being randomly selected and the width
        % of each dead line was randomly generated from 1 to 3
        RL_dl = randperm(L,38);
        for i=1:38
            indp=randperm(8,1)+2;
            ind=randperm(nc-1,indp);
            an = zeros(length(ind),1);
            for k=1:length(ind)
                an(k)=randperm(3,1);
            end
            % searching the location of an which value is 1,2,3
            loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
            Noisy_Img(:,ind(loc1),RL_dl(i))=1; 
            Noisy_Img(:,ind(loc2):ind(loc2)+1,RL_dl(i))=1;
            Noisy_Img(:,ind(loc3)-1:ind(loc3)+1,RL_dl(i))=1;
        end
        % dead lines by another way
%         for i=126:145
%             loc_dl = ceil(rand(1,20)*nc);
%             loc_dl = [loc_dl, 20:22];
%             loc_dl = [loc_dl, 142:143];
%             Noisy_Img(:,loc_dl,i) = ones(nr,size(loc_dl,2));
%         end  
        % stripes: 40% band and each band was randomly generated from 6 to 15
        RL_s = randperm(L, 76);
        for i = 1:length(RL_s)
            stripenum = randperm(10,1)+5;
            locolumn    = randperm(nc,stripenum);
            Noisy_Img(:,locolumn,RL_s(i))=0.2*rand(1)+0.6;
        end
        save Noisy_indian_CASE4 Img Noisy_Img RL_sp RL_dl RL_s sigma;
end