function [ aaa, bbb ] = removeUnwantedICs( keeplist, PCAquant, niifiles, validVoxels, X)
%removeUnwantedICs - Function to remove the effect of the user-specified
% ICs that are not of interest. tcgicaHINT is run immediately 
% after this function to obtain a new set of initial parameter estimates
% for the EM algorithm.
%
%Syntax:  [ XXX ] =
%removeUnwantedICs( XXX)
%
%Inputs:
%    keeplist - boolean list of the ICs to include in the analysis
%    A - ICA level mixing matrix from the group ica decomposition
%    stage1White - the PCA whitening matrix WHICH 
% .  stage2White
%    stage1dewhite - 
% .  stage2dewhite
%

%Outputs:
%    XXX    - Object containing initial estimates for the EM algorithm
%    
%
%See also: tcgicaHINT.m

% Get the list of ICs to remove

aaa=1; bbb=2;

%% Recover the data sizes used to obtain the initial guess;
nBlocks = PCAquant.nBlocks;
subjPerBlock = PCAquant.subjPerBlock;
nSubICBlocks = PCAquant.nSubICBlocks;
voxPerBlock = PCAquant.voxPerBlock;
numPCA = size(PCAquant.stage1dewhite,2);
numPCA2 = size(PCAquant.stage2dewhite,2);


%% Other needed quantities
V = sum(voxPerBlock);
N = length(niifiles);
qstar = sum(keeplist);
q = length(keeplist);
T = size(PCAquant.redICYWhite, 2);

%% Storage for things to keep
YtildeStar = zeros(qstar, V);
CmatStar = zeros(N * qstar, 1);

%% Storage for refining initial guess - these have the same names as the 
% whiting and dewhitening matrix from tcgicaHINT, however they are
% different. The matrices from the tcgicaHINT are stored under the PCAquant
% object
stage1White = zeros(numPCA, T, N); stage1deWhite = zeros(T, numPCA, N);
stage2White = zeros(numPCA2, numPCA*max(subjPerBlock), nBlocks); stage2deWhite = zeros(numPCA*max(subjPerBlock), numPCA2, nBlocks);
% Matrices used in the original data cleaning, need them again when
% reconstructing cleaned data to estimate S0, beta, sig2sq
HDGAiStore = zeros(T, q, N);
HDGAiremStore = zeros(T, q, N);
% Storage for second stage reduction
stackedStage2Data = zeros(numPCA2 * nBlocks, V);

%% Step 1 - Back reconstruct a "removal dataset" of the original subject data from the dewhitening
% matrices using the columns of the mixing matrix corresponding to ICs to
% be removed. Then load the original subject data and subtract off this
% removal dataset.
disp('Removing unwanted components from the data...')

% Create a mixing matrix with only the UNWANTED columns
Aremove = PCAquant.Agrp;
Aremove(:, find(keeplist == 1)) = 0;

% Back reconstruct a "unwanted" dataset, Yremove
subjIndex = 0;
glength = size(PCAquant.G, 1) / nBlocks;
for iGroup = 1:nBlocks
        stackedPCAData = zeros(numPCA * subjPerBlock(iGroup), V);
        gstart = glength*(iGroup-1) + 1; gend = glength*iGroup; % indices of G
        istart = -numPCA + 1; iend = 0;
        for iSubj = 1:subjPerBlock(iGroup)
            % adjust the storage indices
            istart = istart + numPCA; iend = iend + numPCA;
            subjIndex = subjIndex + 1;
            % Load the subjects data
            image = load_nii(niifiles{subjIndex}); res = reshape(image.img,[], T)';
            % Back reconstruct a "removal" dataset
            X_tilde_all = res(:,validVoxels); [X_tilde_all, ] = remmean(X_tilde_all);
            sistart = (iSubj-1) * numPCA + 1; siend = iSubj * numPCA ;
            HDGAi =  squeeze(PCAquant.stage1dewhite(:,:,subjIndex)) * ...
                     squeeze(PCAquant.stage2dewhite(sistart:siend,:,iGroup)) *...
                     PCAquant.G(gstart:gend, :) * PCAquant.Agrp;
            HDGAiStore(:, :, subjIndex) = HDGAi;
            subIC = pinv(HDGAi) * X_tilde_all;
            % Create the transformation matrix with columns of A
            % corresponding to ICs to keep in the analysis removed.
            HDGAirem = squeeze(PCAquant.stage1dewhite(:,:,subjIndex)) * ...
                     squeeze(PCAquant.stage2dewhite(sistart:siend,:,iGroup)) *...
                     PCAquant.G(gstart:gend, :) * Aremove;
            HDGAiremStore(:, :, subjIndex) = HDGAirem;
            X_tilde_remove = HDGAirem * subIC;
            % Clean the dataset by removing the unwanted ICs
            X_tilde_cleaned = X_tilde_all - X_tilde_remove;
            
            % Obtain YtildeStar, the reduced subject level data based on
            % the reduced number of ICs
            [U_incr, D_incr] = pcamat(X_tilde_cleaned);
            lambda = sort(diag(D_incr),'descend');
            sigma2_ML = sum(lambda(qstar+1:length(lambda))) / (length(lambda)-qstar);   
            U_q = U_incr(:,(size(U_incr,2)-qstar+1):size(U_incr,2));
            D_q = diag(D_incr((size(U_incr,2)-qstar+1):size(U_incr,2),...
                (size(U_incr,2)-qstar+1):size(U_incr,2)));
            whitenMat = diag((D_q-sigma2_ML).^(-1/2)) * U_q';
            dewhitenMat = U_q * diag((D_q-sigma2_ML) .^ (1/2));
            YtildeStar = whitenMat * X_tilde_cleaned;
            % Store the corresponding part of C_matrix_diag
            cIndStart = 1+(subjIndex-1)*qstar; cIndEnd = (subjIndex)*qstar;
            CmatStar(cIndStart:cIndEnd,:) =...
                diag( ( whitenMat * (eye(T)-1/T * ones(T)) ) *...
                ( whitenMat * (eye(T)-1/T * ones(T)) )' );
            
            % Now we need to perform the same types of PCA reduction on the
            % cleaned data X_tilde_cleaned
            lambda = sort(diag(D_incr),'descend');
            sigma2_ML = sum(lambda(numPCA+1:length(lambda))) / (length(lambda)-numPCA);   
            U_q = U_incr(:,(size(U_incr,2)-numPCA+1):size(U_incr,2));
            D_q = diag(D_incr((size(U_incr,2)-numPCA+1):size(U_incr,2),...
                (size(U_incr,2)-numPCA+1):size(U_incr,2)));
            % Store the important pieces
            stage1White(:, :, subjIndex) = diag((D_q-sigma2_ML).^(-1/2)) * U_q';
            stage1deWhite(:, :, subjIndex) = U_q * diag((D_q-sigma2_ML) .^ (1/2));
            stackedPCAData(istart:iend, :) = stage1White(:, :, subjIndex) * X_tilde_all;           
        end
        
        %% Stage 2 - Perform the second stage of reduction on the block of data
        if nBlocks > 1
            [U_incr_2, D_incr_2] = pcamat(stackedPCAData);
            lambda = sort(diag(D_incr_2),'descend');
            if numPCA2 < numPCA
               sigma2_ML = sum(lambda(numPCA2+1:length(lambda))) / (length(lambda)-numPCA2);   
            else
               sigma2_ML = 0;
            end
            U_q = U_incr_2(:,(size(U_incr_2,2)-numPCA2+1):size(U_incr_2,2));
            D_q = diag(D_incr_2((size(U_incr_2,2)-numPCA2+1):size(U_incr_2,2),...
               (size(U_incr_2,2)-numPCA2+1):size(U_incr_2,2)));

            ng = numPCA * subjPerBlock(iGroup); % Account for possible different group size
            stage2White(:,1:ng,iGroup) = diag((D_q-sigma2_ML).^(-1/2)) * U_q';
            stage2deWhite(1:ng,:,iGroup) = U_q * diag((D_q-sigma2_ML) .^ (1/2));

            % Store the 2nd stage results
            storeStart = (iGroup-1)*numPCA2 + 1;
            storeEnd = numPCA2 * iGroup;
            stackedStage2Data(storeStart:storeEnd, :) = stage2White(:,1:ng,iGroup) * stackedPCAData;
        else % handle typical case where blocks not required
            stage2White = eye(numPCA2);
            stage2deWhite = eye(numPCA2);
            stackedStage2Data = stackedPCAData;
        end
end

%% Stage 3 - Perform PCA down to q components
[U_incr, D_incr] = pcamat(stackedStage2Data);
lambda = sort(diag(D_incr),'descend');
sigma2_ML = sum(lambda(qstar+1:length(lambda))) / (length(lambda)-qstar);   
% Recover the whitening matrix (G inverse)
U_q = U_incr(:,(size(U_incr,2)-qstar+1):size(U_incr,2));
D_q = diag(D_incr((size(U_incr,2)-qstar+1):size(U_incr,2),...
        (size(U_incr,2)-qstar+1):size(U_incr,2)));
Ginv = diag((D_q).^(-1/2)) * U_q';
groupPCA_X = Ginv * stackedStage2Data; % this is X that the ICA is performed on

%% Stage 4 - ICA
[icasig, Amix, W] = fastica_full(groupPCA_X, 'numOFIC', qstar);
Agrp = Amix;

%% Reconstruct the single subject matrices - Over Blocks of Voxels

H = stage2deWhite;
G = pinv(Ginv);

%%%%%%%
for iVoxBlock = 1:nSubICBlocks 
    % Create storage for this round of subject level ICs
    singleSubjectICs = zeros(qstar, voxPerBlock(iVoxBlock), N);
    voxStart = 1;
    if iVoxBlock > 1
        voxStart = 1 + sum(voxPerBlock(1:(iVoxBlock-1)));
    end
    voxEnd = sum(voxPerBlock(1:iVoxBlock));
    disp(['Estimating for voxels ' num2str(voxStart) ' through ' num2str(voxEnd)])

    % Loop through and create the subject level ICs for the relevant voxels
    glength = size(G, 1) / nBlocks;
    subjIndex = 0; % track the current subject
    for iGroup = 1:nBlocks
        gstart = glength*(iGroup-1) + 1; gend = glength*iGroup; % indices of G
        for iSubj = 1:subjPerBlock(iGroup)
            subjIndex = subjIndex + 1;
            % index of the dewhiten matrix for this subj
            image = load_nii(niifiles{subjIndex}); res = reshape(image.img,[], T)';
            X_tilde_all = res(:,validVoxels); [X_tilde_all, ] = remmean(X_tilde_all);
            disp('store pinv instead so dont have to do this twice')
            subIC =  pinv(HDGAiStore(:, :, subjIndex)) * X_tilde_all;
            X_tilde_remove = HDGAiremStore(:, :, subjIndex) * subIC;
            X_tilde_cleaned = X_tilde_all - X_tilde_remove;
            % Create the subject-specific matrix based on the NEW whitening
            % mats
            sistart = (iSubj-1) * numPCA + 1; siend = iSubj*numPCA ;
            HDGAi =  squeeze(stage1deWhite(:,:,subjIndex)) * ...
                squeeze(stage2deWhite(sistart:siend,:,iGroup)) *...
                G(gstart:gend, :) * Agrp;
            singleSubjectICs(:,:,subjIndex) = pinv(HDGAi) * X_tilde_cleaned(:, voxStart:voxEnd);
        end
    end
    
    disp('remove scaling')
    singleSubjectICs = singleSubjectICs ./ max(abs(singleSubjectICs(:)));
    
    % Create a design matrix that includes an intercept (this is S0)
    Xfull = [ ones(N, 1), X];
    
    % Start getting estimates
    epsilon2temp = zeros(qstar, V, N);
    voxRange = voxStart:voxEnd;
    for v = 1:voxPerBlock(iVoxBlock)
        estimate = (Xfull'*Xfull)^(-1) * Xfull' * squeeze(singleSubjectICs(:,v,:))';
        beta(:,:,voxRange(v)) = estimate(2:(p), :);
        S0(:,voxRange(v)) = estimate(1,:)';
        aggregate_IC_map(:, v) = mean(Xfull * estimate); % do this here to avoid reconstruction
        epsilon2temp(:,voxRange(v),:) = singleSubjectICs(:,v,:) - reshape( (Xfull * estimate)', [q,1,N] );
    end
    sigma2_sq = var(reshape(epsilon2temp, [q,V*N]), 0, 2);
end    

xxx=1;
end

