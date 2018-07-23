function [ theta, beta, grpSig, s0_agg, PCAquant ] = tcgicaHINT_ext(niifiles, maskf, validVoxels,...
    X, numPCA, q, N, time_num, V, prefix)

%% Stages:
% 1. PCA dimension reduction of each subject down to nPC x V
% 2. PCA dimension reduction of each block of nPC x V datasets down to nPC2
%             block = 1, ..., nBlock
% 3. PCA dimension reduction of stacked block data from
%        nPC2*nBlock x V   to   q x V
% 4. ICA on the q x V dataset to get initial values

PCAquant = struct();
T = time_num;
p = size(X, 2) + 1;
%maxBytes = 1e9; disp('change max bytes to 8e9!!!')
maxBytes = 10e9; disp('change max bytes to 8e9!!!')

%% Preallocate for final things to return
beta = zeros(p-1, q, V);
S0 = zeros(q, V);
A = zeros(q,q,N);

%% Determine the required number of blocks based on the requested number of PCs
bytesRequired = 8 * numPCA * V * N;
nBlocks = ceil(bytesRequired / maxBytes);
%nBlocks = ceil(bytesRequired / 8e9);
subjPerBlock = floor(N / nBlocks) * ones(nBlocks, 1);
subjPerBlock(nBlocks) = N - (nBlocks-1)*floor(N / nBlocks);
% determine number of PCAs for the stacked data, goal is to be as large as
% possible -> so as close to subjPerBlock * numPCA as possible
maxBytesPerBlock = maxBytes ./ nBlocks;
numPCA2 = min( floor(maxBytesPerBlock ./ (8 * V)), floor(numPCA * subjPerBlock(1)) );

% Determine how large the stacked subject ica data would be (split on
% voxels)
bytesRequiredSubjICs = 8 * q * V * N;
disp('change this to 8e9')
%voxBytes = 0.2e9;
voxBytes = maxBytes;
nSubICBlocks = ceil(bytesRequiredSubjICs / voxBytes);
%nSubICBlocks = ceil(bytesRequiredSubjICs / 8e9);
voxPerBlock = floor(V / nSubICBlocks) * ones(nSubICBlocks, 1);
voxPerBlock(nSubICBlocks) = V - (nSubICBlocks-1)*floor(V / nSubICBlocks);

% Setup whitening and dewhiten storage for the intermediate steps
redICYWhite = zeros(q, T, N); redICYdewhite = zeros(T, q, N);
stage1White = zeros(numPCA, T, N); stage1deWhite = zeros(T, numPCA, N);
stage2White = zeros(numPCA2, numPCA*max(subjPerBlock), nBlocks); stage2deWhite = zeros(numPCA*max(subjPerBlock), numPCA2, nBlocks);
%stage3White = zeros(q, numPCA2*nBlocks); stage3dewhite = zeros(numPCA2*nBlocks, q);

stackedStage2Data = zeros(numPCA2 * nBlocks, V);

% Loop over each block and do the first and second stages of dimension reduction
for iGroup = 1:nBlocks

    %% Stage 1: PCA dimension reduction
    stackedPCAData = zeros(numPCA * subjPerBlock(iGroup), V);
    istart = -numPCA + 1; iend = 0; k = time_num;
    for i = 1:subjPerBlock(iGroup)
        subjInd = i + sum(subjPerBlock(1:(iGroup-1)));
        disp(['Stage 1 PCA for subject ' num2str(subjInd)]);
        % adjust the storage indices
        istart = istart + numPCA; iend = iend + numPCA;
        % Load the data, shape it to the right dimension
        image = load_nii(niifiles{subjInd}); res = reshape(image.img,[], k)';
        % Center the data
        X_tilde_all = res(:,validVoxels); [X_tilde_all, ] = remmean(X_tilde_all);
        % run pca on X_tilde_all
        [U_incr, D_incr] = pcamat(X_tilde_all);
        lambda = sort(diag(D_incr),'descend');
        sigma2_ML = sum(lambda(numPCA+1:length(lambda))) / (length(lambda)-numPCA);   
        U_q = U_incr(:,(size(U_incr,2)-numPCA+1):size(U_incr,2));
        D_q = diag(D_incr((size(U_incr,2)-numPCA+1):size(U_incr,2),...
            (size(U_incr,2)-numPCA+1):size(U_incr,2)));
        % Store the important pieces
        stage1White(:, :, subjInd) = diag((D_q-sigma2_ML).^(-1/2)) * U_q';
        stage1deWhite(:, :, subjInd) = U_q * diag((D_q-sigma2_ML) .^ (1/2));
        stackedPCAData(istart:iend, :) = stage1White(:, :, subjInd) * X_tilde_all;
        % Store the other sized whitening and dewhitening matrices for the 
        % reduction to q components. This will be used later to construct
        % Ytilde when getting starting values for the A matrix.
        U_q = U_incr(:,(size(U_incr,2)-q+1):size(U_incr,2));
        D_q = diag(D_incr((size(U_incr,2)-q+1):size(U_incr,2),...
            (size(U_incr,2)-q+1):size(U_incr,2)));
        % sum across all the remaining eigenvalues
        sigma2_ML = sum(lambda(q+1:length(lambda))) / (length(lambda)-q);
        % whitening, dewhitening matrix;
        my_whiteningMatrix = diag((D_q-sigma2_ML).^(-1/2)) * U_q';
        my_dewhiteningMatrix = U_q * diag((D_q-sigma2_ML) .^ (1/2));
        redICYWhite(:,:,subjInd) = my_whiteningMatrix;
        redICYdewhite(:,:,subjInd) = my_dewhiteningMatrix;
    end

    %% Stage 2 - Perform the second stage of reduction on the block of data

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
    
    stackedStage2Data = stackedPCAData;

end % end of loop over blocks

%% Stage 3 - Perform PCA down to q components
[U_incr, D_incr] = pcamat(stackedStage2Data);
lambda = sort(diag(D_incr),'descend');
sigma2_ML = sum(lambda(q+1:length(lambda))) / (length(lambda)-q);   
% Recover the whitening matrix (G inverse)
U_q = U_incr(:,(size(U_incr,2)-q+1):size(U_incr,2));
D_q = diag(D_incr((size(U_incr,2)-q+1):size(U_incr,2),...
        (size(U_incr,2)-q+1):size(U_incr,2)));
Ginv = diag((D_q).^(-1/2)) * U_q';
groupPCA_X = Ginv * stackedStage2Data; % this is X that the ICA is performed on

%% Stage 4 - ICA
[icasig, Amix, W] = fastica_full(groupPCA_X, 'numOFIC', q);
Agrp = Amix;

% view icasig


%% Reconstruct the single subject matrices - Over Blocks of Voxels

H = stage2deWhite;
G = pinv(Ginv);

%%%%%%%
for iVoxBlock = 1:nSubICBlocks 

    % Create storage for this round of subject level ICs
    disp('might should be n in group')
    singleSubjectICs = zeros(q, voxPerBlock(iVoxBlock), N);
    voxStart = 1;
    if iVoxBlock > 1
        voxStart = 1 + sum(voxPerBlock(1:(iVoxBlock-1)));
    end
    voxEnd = sum(voxPerBlock(1:iVoxBlock));
    disp(['Estimating for voxels ' num2str(voxStart) ' through ' num2str(voxEnd)])

    % Loop through and create the subject level ICs for the relevant voxels
    glength = size(G, 1) / nBlocks;
    for iGroup = 1:nBlocks
        istart = -numPCA + 1; iend = 0;
        disp('this is wrong index for stacked stage 2')
        gstart = glength*(iGroup-1) + 1; gend = glength*iGroup;
        for iSubj = 1:subjPerBlock(iGroup)
            istart = istart + numPCA; iend = iend +  numPCA;
            %HGAi = squeeze(H(istart:iend, :, iGroup)) * G(gstart:gend, :) * Agrp;
            %GAi = G(gstart:gend, :) * Agrp;
            GAi = G(istart:iend, :) * Agrp;
            %singleSubjectICs(:,:,iSubj) = pinv(GAi) *...
                %squeeze(stage2White(istart:iend, :, iGroup))'  *...
                %stackedStage2Data(istart:iend, voxStart:voxEnd);
            disp('i think this is also the wrong stage 2 index, pca2?')
            singleSubjectICs(:,:,iSubj) = pinv(GAi) * stackedStage2Data(istart:iend, voxStart:voxEnd);
        end
    end
    
    disp('remove scaling')
    singleSubjectICs = singleSubjectICs ./ max(abs(singleSubjectICs(:)));
    
    % Create a design matrix that includes an intercept (this is S0)
    Xfull = [ ones(N, 1), X];
    
    % Start getting estimates
    %epsilon2temp = zeros(q, V, N);
    voxRange = voxStart:voxEnd;
    for v = 1:voxPerBlock(iVoxBlock)
        estimate = (Xfull'*Xfull)^(-1) * Xfull' * squeeze(singleSubjectICs(:,v,:))';
        beta(:,:,voxRange(v)) = estimate(2:(p), :);
        S0(:,voxRange(v)) = estimate(1,:)';
        %epsilon2temp(:,voxRange(v),:) = singleSubjectICs(:,v,:) - reshape( (Xfull * estimate)', [q,1,N] );
    end
    %sigma2_sq = var(reshape(epsilon2temp, [q,V*N]), 0, 2);

end    

disp('Sigma2_sq not estimated')
%% Reconstruct for each N
singleSubjectIC = zeros(q, voxPerBlock(iVoxBlock));
% Storage for the aggregate map
aggregate_IC_map = zeros(q, voxPerBlock(iVoxBlock));
index = 0; % index controlling subjects
for iGroup = 1:nBlocks
    istart = -numPCA + 1; iend = 0;
    gstart = glength*(iGroup-1) + 1; gend = glength*iGroup;
    for iSubj = 1:subjPerBlock(iGroup)
        % Recover the subjects IC
        index = index + 1; % increment subject
        istart = istart + numPCA; iend = iend +  numPCA;
        GAi = G(istart:iend, :) * Agrp;
        disp('this indexing is for wrong level')
        singleSubjectIC = pinv(GAi) * stackedStage2Data(istart:iend, :);
        aggregate_IC_map = aggregate_IC_map + 1/N * singleSubjectIC ;
        % Generate the corresponding Ytilde variable, which is the single subject
        %     data reduced to q ICs.
        image = load_nii(niifiles{index}); res = reshape(image.img,[], k)';
        % center the data
        X_tilde_all = res(:,validVoxels); [X_tilde_all, ] = remmean(X_tilde_all);
        Ytilde = redICYWhite(:,:, index) * X_tilde_all;
        % Calculate the estimate for the mixing matrix
        cS_i = singleSubjectIC; %sInd = q*(iSubj-1)+1; eInd = iSubj*q;
        A_tempi = (cS_i * cS_i')^(-1) * cS_i * Ytilde';
        Asym = A_tempi';
        A(:,:,index) = Asym*real(inv(Asym'*Asym)^(1/2)); 
    end
end

disp('single subject error not calculated')

% Initial Guess: fit a Gaussian mixture
m=2;

for j =1:q
    GMModel = fitgmdist(S0(j,:)' ,m+1);
    id = find(abs(GMModel.mu) == max(abs(GMModel.mu)));
    theta.miu3(1+m*(j-1): m*j, 1) =[GMModel.mu(id), 0];
    idzero = abs(GMModel.mu) == min(abs(GMModel.mu));
    theta.sigma3_sq(1+m*(j-1): m*j, 1) = [GMModel.Sigma(id), GMModel.Sigma(idzero)];
    theta.pi(1+m*(j-1): m*j, 1)  =[GMModel.PComponents(id), 1-GMModel.PComponents(id)];
end

% create the final variables to return (beta already created)
theta.sigma1_sq = 9999;
theta.sigma2_sq = 9999;
theta.A = A;
grpSig = S0;

disp('Change to icasig')
%s0_agg = S0;
s0_agg = aggregate_IC_map;

% PCA whitening and dewhitening matrices for working backwards to ytilde
PCAquant.redICYWhite = redICYWhite;
PCAquant.redICYdewhite = redICYdewhite;
PCAquant.stage1White = stage1White;
PCAquant.stage2White = stage2White;
PCAquant.stage1dewhite = stage1deWhite;
PCAquant.stage2dewhite = stage2deWhite;
PCAquant.Agrp = Agrp;
PCAquant.numPCA2 = numPCA2;
PCAquant.G = G;
PCAquant.nBlocks = nBlocks;
PCAquant.subjPerBlock = subjPerBlock;
PCAquant.nSubICBlocks = nSubICBlocks;
PCAquant.voxPerBlock = voxPerBlock;

end


