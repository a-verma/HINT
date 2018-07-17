function [ theta, beta, grpSig, s0_agg ] = tcgicaHINT_ext(niifiles, maskf, validVoxels,...
    Ytilde, X, numPCA, q, N, time_num, V, prefix)

%% Stages:
% 1. PCA dimension reduction of each subject down to nPC x V
% 2. PCA dimension reduction of each block of nPC x V datasets down to nPC2
%             block = 1, ..., nBlock
% 3. PCA dimension reduction of stacked block data from
%        nPC2*nBlock x V   to   q x V
% 4. ICA on the q x V dataset to get initial values

T = time_num;
p = size(X, 2) + 1;
maxBytes = 1e9; disp('change max bytes to 8e9!!!')

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
numPCA2 = floor(maxBytesPerBlock ./ (8 * V));

% Determine how large the stacked subject ica data would be (split on
% voxels)
bytesRequiredSubjICs = 8 * q * V * N;
disp('change this to 8e9')
nSubICBlocks = ceil(bytesRequiredSubjICs / 0.2e9);
%nSubICBlocks = ceil(bytesRequiredSubjICs / 8e9);
voxPerBlock = floor(V / nSubICBlocks) * ones(nSubICBlocks, 1);
voxPerBlock(nSubICBlocks) = V - (nSubICBlocks-1)*floor(V / nSubICBlocks);

% Setup whitening and dewhiten storage for the intermediate steps
redICYWhite = zeros(q, T, N); redICYdewhite = zeros(T, q, N);
stage1White = zeros(numPCA, T, N); stage1dewhite = zeros(T, numPCA, N);
stage2White = zeros(numPCA2, numPCA*max(subjPerBlock), nBlocks); stage2dewhite = zeros(numPCA*max(subjPerBlock), numPCA2, nBlocks);
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
        stage1dewhite(:, :, subjInd) = U_q * diag((D_q-sigma2_ML) .^ (1/2));
        stackedPCAData(istart:iend, :) = stage1White(:, :, subjInd) * X_tilde_all;
    end

    %% Stage 2 - Perform the second stage of reduction on the block of data

    [U_incr_2, D_incr_2] = pcamat(stackedPCAData);
    lambda = sort(diag(D_incr),'descend');
    sigma2_ML = sum(lambda(numPCA2+1:length(lambda))) / (length(lambda)-numPCA2);   
    U_q = U_incr_2(:,(size(U_incr_2,2)-numPCA2+1):size(U_incr_2,2));
    D_q = diag(D_incr_2((size(U_incr_2,2)-numPCA2+1):size(U_incr_2,2),...
        (size(U_incr_2,2)-numPCA2+1):size(U_incr_2,2)));

    ng = numPCA * subjPerBlock(iGroup); % Account for possible different group size
    stage2White(:,1:ng,iGroup) = diag((D_q-sigma2_ML).^(-1/2)) * U_q';
    stage2dewhite(1:ng,:,iGroup) = U_q * diag((D_q-sigma2_ML) .^ (1/2));

    % Store the 2nd stage results
    storeStart = (iGroup-1)*numPCA2 + 1;
    storeEnd = numPCA2 * iGroup;
    stackedStage2Data(storeStart:storeEnd, :) = stage2White(:,1:ng,iGroup) * stackedPCAData;

end % end of loop over blocks

%% Stage 3 - Perform PCA down to q components
[U_incr, D_incr] = pcamat(stackedStage2Data);
% Recover the whitening matrix (G inverse)
U_q = U_incr(:,(size(U_incr,2)-q+1):size(U_incr,2));
D_q = diag(D_incr((size(U_incr,2)-q+1):size(U_incr,2),...
        (size(U_incr,2)-q+1):size(U_incr,2)));
Ginv = diag((D_q).^(-1/2)) * U_q';
groupPCA_X = Ginv * stackedStage2Data; % this is X that the ICA is performed on

%% Stage 4 - ICA
[icasig, Amix, W] = fastica_full(groupPCA_X, 'numOFIC', q);
Agrp = Amix';

%% Reconstruct the single subject matrices - Over Blocks of Voxels

H = stage2dewhite;
G = pinv(Ginv);

%%%%%%%
for iVoxBlock = 1:nSubICBlocks 

    % Create storage for this round of subject level ICs
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
        gstart = glength*(iGroup-1) + 1; gend = glength*iGroup;
        for iSubj = 1:subjPerBlock(iGroup)
            istart = istart + numPCA; iend = iend +  numPCA;
            HGAi = squeeze(H(istart:iend, :, iGroup)) * G(gstart:gend, :) * Agrp;
            singleSubjectICs(:,:,iSubj) = pinv(HGAi) * stackedPCAData(istart:iend, voxStart:voxEnd);
        end
    end
    
    %disp('REMOVE THIS SCALING STEP')
    %singleSubjectICs = singleSubjectICs ./ max(abs(singleSubjectICs(:)));

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
index = 0; % index controlling subjects
for iGroup = 1:nBlocks
    istart = -numPCA + 1; iend = 0;
    gstart = glength*(iGroup-1) + 1; gend = glength*iGroup;
    for iSubj = 1:subjPerBlock(iGroup)
        % Recover the subjects IC
        index = index + 1;
        istart = istart + numPCA; iend = iend +  numPCA;
        HGAi = squeeze(H(istart:iend, :, iGroup)) * G(gstart:gend, :) * Agrp;
        singleSubjectIC = pinv(HGAi) * stackedPCAData(istart:iend, :);
        
        % Generate the corresponding Ytilde variable, which is the single subject
        %     data reduced to q ICs.
        image = load_nii(niifiles{iSubj}); res = reshape(image.img,[], k)';
        % center the data
        X_tilde_all = res(:,validVoxels); [X_tilde_all, ] = remmean(X_tilde_all);
        [U_incr, D_incr] = pcamat(X_tilde_all);
        % sort the eig values, IX:index
        lambda = sort(diag(D_incr),'descend');
        U_q = U_incr(:,(size(U_incr,2)-q+1):size(U_incr,2));
        D_q = diag(D_incr((size(U_incr,2)-q+1):size(U_incr,2),...
            (size(U_incr,2)-q+1):size(U_incr,2)));
        % sum across all the remaining eigenvalues
        sigma2_ML = sum(lambda(q+1:length(lambda))) / (length(lambda)-q);
        % whitening, dewhitening matrix and whitened data;
        my_whiteningMatrix = diag((D_q-sigma2_ML).^(-1/2)) * U_q';
        my_dewhiteningMatrix = U_q * diag((D_q-sigma2_ML) .^ (1/2));
        Ytilde = my_whiteningMatrix * X_tilde_all;
        
        % Calculate the estimate for the mixing matrix
        cS_i = singleSubjectIC; sInd = q*(iSubj-1)+1; eInd = iSubj*q;
        A_tempi = (cS_i * cS_i')^(-1) * cS_i * Ytilde';
        Asym = A_tempi';
        A(:,:,index) = Asym*real(inv(Asym'*Asym)^(1/2));
        
        % Calculate sigma1_squared (subject level error)
        
    end
end



% end of loop

% Calculate sigma1_squared (subject level error)
% errors = zeros(q, V, N);
% for i=1:N
%     sInd = q*(i-1)+1; eInd = i*q;
%     errors(:,:,i) = Ytilde(sInd:eInd,:) - A(:,:,i) * singleSubjectICs(:,:,i);
% end
% sigma1_sq = var( reshape(errors, [1, q * V * N] ) );
disp('single subject error not calculated')


%% \:
% 1: Creates Ytilde, a q x V data set of reduced dimension - TODO remove
% maybe?
% 2: Get all the whitening and dewhitening matrices for stage 1?

% stackedPCAData = zeros(numPCA*N, V);
% istart = -numPCA+1;
% iend = 0;
% k=time_num;
% deWhite = zeros(T, numPCA, N);
% for i=1:N
%     disp(['Stage 1 PCA for subject ' num2str(i)])
%     % adjust the storage indices
%     istart = istart + numPCA;
%     iend = iend + numPCA;
%     % Load the data
%     image = load_nii(niifiles{i});
%     res = reshape(image.img,[], k)';
%     % X tilde all is raw T x V subject level data for subject i
%     X_tilde_all = res(:,validVoxels);
%     % Center the data
%     [X_tilde_all, ] = remmean(X_tilde_all);
%     % run pca on X_tilde_all
%     [U_incr, D_incr] = pcamat(X_tilde_all);
%     % sort the eig values, IX:index
%     lambda = sort(diag(D_incr),'descend');
%     % Reduce to numPCA
%     U_q = U_incr(:,(size(U_incr,2)-numPCA+1):size(U_incr,2));
%     D_q = diag(D_incr((size(U_incr,2)-numPCA+1):size(U_incr,2),...
%         (size(U_incr,2)-numPCA+1):size(U_incr,2)));
%     % sum across all the remaining eigenvalues
%     sigma2_ML = sum(lambda(numPCA+1:length(lambda))) / (length(lambda)-numPCA);    
%     % whitening, dewhitening matrix and whitened data;
%     my_whiteningMatrix = diag((D_q-sigma2_ML).^(-1/2)) * U_q';
%     my_dewhiteningMatrix = U_q * diag((D_q-sigma2_ML) .^ (1/2));
%     deWhite(:,:,i) = my_dewhiteningMatrix;
%     my_whitesig = my_whiteningMatrix * X_tilde_all;
%     if (i == 1)
%         % transform matrix for the two-stage dim reduction and whitening;
%         H_matrix = my_whiteningMatrix * (eye(T)-1/T * ones(T));
%         H_matrix_inv = my_dewhiteningMatrix;
%         stackedPCAData(istart:iend, :) = my_whitesig;
%     else
%         newHmat = my_whiteningMatrix * (eye(T)-1/T * ones(T));
%         H_matrix = blkdiag(H_matrix, newHmat);
%         H_matrix_inv = blkdiag(H_matrix_inv, my_dewhiteningMatrix);
%         stackedPCAData(istart:iend, :) = my_whitesig;
%     end
%     
% end % end of PCA loop
% 
% disp('Starting new part')
% disp('Performing intermediate PCA reduction, this will need to be absorbed into earlier step')
% intermediatePCAdata1 = zeros(numPCA*N/2, V);
% intermediatePCAdata2 = zeros(numPCA*N/2, V);

% fill in the intermediate data holders, this will happen in a loop to save
% space
% intermediatePCAdata1 = stackedPCAData(1:(numPCA*N/2), :);
% intermediatePCAdata2 = stackedPCAData((numPCA*N/2+1):(numPCA*N), :);
% % perform a pca on each data group (should happen in a loop)
% npca2 = 45;
% [U_incr_2, D_incr_2] = pcamat(intermediatePCAdata1);
% U_q2 = U_incr_2(:,(size(U_incr_2,2)-npca2+1):size(U_incr_2,2));
% D_q2 = diag(D_incr_2((size(U_incr_2,2)-npca2+1):size(U_incr_2,2),...
%     (size(U_incr_2,2)-npca2+1):size(U_incr_2,2)));
% my_whiteningMatrix2 = diag((D_q2).^(-1/2)) * U_q2';
% my_dewhiteningMatrix2 = U_q2 * diag((D_q2) .^ (1/2));
% my_whitesig2 = my_whiteningMatrix2 * intermediatePCAdata1;
% [U_incr_2, D_incr_2] = pcamat(intermediatePCAdata2);
% U_q2 = U_incr_2(:,(size(U_incr_2,2)-npca2+1):size(U_incr_2,2));
% D_q2 = diag(D_incr_2((size(U_incr_2,2)-npca2+1):size(U_incr_2,2),...
%     (size(U_incr_2,2)-npca2+1):size(U_incr_2,2)));
% my_whiteningMatrix3 = diag((D_q2).^(-1/2)) * U_q2';
% my_dewhiteningMatrix3 = U_q2 * diag((D_q2) .^ (1/2));
% my_whitesig3 = my_whiteningMatrix3 * intermediatePCAdata2;
% % Stack the intermediate PCA level data
% stackedIntermediatePCA = [my_whitesig2; my_whitesig3];
% disp('Ending new part')
% 
% disp("Performing group level PCA")

% % This is the "GIFT" type step - Calculating Ginv for X = Ginv Finv Y
% [U_incr, D_incr] = pcamat(stackedIntermediatePCA);
% % Recover the whitening matrix (G inverse)
% U_q = U_incr(:,(size(U_incr,2)-q+1):size(U_incr,2));
% D_q = diag(D_incr((size(U_incr,2)-q+1):size(U_incr,2),...
%         (size(U_incr,2)-q+1):size(U_incr,2)));
% Ginv = diag((D_q).^(-1/2)) * U_q';
% groupPCA_X = Ginv * stackedIntermediatePCA; % this is X that the ICA is performed on

% % [icasig, A, W] = FASTICA (mixedsig); outputs the estimated separating
% % matrix W and the corresponding mixing matrix A.
% % CURRENTLY HWERE, WHY ICA RIGHT AFTER PCA WITH SAME DIM?
% [icasig, Amix, W] = fastica_full(groupPCA_X, 'numOFIC', q);
% Agrp = Amix';

% Recover the single subject matrices - think dimensions might have gone
% wrong somewhere
% singleSubjectICs = zeros(q, V, N);
% H = [my_dewhiteningMatrix2; my_dewhiteningMatrix3];
% istart = -numPCA+1;
% iend = 0;
% G = pinv(Ginv);
% for iSubj = 1:N
%     istart = istart + numPCA;
%     iend = iend +  numPCA;
%     gstart = 1;
%     gend = 45;
%     if iSubj > N/2
%         gstart = 46;
%         gend = 90;
%     end
%     HGAi = H(istart:iend, :) * G(gstart:gend, :) * Agrp;
%     singleSubjectICs(:,:,iSubj) = pinv(HGAi) * stackedPCAData(istart:iend, :);
% end
% 
% % Scale everything TODO figure out how this aligns with the qxV preproc data
% singleSubjectICs = singleSubjectICs ./ max(abs(singleSubjectICs(:)));
% 
% % Get the number of covariates
% p = size(X, 2) + 1;
% 
% % Create a design matrix that includes an intercept (this is S0)
% Xfull = [ ones(N, 1), X];
% 
% % Start getting estimates
% beta = zeros(p-1, q, V);
% S0 = zeros(q, V);
% epsilon2temp = zeros(q, V, N);
% for v = 1:V
%     estimate = (Xfull'*Xfull)^(-1) * Xfull' * squeeze(singleSubjectICs(:,v,:))';
%     beta(:,:,v) = estimate(2:(p), :);
%     S0(:,v) = estimate(1,:)';
%     epsilon2temp(:,v,:) = singleSubjectICs(:,v,:) - reshape( (Xfull * estimate)', [q,1,N] );
% end
% sigma2_sq = var(reshape(epsilon2temp, [q,V*N]), 0, 2);

% % Look at the estimates, how do they look?
% newImage = zeros(91, 109, 91);
%  newImage(validVoxels) = S0(1, :);
%  newImage(find(newImage == 0)) = nan;
% for ic = 1:q
% newImage(validVoxels) = S0(ic, :);
% for i = 1:91
%         disp([ 'IC ' num2str(ic) ' slice ' num2str(i) ])
%         imagesc(squeeze(newImage(:, :, i)))
%         pause(0.1)
% end
% end

% Save to nii files to see if things are working correctly
pth = '/Users/joshlukemire/Desktop/testICAINIGUESS/xstage_iniguess_';
anat = load_nii(maskf);
for ic = 1:q
    fname = [pth num2str(ic) '.nii'];
    newImage = zeros(91, 109, 91);
    newImage(validVoxels) = S0(ic, :);
    newIC = make_nii(newImage);
    newIC.hdr.hist.originator = anat.hdr.hist.originator;
    save_nii(newIC, fname);
end

% Get the subject level mixing matrices
% A = zeros(q, q, N);
% for i = 1:N
%     cS_i = singleSubjectICs(:,:,i);
%     sInd = q*(i-1)+1; eInd = i*q;
%     A_tempi = (cS_i * cS_i')^(-1) * cS_i * Ytilde[sInd:eInd,:]';
%     Asym = A_tempi';
%     A(:,:,i) = Asym * real(inv(Asym'*Asym)^(1/2));
% end
% 
% A = zeros(q,q,N);
% for i = 1:N
%     cS_i = singleSubjectICs(:,:,i); sInd = q*(i-1)+1; eInd = i*q;
%     A_tempi = (cS_i * cS_i')^(-1) * cS_i * Ytilde(sInd:eInd,:)';
%     Asym = A_tempi';
%     A(:,:,i) = Asym*real(inv(Asym'*Asym)^(1/2));
% end
% 
% % Calculate sigma1_squared (subject level error)
% errors = zeros(q, V, N);
% for i=1:N
%     sInd = q*(i-1)+1; eInd = i*q;
%     errors(:,:,i) = Ytilde(sInd:eInd,:) - A(:,:,i) * singleSubjectICs(:,:,i);
% end
% sigma1_sq = var( reshape(errors, [1, q * V * N] ) );
% 
% % Calculate sigma1_squared (subject level error)
% errors = zeros(q, V, N);
% for i=1:N
%     sInd = q*(i-1)+1; eInd = i*q;
%     errors(:,:,i) = Ytilde(sInd:eInd,:) - A(:,:,i)*singleSubjectICs(:,:,i);
% end
% sigma1_sq = var(reshape(errors, [1, q*V*N]));

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

s0_agg = icasig;

 

end


