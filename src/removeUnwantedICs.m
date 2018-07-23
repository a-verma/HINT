function [ aaa, bbb ] = removeUnwantedICs( keeplist, PCAquant, niifiles)
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

N = length(niifiles);
qstar = sum(keeplist);

%% Recover the data sizes used to obtain the initial guess;
nBlocks = PCAquant.nBlocks;
subjPerBlock = PCAquant.subjPerBlock;
nSubICBlocks = PCAquant.nSubICBlocks;
voxPerBlock = PCAquant.voxPerBlock;

%% Figure out sizes based on the new, reduced Q

%% Step 1 - Back reconstruct a "removal dataset" of the original subject data from the dewhitening
% matrices using the columns of the mixing matrix corresponding to ICs to
% be removed. Then load the original subject data and subtract off this
% removal dataset.
disp('Removing unwanted components from the data...')

% Create a mixing matrix with only the UNWANTED columns
Aremove = PCAquant.Agrp;
Aremove(:, find(keeplist == 1)) = 0;

% Back reconstruct a "unwanted" dataset, Yremove
for iSubj = 1:N
    % Load the subjects data
    image = load_nii(niifiles{iSubj}); res = reshape(image.img,[], k)';
    % Back reconstruct a "removal" dataset - first need subject specific
    % dewhitening matrix
    creationMatrix = 
    % Create the removal dataset
    removeImage = Aremove 
end


%% old julia code

% 
%     # Preallocate space
%     stackedPCA = zeros(nPC*N, V)
% 
%     # Calculate the new stackedPCAData
%     G = pinv(Ginv)
%     istart = -nPC+1; iend = 0
%     for iSubj = 1:N
%         istart += nPC
%         iend += nPC
%         GAi = G[istart:iend, :] * Anew
%         stackedPCA[istart:iend, :] = GAi * singleSubjectICs[:,:,iSubj]
%     end
% 
%     # Define some things
%     H_matrix = 0.0
%     H_matrix_inv = 0.0
%     Y_tilde_all = 0.0
% 
%     # Get the cleaned version of the original data and reduce it to qstar ICs
%     istart = -nPC+1; iend = 0
%     print("here")
%     for iSubj = 1:N
%         istart += nPC
%         iend += nPC
%         # TODO should it be transpose then inverse?
%         yraw = pinv(pcaWhiteMat[:,:,iSubj])' * stackedPCA[istart:iend, :]
% 
%         # Calculate X tilde all
%         X_tilde_all = yraw .- mean(yraw, 2)
% 
%         my_whiteningMatrix, my_dewhiteningMatrix, my_whitesig = eig_reduce(X_tilde_all, qstar)
% 
%         if (iSubj == 1)
%             # transform matrix for the two-stage dim reduction and whitening;
%             H_matrix = my_whiteningMatrix * (eye(T) - 1/T * ones(T, T))
%             H_matrix_inv = my_dewhiteningMatrix
%             Y_tilde_all = my_whitesig
%         else
%             # The below line is new from matlab version- Josh
%             H_matrix_new = my_whiteningMatrix * (eye(T) - 1/T * ones(T, T))
%             H_matrix = cat([1,2], H_matrix, H_matrix_new)
%             #H_matrix = blkdiag(H_matrix, H_matrix_new);
%             H_matrix_inv = cat([1,2], H_matrix_inv, my_dewhiteningMatrix)
%             #H_matrix_inv = blkdiag(H_matrix_inv, my_dewhiteningMatrix);
%             Y_tilde_all = vcat(Y_tilde_all, my_whitesig)
%         end
% 
%     end
% 
%     # Now calculate the new ytilde based on the cleaned original data
%     YtildeStar = Y_tilde_all
%     C_matrix = H_matrix * H_matrix'
%     C_matrix_diag = diag(C_matrix)
% 
%     return(stackedPCA, YtildeStar, qstar, C_matrix_diag)

end

