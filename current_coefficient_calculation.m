clc;
clear;
close all;

% This code finds the current source coefficients for the model
% For this code to work, the NOVA data and the Ansys HFSS data must be aligned

%%
% For NOVA results in matrix BE
% Step 1: identify and choose the regions containing the phantom
% Step 2: reshape them into 1D 
% Step 3: only take the points that contain the phantom (using the mask)
% Step 4: x, y, z point vectors are concatenated to form a single vector

% Load NOVA data
load('phantom1S_coils1-8.mat')
B_NOVA = B(:,:,31:113,31:113,19:101);

% Repeat for each channel (can be turned into a loop)
ch=8; % Rerun the code for channels from 1 to 8, change in each iteration

% Initialize the matrices (channels, x_pixels, y_pixels, z_pixels)
B1_est = zeros(8,83,83,83); % For simulation data
Nova_B1s = zeros(8,83,83,83); % For NOVA data

% Separate NOVA data into x, y, z components
BBX(:,:,:)=B(ch,1,31:113,31:113,19:101); % We only take the 3D rectangular volume with the phantom
BBY(:,:,:)=B(ch,2,31:113,31:113,19:101);
BBZ(:,:,:)=B(ch,3,31:113,31:113,19:101);

% For the phantom mask we use the conductivity map
maskbb(:,:,:)=matcond(1,31:113,31:113,19:101); % We only take the 3D rectangular volume with the phantom


% Reshape all of the matrices to 1D vectors
BX=reshape(BBX,[],1);
BY=reshape(BBY,[],1);
BZ=reshape(BBZ,[],1);

maskb=reshape(maskbb,[],1);

% Keep only the data points inside the mask (the data points containing the phantom)
q = find(maskb);
m=0;
for i=1:size(q,1)
    m=m+1;
    BXNNZ(m,1)=BX(q(m));
    BYNNZ(m,1)=BY(q(m));
    BZNNZ(m,1)=BZ(q(m));
end

BE=[BXNNZ;BYNNZ;BZNNZ]; % Combine x,y,z data points into a column vector

%% 
% For Ansys HFSS simulation results in matrix BES
% This loop extracts the data for Ansys simulation results:
% Step 1: turn them into 3D
% Step 2: reshape them into 1D 
% Step 3: only take the points that contain the phantom (using the mask)
% Step 4: x, y, z point vectors are concatenated to form a single vector
% Step 5: combine 14 current sources into 3 groups

BES = [];
for ch_temp =1:8
    % Load Ansys HFSS simulation results for ch_temp
    load("B_asymmetric_sphere"+string(ch_temp)+".mat", "B_asymmetric_sphere");

    % Step 1: reshape Ansys matrices to 3D 
    for n=1:14  % Repeat for each current source in the channel
          m=1;
        for i=1:83 % For each pixel in x
            for j=1:83 % For each pixel in y
                for k=1:83 % For each pixel in z
                    BS3DX(i,j,k)=B_asymmetric_sphere(n,m);
                    BS3DY(i,j,k)=B_asymmetric_sphere(n,m+571787);
                    BS3DZ(i,j,k)=B_asymmetric_sphere(n,m+2*571787);
                    m=m+1;
                        
                end
            end
        end
        % Step 2: reshape matrices to 1D 
        BsX=reshape(BS3DX.*maskbb,[],1);
        BsY=reshape(BS3DY.*maskbb,[],1);
        BsZ=reshape(BS3DZ.*maskbb,[],1);

        % Step 3: only take the points inside the phantom 
        l=0;
        for i=1:size(q,1)
            l=l+1;
            BsXNNZ(l,1)=BsX(q(l));
            BsYNNZ(l,1)=BsY(q(l));
            BsZNNZ(l,1)=BsZ(q(l));
        end
        BEsV(:,n)=[BsXNNZ;BsYNNZ;BsZNNZ]; % Step 4: concatenate x,y,z vectors

    end

    %Step 5: combine current sources to reduce degrees of freedom
    upper1 = BEsV(:,1)+BEsV(:,2)-BEsV(:,3)-BEsV(:,4);
    upper2 = BEsV(:,5)+BEsV(:,7);
    upper3 = BEsV(:,6);
    lower1 = BEsV(:,8)+BEsV(:,9)-BEsV(:,10)-BEsV(:,11);
    lower2 = BEsV(:,12)+BEsV(:,14);
    lower3 = BEsV(:,13);
     
    % temp_BES is a temporary matrix containing magnetic fields for one channel
    temp_BES= [upper1, lower1, upper2+upper3+lower2+lower3];

    BES = [BES,temp_BES]; % matrix containing magnetic fields for all channels
end

%% 

% This part of the code calculates the coefficients for the current sources (alphas) using the matrices BE (NOVA magnetic or electric field) and BES (simulation magnetic or electric field)

% Step 1: Normalize the matrices
% Step 2: Singular value decomposition for the Ansys HFSS simulation data (matrix BES)
% Step 3: Calculate cumulative energy
% Step 4: Keep the parts of the matrix that retain 99.9% of the cumulative energy
% Step 5: Calculate the current source coefficients from the remaining matrix with SVD based regularization
% Step 6: Check the effective condition number
% Step 7: Rescale the matrix

% Step 1: normalize the matrices

scale_A = norm(BES, 'fro');
scale_A_objective = norm(BE);

BES_scaled = BES / scale_A;
BE_scaled  = BE  / scale_A_objective;

% Step 2: SVD for BES

[U, S, V] = svd(BES_scaled, "econ");

singular_values = diag(S); % Singular values

% Step 3: Calculate cumulative energy

cumulative_energy = cumsum(singular_values.^2) / sum(singular_values.^2);

% Step 4: keep the parts of the matrix that retain 99.99% of the cumulative energy

retain_percentage = 0.9999; % you can change this according to how much you want to keep

k = find(cumulative_energy >= retain_percentage, 1); %Remove singular values after index k
idx_keep = 1:k;
U_k = U(:, idx_keep);
S_k = S(idx_keep, idx_keep);
V_k = V(:, idx_keep);


column_contributions = sum(abs(V_k).^2, 2); % This is the sum of squared contributions across the removed right singular vectors
[contr, sorted_column_indices] = sort(column_contributions, 'descend');

% Display
disp('Columns contributing to removed components:');
disp([contr, sorted_column_indices]);


% Step 5: Calculate the current source coefficients from the remaining matrix with SVD based regularization

lambda = 10^(-10);
s = diag(S_k);
alpha_scaled = V_k * diag(s ./ (s.^2 + lambda)) * (U_k' * BE_scaled);

% Step 6: Check the effective condition number

filter_factors = s ./ (s.^2 + lambda);  % Regularization filter
cond_eff = max(filter_factors) / min(filter_factors);
disp(['Effective condition number (SVD-based): ', num2str(cond_eff)]);

% Step 7: Rescale the matrix of alphas back to original units
alpha = alpha_scaled * (scale_A_objective/scale_A);

%%

%Step 1: Take the magnitudes and angles of the current coefficients (alphas)
%Step 2: Save the current source coefficients (alphas) in matrix CC for each channel
%Step 3: Save CC for future use in a *.mat file


CC = [];

for temp_ch= 1:8

    a1=alpha(1+(temp_ch-1)*3);
    a2=alpha(2+(temp_ch-1)*3);
    a3=alpha(3+(temp_ch-1)*3);

    %Step 1: Take the magnitudes and angles of the current coefficients (alphas)

    a1abs= abs(a1);
    a1ang= rad2deg(angle(a1));
    a2abs= abs(a2);
    a2ang= rad2deg(angle(a2));
    a3abs= abs(a3);
    a3ang= rad2deg(angle(a3));


    temp_CC=[a1abs a2abs a3abs; a1ang a2ang a3ang];

    %Step 2: Save the current source coefficients (alphas) in matrix CC for each channel

    CC = [CC;temp_CC];
end

%Step 3: Save CC for future use in a *.mat file

save("current_coefficients_for_channel_"+string(ch)+".mat", "CC")
disp("Channel "+string(ch)+" Finished :)")





