%Reads the B field for the phantom from the Ansys HFSS data and stores the
%fields for each channel in *.mat files

clc;
clear;
close all;

for ch = 1:8

    path_directory = ''; % enter your path directory
    addpath(genpath(path_directory));
    
    
    % Load HFSS simulation data results
    
    delimiterIn = ' ';
    headerlinesIn = 2;
    for current_source = 1:14
        current_source
        if current_source <8
            filename = 'UM'+ string(ch) + string(current_source)+'.fld'; %change the name of your file accordingly
        else
            filename = 'LM'+ string(ch) + string(current_source-7)+'.fld'; %change the name of your file accordingly
        end
        Bfield{current_source} = importdata(filename,delimiterIn,headerlinesIn);
        Bfielddata{current_source}=Bfield{current_source}.data;
    
    end
    
    rmpath(path_directory);
    
    %%
    
    U_BS1 = Bfield{1};
    U_BS2 = Bfield{2};
    U_BS3 = Bfield{3};
    U_BS4 = Bfield{4};
    U_BS5 = Bfield{5};
    U_BS6 = Bfield{6};
    U_BS7 = Bfield{7};
    
    L_BS1 = Bfield{8};
    L_BS2 = Bfield{9};
    L_BS3 = Bfield{10};
    L_BS4 = Bfield{11};
    L_BS5 = Bfield{12};
    L_BS6 = Bfield{13};
    L_BS7 = Bfield{14};
    
    
    U_Bs1=U_BS1.data;
    U_Bs2=U_BS2.data;
    U_Bs3=U_BS3.data;
    U_Bs4=U_BS4.data;
    U_Bs5=U_BS5.data;
    U_Bs6=U_BS6.data;
    U_Bs7=U_BS7.data;
    
    L_Bs1=L_BS1.data;
    L_Bs2=L_BS2.data;
    L_Bs3=L_BS3.data;
    L_Bs4=L_BS4.data;
    L_Bs5=L_BS5.data;
    L_Bs6=L_BS6.data;
    L_Bs7=L_BS7.data;
    
    Bs(1,:,:)=U_Bs1(:,:);
    Bs(2,:,:)=U_Bs2(:,:);
    Bs(3,:,:)=U_Bs3(:,:);
    Bs(4,:,:)=U_Bs4(:,:);
    Bs(5,:,:)=U_Bs5(:,:);
    Bs(6,:,:)=U_Bs6(:,:);
    Bs(7,:,:)=U_Bs7(:,:);
    
    Bs(8,:,:)=L_Bs1(:,:);
    Bs(9,:,:)=L_Bs2(:,:);
    Bs(10,:,:)=L_Bs3(:,:);
    Bs(11,:,:)=L_Bs4(:,:);
    Bs(12,:,:)=L_Bs5(:,:);
    Bs(13,:,:)=L_Bs6(:,:);
    Bs(14,:,:)=L_Bs7(:,:);
    
     
    for i=1:14
    
        B_asymmetric_sphere(i,:)=[Bs(i,:,4)+1j.*(Bs(i,:,5)) Bs(i,:,6)+1j.*(Bs(i,:,7)) Bs(i,:,8)+1j.*(Bs(i,:,9))];
    
    end
    
    nameB = "B_asymmetric_sphere"+string(ch)+".mat";
    save(nameB,"B_asymmetric_sphere")

end