
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
addpath('C:\Users\petnov\Dropbox\shared - copy\MESHES\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');


enableVTK;
pats={'06'};
    patient_nr=pats{1};
    %mkdir 'D:\ARVC meshing automatic\patients\patient01\regions';
H1 = read_VTK(  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\HEART_0.70.vtk')); %our target heart
    [EPI1,LV1,RV1,~,MAT1] = HEARTparts( H1 );
    patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\');

    %H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart!
    H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste06_coarse','\HEART'));
    H_mesh.xyz = H_mesh.xyz*10;
    
  
    mesh_cell_property2 =  load(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'AHA_regions.mat'));
        
    mesh_cell_property2= mesh_cell_property2.mesh_cell_property2;
    
%     %scar type 5 - 1 conducting channel with
%     non uniform density borders and varying endo-epi surfaces
%      specify the params of the algorithm in terms of 
   % position of the scar, select  using cursor selection tool (select data tips (2nd row on the top right plot menu in Matlab
   % 2019, select point and right click on export cursor Data to workspace)
   % in matlab)
   % loop variable defines conduction slowing in surface of border of scar region 
     
 figure()
 headlight
 hold on
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','facecolor','b')
        
 distance= @(vec,vec1) ((vec(:,1)-vec1(:,1)).^2+(vec(:,2)-vec1(:,2)).^2+(vec(:,3)-vec1(:,3)).^2).^0.5;

 
        for j=[3,4,5] %scar size
            for i=9 %rest of scar
                for t=[1,1.5,2] %epi border scar
                    for r=[1,1.5,2] %endo border scar
               k=0;
             
                       %  scar_properties ={};

                x= cursor_info.Position(1);
                     y=    cursor_info.Position(2);
                     z= cursor_info.Position(3);

             ts_scaling_no_fib=ones(23,1);

             %make large LGE patch homogenious outside
                
                %get distance to the nearst epicardium and endocardial
                %surface points from the initial point selected
               epi1= knnsearch(H_mesh.xyz(  H_mesh.epi,:),[x,y,z]);
         endo1=knnsearch(H_mesh.xyz(H_mesh.rv,:),H_mesh.xyz(epi1,:));
         
         %find the mid point of the scar in the centre of the myocardium
         mid1= (H_mesh.xyz(H_mesh.epi(epi1),:)+H_mesh.xyz(H_mesh.rv(endo1),:))/2;
         %find the distance to the centre of the myocardium
         scar_dists=distance(H_mesh.xyz,mid1);
         
         %whole scar is between 3 and 5mm * scaling factor away from the
         %mid point (sperical shell)
             border_ids=find(scar_dists>3*j & scar_dists <5*j);
             %middle of it will be a scar border
             densest_ids = find(scar_dists<(4*j)+0.35 & scar_dists >4*j-0.35); %new addition, generate scar border from here
             
            %now find distance to nearest point in densest_ids
            %border ids is the whole fibrotic scar
            
            indexes_closest_mid=knnsearch(H_mesh.xyz(densest_ids,:),H_mesh.xyz(border_ids,:));
            
            %find distance to closest mid_point of scar
            dists_to_centre_scar=distance(H_mesh.xyz(border_ids,:),H_mesh.xyz(densest_ids(indexes_closest_mid),:));
            
            scalings_prev=load('D:/ARVC meshing automatic/patients/patient06/scalings/eikonal06_coarse_fib_0_scar_inhomog_grad_4_size_3_epi_intact_0.csv');
            
            %start by assigning an array of ones 
           no_lge= ts_scaling_no_fib(mesh_cell_property2);
           
           %normalise the distance 
           dists_to_centre_scar = dists_to_centre_scar/max(dists_to_centre_scar);
           %assign the scar scalings according to normalised distance
           
           
           dists_to_centre_scar_inv = abs(1-dists_to_centre_scar);
           
           %   no_lge(border_ids)=1./dists_to_centre_scar;
           
              y= @ (m,lin_dists)(m*lin_dists+1);
              
              
              scar_tissue=y(4,dists_to_centre_scar_inv);
              no_lge(border_ids)=scar_tissue;
%             for l=1:size(border_ids)
%                 dummy=find(H_mesh.epi == border_ids(l));
%                 if isempty(dummy)==0
%                     no_lge(border_ids(l))=1;
%                 end
%             end

          
            cum_endo=[];
            cum_epi=[];
            
            %this section should make the endocardium/epicardium change
            %conduction properties 
            
             %find distance to centre of scar from 
          
             for l=1:size(border_ids)
                dummy=find(H_mesh.epi == border_ids(l));
                if isempty(dummy)==0
                    no_lge(border_ids(l))=t*scar_tissue(l);
                    cum_epi=[cum_epi,t*scar_tissue(l)];
                end
            end
               for l=1:size(border_ids)
                dummy=find(H_mesh.rv == border_ids(l));
                if isempty(dummy)==0
                    no_lge(border_ids(l))=r*scar_tissue(l);
                    cum_endo=[cum_endo,r*scar_tissue(l)];
                    
                end
               end

  
              
              
%             dummy=find(cum_epi<1);
%             cum_epi(dummy)=1;
%             dummy=find(cum_endo<1);
%             cum_endo(dummy)=1;
%             
%             scar_properties{r,t,1}=cum_epi;
%             scar_properties{r,t,2}=cum_endo;
%             
           
           % dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(0),'.csv'), no_lge,'precision',10);
           dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(0),'_RV_anterior_scar_border_slowing_',num2str(i),'epi_scar_',num2str(t),'endo_scar',num2str(r),'_size_',num2str(j),'.csv'), no_lge,'precision',10);
         %   dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(2),'_scar.csv'),  int_lge,'precision',10);
          %  dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(3),'_scar.csv'),  adv_lge,'precision',10);
           end
            end
        
        end
        end
        
        

 figure()
 headlight
 hold on
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',no_lge,'facecolor','interp')
 colorbar
 axis off
 
 xs=linspace(0,1,50);
 figure()
 plot(xs,y(xs,2),'LineWidth',3)
 hold on
  plot(xs,y(xs,4),'LineWidth',3)
 plot(xs,y(xs,8),'LineWidth',3)
 legend({'typical ACM scar','pronounced ACM scar', 'severe ACM scar'},'Fontsize',35)
 xlabel('normalised distance from middle of scar to border','Fontsize',35)
 ylabel('conduction scaling factor','Fontsize',35)
 set(gca,'Fontsize',25)
 
 
 scar_old_long=load('D:/ARVC meshing automatic/patients/patient06/scalings/eikonal06_coarse_fib_0_scar_along_fiber_slowing_4_length_10_radius_9.csv');
 
 positions=load('D:/ARVC meshing automatic/patients/patient06/mpp/ECG_ELECTRODES.mat');
 
 pos2=load('D:/ARVC meshing automatic/patients/patient06/DTI003_electrodePositions.csv');
 
 
 
 figure()
 headlight
 hold on
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',scar_old_long,'facecolor','interp')
 colorbar

       