function [scarVolumeFraction]=scar_volumes_H_mesh_fcn(H)
%location is the name of the folder containing the scars for a specific
%location
%ventricle is a double, either 1 if the scar is in the LV or other if the
%RV. The septum is assumed to be part of the LV.
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts')) %include the distance function
addpath('C:\Users\petnov\Dropbox\shared - copy\MESHES\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\'); 

enableVTK;

    %mkdir 'D:\ARVC meshing automatic\patients\patient01\regions';
    patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\DTI003\');

    %H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart

    H=H;
  
    if max(H.lge)==1
        scarVolumeFraction =0;
        return 
    end
figure()

scar = H.lge;        
   
%lid (cover of ventricles) surfaces and the transformation matrix
% 
% figure()
% plotMESH( transform( EPI , MAT ) ,'ne'); headlight %use the transform tool to bring the mesh to the normalised orientation (facing upwards rel to the xy plane)
% hplotMESH( transform( LV  , MAT ) ,'r','ne'); 
% hplotMESH( transform( RV  , MAT ) ,'g','ne'); 
% hplotMESH( transform( LID , MAT ) ,'c','ne'); 
% 
% 
% plotMESH( MeshFillHoles( LV ) ,'ne' ); headlight %in order to use Gauss' formula to compute volumes, we need to fill any holes left


total_muscle_mass=sum( meshQuality( MeshFillHoles( H ) , 'volume' ) ); %compute the mesh volume by summing the volumes of each tetrahedron in the mesh
%extract mesh of scar
%find tetra which are made up entirely of points from the scar region
ids_scar=find(scar>1);

closest_point_to_tetra_centre = knnsearch(H.xyz,H.tricentres);
tetra_scar_ids=[];
for i= 1:size(ids_scar)
    idx_scar=ids_scar(i);
    dummy=find(closest_point_to_tetra_centre==idx_scar);
    if isempty(dummy)==0
        tetra_scar_ids=vertcat(tetra_scar_ids,dummy); %id of the tetra centres
    end
end

H_scar =struct;
H_scar.xyz=H.xyz;
H_scar.tri = H.tri(tetra_scar_ids,:);

hplotMESH(MeshFillHoles( H_scar),'ne')
headlight
scarVolumeFraction= sum( meshQuality( MeshFillHoles( H_scar ) , 'volume' ) )*100/total_muscle_mass;

end
