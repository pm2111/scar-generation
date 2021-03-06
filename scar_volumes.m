
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
addpath('C:\Users\petnov\Dropbox\shared - copy\MESHES\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\'); 

enableVTK;

    %mkdir 'D:\ARVC meshing automatic\patients\patient01\regions';
    patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\DTI003\');

    %H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart
H=struct;
H.xyz=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_xyz.csv');
H.tri=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_tri.csv');
H.rv=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_rvface.csv');
H.lv=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_lvface.csv');
H.epi=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_epiface.csv');
H.tricentres=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_tetrahedronCenters.csv');

H.face=vertcat(H.epi,H.lv,H.rv);
  H.xyz=H.xyz;

names =dir(strcat(patient_sim_folder,'scars_LV_anterior\'));
scar = dlmread(strcat(patient_sim_folder,'scars_LV_anterior\',names(3).name));               

[EPI,LV,RV,LID,MAT] = HEARTparts( H ); %separate into epicardial, left, right,
%lid (cover of ventricles) surfaces and the transformation matrix

figure()
plotMESH( transform( EPI , MAT ) ,'ne'); headlight %use the transform tool to bring the mesh to the normalised orientation (facing upwards rel to the xy plane)
hplotMESH( transform( LV  , MAT ) ,'r','ne'); 
hplotMESH( transform( RV  , MAT ) ,'g','ne'); 
hplotMESH( transform( LID , MAT ) ,'c','ne'); 


plotMESH( MeshFillHoles( LV ) ,'ne' ); headlight %in order to use Gauss' formula to compute volumes, we need to fill any holes left


muscle_LV=meshVolume( MeshFillHoles( LV ) ) %use gauss' thm to compute the mesh volume (NB all normals must point in the same direction)
muscle_RV=meshVolume( MeshFillHoles( RV ) ) %use gauss' thm to compute the mesh volume (NB all normals must point in the same direction)

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

figure()
hplotMESH(MeshFillHoles( H_scar),'ne')
headlight

meshVolume(MeshFillHoles(H_scar));


