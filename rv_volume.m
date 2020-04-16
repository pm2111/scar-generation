%function to compute the RV volume
function [rv_vol] = rv_volume(H_mesh)

expected_RV_wall_thickness = 0.4;


[nearest_endo_rv_point,dist_to_nearest_endo_point] = knnsearch(H_mesh.xyz(H_mesh.rv,:),H_mesh.xyz(H_mesh.epi,:));


rv_epicardial_node_indexes = H_mesh.epi(dist_to_nearest_endo_point<expected_RV_wall_thickness);



distance_epi_all_endo_rv = euclidian_distance(H_mesh.xyz(H_mesh.epi,:),H_mesh.xyz(nearest_endo_rv_point,:));






rv_endocardial_node_pair_indexes = knnsearch(H_mesh.xyz(H_mesh.rv,:),H_mesh.xyz(rv_epicardial_node_indexes,:));
rv_endocardial_node_pair_indexes =H_mesh.rv(rv_endocardial_node_pair_indexes); %keep this in whole mesh indexing


OA = H_mesh.xyz(rv_epicardial_node_indexes,:);
OB = H_mesh.xyz(rv_endocardial_node_pair_indexes,:);
AB = -OA + OB;
AC = AB/2;
OC = OA + AC; %this is the mid point in the RV

plottetramesh(H_mesh);
scatter3(OA(:,1),OA(:,2),OA(:,3),105,'filled');
scatter3(OB(:,1),OB(:,2),OB(:,3),105,'filled');
scatter3(OC(:,1),OC(:,2),OC(:,3),105,'filled');

rs = 1.0; %maxium distance between OC and furthest point
rv_nodes=[];
for i=1:size(OC,1)
    OCdummy=OC(i,:);
    [rv_node_candidates,OCXdist] = knnsearch(H_mesh.xyz,OCdummy,'K',1000,'Distance','euclidean');

    rv_nodes =horzcat(rv_nodes,rv_node_candidates(OCXdist<rs));
    if rem(i,100)==0
    disp(['i have searched',' ',num2str(i), ' ' ,'spheres']);
    rv_nodes=unique(rv_nodes);

    end
end

rv_nodes_xyz=H_mesh.xyz(rv_nodes',:);

figure()
plottetramesh(H_mesh,1);
alpha(0.3)
scatter3(rv_nodes_xyz(:,1),rv_nodes_xyz(:,2),rv_nodes_xyz(:,3),25,'filled')

%find tetra which are made up entirely of points from the RV region
closest_point_to_tetra_centre = knnsearch(H_mesh.xyz,H_mesh.tricentres);
rv_tri_ids=[];
for i= 1:size(rv_nodes,2)
    dummy1=rv_nodes(i);
    dummy=find(closest_point_to_tetra_centre==dummy1);
    if isempty(dummy)==0
        rv_tri_ids=vertcat(rv_tri_ids,dummy); %id of the tetra centres
    end
end


%find surface tris which are made up the RV epicardial surface

x1=H_mesh.xyz(H_mesh.rv(1,:),1)';
y1=H_mesh.xyz(H_mesh.rv(1,:),2)';
z1=H_mesh.xyz(H_mesh.rv(1,:),3)';

tri=H_mesh.face;
x=H_mesh.xyz(:,1);
y=H_mesh.xyz(:,2);
z=H_mesh.xyz(:,3);


TR=triangulation(H_mesh.epi,x,y,z);
C=circumcenter(TR);

figure()
plottetramesh(H_mesh,1)
hold on 
plot3(C(:,1),C(:,2),C(:,3),'r.','MarkerSize',10)
hold off
axis off


rv_face_tri_ids=[];
rv_epi_faces=[];

dummy1=rv_nodes;
dummy=knnsearch(C,H_mesh.xyz(rv_nodes,:),'K',10);
if isempty(dummy)==0
   % rv_face_tri_ids=vertcat(rv_tri_ids,dummy); %id of the tetra centres
    rv_epi_faces=vertcat(rv_epi_faces,H_mesh.epi(dummy,:));
end

H_rv=struct;
H_rv.xyz=H_mesh.xyz;
H_rv.face=rv_epi_faces;
H_rv.tri=H_mesh.tri(rv_tri_ids,:);
H_rv.face = vertcat(H_rv.face,H_mesh.rv);
figure()
plottetramesh(H_rv,1)

hplotMESH(H_rv)

RVVolume=sum( meshQuality( MeshFillHoles( H_rv ) , 'volume' ) );
TotVolume= sum( meshQuality( MeshFillHoles( H_mesh ) , 'volume' ) );



end