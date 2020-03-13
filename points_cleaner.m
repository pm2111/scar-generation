function [lge]=points_cleaner(H_mesh,cursor,lge)

%compute distances to all points
 distance= @(vec,vec1) ((vec(:,1)-vec1(:,1)).^2+(vec(:,2)-vec1(:,2)).^2+(vec(:,3)-vec1(:,3)).^2).^0.5;

 dists = distance(H_mesh.xyz,cursor.Position);
 
 lge(find(dists<1.5))=1;
end 
