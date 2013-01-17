function matrix=removePadding(matrix)

matrix(:,:,end)  = [];
matrix(:,end,:)  = [];
matrix(end,:,:)  = [];
matrix(:,:,1)    = [];
matrix(:,1,:)    = [];
matrix(1,:,:)    = [];