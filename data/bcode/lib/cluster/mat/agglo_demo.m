function agglo_demo()

%Generating 100 random features of dimension 10, 1 feature per column 
  F=rand(10,100);


cl_nb=20;%maximum number of clusters
thres=-10;%minimum similarity between features (negative euclidean distance)

%set threshold to -1000000 if a given cl_nb is required
%set cl_nb to 0 if similarity threshold is more important

[tr ass]=clrnncagglo('agglom',F,[cl_nb thres])
%tr is a clustering trace and ass is a cluster assignment
