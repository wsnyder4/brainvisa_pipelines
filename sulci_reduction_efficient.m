%% read each subject OFC image as a 3D block in a 4D matrix holding all subjects
% identify our sulci of interst, convert to point coords 
% reverse x coords of $ hemisphere, convert point coord lists to point clouds

filelist = dir("/Users/willsnyder/Downloads/images/itksnap*.nii");
prefix = "/Users/willsnyder/Downloads/images/";
clouds.all = cell(1, 2*length(filelist));
for i = 1:length(filelist)
    disp(i)
    fname = filelist(i).name;
    %folder = fname(1:6);
    fname = char(prefix + fname);
    V = spm_vol(fname);
    [P, XYZ] = spm_read_vols(V);
    left = P(1:(end/2), :, :);
    right = P((end/2+1:end), : , :);
    [xl, yl, zl ] = ind2sub(size(left), find(left == 3));
    [xr, yr, zr ] = ind2sub(size(right),find(right == 3));
    
    left_coords = [xl yl zl];
    right_coords = [-xr yr zr];
    
    ptcloudl = pointCloud(left_coords);
    ptcloudr = pointCloud(right_coords);
    %disp((2*(i-1))+1)
    %disp(2*i)
    
    clouds.all{(2*(i-1))+1} = {ptcloudl};
    clouds.all{2*i} = {ptcloudr};
end




%% use parfor loop to, pairwise, compute RMSE, to be stored in pre-allocated matrix
tic;
distance_mat = zeros(length(clouds.all),length(clouds.all));
for i = 1:length(clouds.all)
    disp(i)
    cloudi = clouds.all(i);
    cloudi = cloudi{1,1};
    if class(cloudi) ~= "pointCloud"
        cloudi = cloudi{1,1};
    end
    for j = 1:length(clouds.all)
        cloudj = clouds.all(j);
        cloudj = cloudj{1,1};
        if class(cloudj) ~= "pointCloud"
            cloudj = cloudj{1,1};
        end
        
        
        if i <= j
            [tform, movingreg, rmse] = pcregistericp(cloudi, cloudj);
            distance_mat(i,j) = rmse;
            distance_mat(j,i) = rmse;
        end
        
        if i==7 && j ==8
            cloudy = movingreg;
            cloudz = cloudj;
        end
        
    end
end
toc;
ptconcat = pcmerge(cloudy, cloudz,1);
figure; pcshow(ptconcat); figure; pcshow(cloudy); figure; pcshow(cloudz)

%% Perform isomap learning on RMSE matrix, retain primary axis

vals = zeros(1,16);
for k = 5:20
    [mappedX, mapping] = isomap(distance_mat, 1, k);

    if length(mappedX) == size(distance_mat,1)
        disp("Mapping Successful")
    else
        disp("Error in nearest neighbor calculations. Please select a different k")
    end
    vals(k) = mapping.val;
end
plot(vals)

[~,ind] = min(vals);

[mappedX, mapping] = isomap(distance_mat, 1, ind);

    if length(mappedX) == size(distance_mat,1)
        disp("Mapping Successful")
    else
        disp("Error in nearest neighbor calculations. Please select a different k")
    end

  %% other isomap approach
  
  residuals = zeros(1,20);
  for k = 5:20    
    [Y, R, E] = Isomap(distance_mat, 'k', k);
    residuals(k) = R(1);
  end
  plot(residuals)
    
  
  [~,k] = min(residuals)
  
  %%
  [Y, R, E] = Isomap(distance_mat, 'k', 10);
  
    
  list = Y.coords{1};
  
  
%% align all sulci to template, "neutral" sulcus
clouds.template = cell(1,2*length(filelist));

[~, neutral_loc ] = min(sum(distance_mat));
for i = 1:length(clouds.all)
    disp(i)
    cloudi = clouds.all(i);
    cloudi = cloudi{1,1};
    if class(cloudi) ~= "pointCloud"
        cloudi = cloudi{1,1};
    end
    
    %%change to not be hardcoded
    cloudj = clouds.all(neutral_loc);
    %cloudj = clouds.all(17);
    
    cloudj = cloudj{1,1};
    if class(cloudj) ~= "pointCloud"
        cloudj = cloudj{1,1};
    end
     
    [tform, movingreg, rmse] = pcregistericp(cloudi, cloudj);
    clouds.template{i} = {movingreg};
    
        
end
    
    
%% look at distribution along primary axis, set number of evenly spaced bins
%-100 to 150 is the range here -- length = 250, make 15 bins

%hist(mappedX)
list =  list';
imgs = zeros(size(P,1) + 20 ,size(P,2) + 20,size(P,3),9);
img_count = 1;
%for i = -100:25:150
for i = -1.5:.5:2.5
    disp(i)
    dist = power((power(abs(list-i), 2) /1000),-1);
    %dist = power((power(abs(mappedX-i), 2) /1000),-1);
    
    
    for j = 1:length(clouds.all)
        cloudj = clouds.template(j);
        cloudj = cloudj{1,1};
        if class(cloudj) ~= "pointCloud"
            cloudj = cloudj{1,1};
        end
        xyz = round(cloudj.Location);
        for s = 1:size(xyz,1)
            x = xyz(s,1);
            y = xyz(s,2);
            z = xyz(s,3);
            imgs(x,y,z,img_count) = imgs(x,y,z,img_count) + 1*dist(j);
        end
        
    end
    img_count = img_count + 1;
end



%% compute weighted average image of each bin
% the average isomap point in each bin will be the center, and weight of
% each image will be inversely proportional to the squared distance from
% each center. Each image must first be registered to the center before
% this is computed


% nearest neighbors approach


imgs = zeros(size(P,1) + 20 ,size(P,2) + 20,size(P,3),9);
img_count = 1;
for i = -1.5:.5:2.5
    
    dist = power((power(abs(list-i), 2) /1000),-1);
    [~ , locs] = maxk(dist,20);
    
    for j = 1:length(locs)
        cloudj = clouds.template(locs(j));
        cloudj = cloudj{1,1};
        if class(cloudj) ~= "pointCloud"
            cloudj = cloudj{1,1};
        end
        xyz = round(cloudj.Location);
        for s = 1:size(xyz,1)
            x = xyz(s,1);
            y = xyz(s,2);
            
            z = xyz(s,3);
            imgs(x,y,z,img_count) = imgs(x,y,z,img_count) + 1;
        end
        
    end
    img_count = img_count + 1;
end



%% add the highest threshold that does not allow for holes, save averaged images, convert to mesh


imgs_gauss_thresh = imgs;
for i = 1:size(imgs,4)
    disp(i)
    imgs_gauss_thresh(:,:,:,i) = imgaussfilt3(imgs_gauss_thresh(:,:,:,i),2);
    %imgs_gauss_thresh(:,:,:,i) = imgaussfilt3(imgs_gauss_thresh(:,:,:,i));
    %icatb_write_nifti_data(char("moving_avg_sulc" + num2str(i) + ".nii"),V, imgs_gauss_thresh(:,:,:,i))
    
end

%% testing

 for j = flip(.1:.1:7)
        disp("j")
        disp(j)
        im_thresh = imgs_gauss(:,:,:,i) > j;
        invert = imcomplement(im_thresh);
        conn = bwconncomp(invert);
        disp("conn comps")
        disp(conn.NumObjects)
 end
    
