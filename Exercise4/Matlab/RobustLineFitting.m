%% Load and plot points
load('points.mat','x','y');
figure;hold on;
plot(x,y,'kx');
axis equal

%% RANSAC parameters
% m is the number of data points
m=length(x);
% s is the size of the random sample
s=2;
% t is the inlier distance threshold
t=sqrt(3.84)*2;
% e is the expected outlier ratio
e=0.8;
% at least one random sample should be free from outliers with probability p 
p=0.999;
% required number of samples
N_estimated=log(1-p)/log(1-(1-e)^s);

inliers=[]
%% RANSAC loop
% First initialize some variables
N=inf;
sample_count=0;
max_inliers=0;
best_line=zeros(3,1);
x_inliers=[]
y_inliers=[]

while(N>sample_count)
    % Pick two random samples
    id1=ceil(m*rand(1));  % sample id 1
    id2=ceil(m*rand(1));  % sample id 2    
    if id1==id2 %if the same point is drawn twice, must draw again
        continue;
    end
    
    % Determine the line crossing the points with the cross product of the points (in homogeneous coordinates).
    % Also normalize the line by dividing each element by sqrt(a^2+b^2), where a and b are the line coefficients
    
    %%-your-code-starts-here%%
    a=[];
    d=[];
    a1=x(id1);a2=y(id1);
    a=[a,a1];a=[a,a2];a=[a,1];
    
    d1=x(id2);d2=y(id2);
    d=[d, d1]; d=[d, d2]; d=[d, 1];
    
    cross_product = cross(a,d);
    
    l=cross_product./sqrt((cross_product(1))^2+(cross_product(2))^2);

    %%-your-code-ends-here%%
    
    % Determine the inliers by finding the indices for the line and data
    % point dot products (absolute value) which are less than inlier threshold.
    
    %%-your-code-starts-here%%
    
    for i=1:m
        pt=[x(i) y(i) 1];
        distance1=abs(dot(l,pt))
        if distance1<t
            inliers=[inliers i];
        end   
    end
    %%-your-code-ends-here%%
    
    % keep the hypothesis giving most inliers so far
    inlier_count=length(inliers);
    if inlier_count>max_inliers
        best_line=l(:);
    end
    
    % update the estimate of the outlier ratio
    e=1-inlier_count/m;
    % update the estimate for the required number of samples
    N=log(1-p)/log(1-(1-e)^s);
    
    sample_count=sample_count+1;
end

% Least squares fitting to the inliers of the best hypothesis, i.e.
% find the inliers similarly as above but this time for the best line.

%%-your-code-starts-here%%
inliers=[];
for i=1:m
    pt2=[x(i) y(i) 1];
    distance2=abs(dot(best_line,pt2));
    if distance2<t
        inliers=[inliers i];
    end
end
x_inliers=x(inliers);
y_inliers=y(inliers);
%%-your-code-ends-here%%

% Fit a line to the given points (non-homogeneous)
l=linefitlsq(x_inliers, y_inliers);

% plot the resulting line and the inliers
k=-l(1)/l(2);
b=-l(3)/l(2);
plot(1:100,k*[1:100]+b,'m-');
plot(x(inliers),y(inliers),'ro');
