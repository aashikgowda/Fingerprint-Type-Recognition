% Read input image
% Filter with Gaussian Filter of size 7 and sigma 1
im = imread('0402.pgm');
im_f = imgaussfilt(im,1,'FilterSize',7);
% Compute gradient along x-axis and y-axis
% Compute gx^2, gy^2 and 2*gx*gy
% Apply Gaussian Filter of size 13 and sigma 3 to above values
[gx, gy]  = imgradientxy(im_f);
gx2 = gx.*gx ;gy2 = gy .* gy;
g2xy = 2.*gx .* gy;
Gx2 = imgaussfilt(gx2,3,'FilterSize',13);
Gy2 = imgaussfilt(gy2,3,'FilterSize',13);
Vx =  imgaussfilt(g2xy,3,'FilterSize',13);
Vy = Gx2 - Gy2;
Vy = imgaussfilt(Vy,3,'FilterSize',13);
% Determine theta value in degree and then convert to radian
theta_degree = 90 + (0.5 .* atan2d(Vx,Vy));
theta = theta_degree .* (pi/180);
% Code to display singularity points
    % Add zero padding
    theta_pad = padarray(theta,[1 1],'replicate');
    [row,col] = size(theta_pad);
    poncar(1:1:row,1:1:col) = 0;
    poncar_ind(1:1:row,1:1:col) = 0;
    poncar_tot(1:1:row,1:1:col) = 0;
    k = [0 -1;1 -1;1 0;1 1;0 1;-1 1;-1 0;-1 -1];
    l = [-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1;-1 0];
    % The entire loop calculates PI(i,j)for (i,j) pixel for its surrounding
    % 8 neighbours
    for i = 2:1:row-1
        for j = 2:1:col-1
            for o = 1:1:8
            poncar_ind(i,j) = theta_pad((i+k(o,1)),(j+k(o,2))) - theta_pad((i+l(o,1)),(j+l(o,2)));
               if poncar_ind(i,j) < -(pi/2)
                poncar_ind(i,j) = pi + poncar_ind(i,j);
            elseif poncar_ind(i,j) > pi/2
                poncar_ind(i,j) = pi - poncar_ind(i,j);
            elseif -(pi/2) <= poncar_ind(i,j) <= (pi/2)
                poncar_ind(i,j) = poncar_ind(i,j);
               end
               poncar_tot(i,j) = poncar_tot(i,j) + poncar_ind(i,j);
            end
            poncar(i,j) = round(poncar_tot(i,j)/(pi))/2;
        end
    end
    % Remove the zero padding
    poncarie = poncar(2:1:row-1,2:1:col-1);
    [l_xaxis,l_yaxis] = find(poncarie == -0.5);
    [l_row,l_col] = size(l_xaxis);
    [d_xaxis,d_yaxis] = find(poncarie == 0.5);
    [d_row,col2] = size(d_xaxis);
  

% THE FOLLOWING CODE FOR RIDGE ORIENTATION DISPLAY WAS ADPOTED FROM:
% Peter Kovesi  
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
% January 2005
    spacing = 7;
   
    figure(1),imshow(im)
    
    [rows, cols] = size(theta);
    len = 0.8*spacing;  % length of orientation lines
    % Subsample the orientation data according to the specified spacing
    theta_orient = theta(spacing:spacing:rows-spacing, ...
		      spacing:spacing:cols-spacing);

    xoff = len/2*cos(theta_orient);
    yoff = len/2*sin(theta_orient);    
    % Determine placement of orientation vectors
    [x,y] = meshgrid(spacing:spacing:cols-spacing, ...
		     spacing:spacing:rows-spacing);
    
    x = x-xoff;
    y = y-yoff;
    
    % Orientation vectors
    u = xoff*2;
    v = yoff*2;

    figure(2),quiver(x,y,u,v,0,'.','linewidth',1);
    
    axis equal, axis ij,  hold off
   
    figure(3),imshow(im)
    hold on
    quiver(x,y,u,v,0,'.','linewidth',0.7);
    
    figure(4), imshow(im)
    hold on
    for i =1:1:l_row
    plot(l_yaxis,l_xaxis,'bo');
    end
    hold on
    for i =1:1:d_row
    plot(d_yaxis,d_xaxis,'r^');
    end