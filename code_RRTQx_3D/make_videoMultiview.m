% The MIT License (MIT)
%
% Copyright (c) January, 2014 michael otte
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.


% this reads a series of saved files to help visualize search progress vs
% run time
clear all
close all

fig = figure(1);
fig.Position = [200 150 1600 1200]

%costDim = 6; % this is the dimension of the nodes files that has cost
%costDim = 5; % this is the dimension of the nodes files that has cost
costDim = 5; % this is the dimension of the nodes files that has cost

dir = 'temp/';
%dir = 'temp_dynamic_2d/'
%dir = 'temp_2d_forest/';
%dir = 'Video_RRTMine_0.1/';
file_ctr = 1;
max_file_ctr = 177; %298, 314

start_move_at_ctr = 30;

minXval = -22;
minYval = -22;
minZval = -18;
maxXval = 22;
maxYval = 22;
maxZval = 2;
tickInterval = 10;

contourGranularity = 2.5; % x, y granularity for cost contour plot
countorLevels = 0:2:200; % the countour levels
%countorLevels = 0:3:300; % the countour levels

%sensor_radius = 20/contourGranularity;
%sensor_radius = 10/contourGranularity;
sensor_radius = 5/contourGranularity;

% make a more accurate colormap, where costs above goal are grey
mycolormapsmall = autumn;
mycolormap = zeros(length(countorLevels),3);
originalsamplevals = 0:length(mycolormapsmall)-1;
myresamplevalues = (0:length(countorLevels)-1)/(length(countorLevels)-1)*(length(mycolormapsmall)-1);
for c = 1:3
  mycolormap(:,c) = interp1(originalsamplevals, mycolormapsmall(:,c), myresamplevalues, 'pchip');
end




videoFill = true; % only uses cost in the negative y direction in center of image 
                  % to fill in missing cost so that don't get wierd stuff that 
                  % looks funnly due to contour interpolation

% open the video file
writerObj = VideoWriter( [dir 'Movie.avi']);
writerObj.FrameRate = 10;
open(writerObj);

while exist([dir 'robotMovePath_1_' num2str(file_ctr) '.txt'], 'file')  && file_ctr < max_file_ctr
    display(file_ctr)
  
    x_start = 0;
    y_start = -5;
    z_start = 0;
    
    EdgeData = load([dir 'edges_1_' num2str(file_ctr) '.txt']);
    
    if isempty(EdgeData)
      i = 0
      raw_x = [];
      raw_y = [];
      raw_z = [];
    else
    
    i = size(EdgeData,1);
    raw_x = EdgeData(1:i,1);
    raw_y = EdgeData(1:i,2);
    raw_z = EdgeData(1:i,3);
    end
    
    x = nan(i*1.5, 1);
    y = nan(i*1.5, 1);
    z = nan(i*1.5, 1);
    
    x(1:3:end-2) = raw_x(1:2:end-1);
    x(2:3:end-1) = raw_x(2:2:end);
    
    y(1:3:end-2) = raw_y(1:2:end-1);
    y(2:3:end-1) = raw_y(2:2:end);
    
    z(1:3:end-2) = raw_z(1:2:end-1);
    z(2:3:end-1) = raw_z(1:2:end-1);
    
    
%     % do not plot edges with inf cost
%     infinds = find(isinf(z));
%     x(infinds) = nan;
%     y(infinds) = nan;
    
    %figure(1)
    %plot3(x,y,z)
    
    
    NodeData = load([dir 'nodes_1_' num2str(file_ctr) '.txt']);
    
    j = size(NodeData,1);
    node_x = NodeData(1:j,1);
    node_y = NodeData(1:j,2);
    node_z = NodeData(1:j,3);
    node_cost = NodeData(1:j,costDim);

    if exist([dir 'CnodmoveGoales_1_' num2str(file_ctr) '.txt'], 'file') 
        CNodeData = load([dir 'Cnodes_1_' num2str(file_ctr) '.txt']);
    else
        CNodeData = [];
    end
    
    j = size(CNodeData,1);
    if j > 0
        Cnode_x = CNodeData(1:j,1);
        Cnode_y = CNodeData(1:j,2);
        Cnode_z = CNodeData(1:j,3);
    else
        Cnode_x = [nan nan];
        Cnode_y = [nan nan];
        Cnode_z = [nan nan];
    end
    
    if isfile([dir 'obstacles_3_' num2str(file_ctr) '.txt'])
        ObstacleData = load([dir 'originalObs_1_' num2str(file_ctr) '.txt']);
    else
        ObstacleData = load([dir 'originalObs_1_' num2str(file_ctr) '.txt']);
    end
    if ~isempty(ObstacleData)
        obs_x = ObstacleData(:,1);
        obs_y = ObstacleData(:,2);
        obs_z = ObstacleData(:,3);
        obs_r = ObstacleData(:,4);
    else
        obs_x = [];
        obs_y = [];
        obs_z = [];
        obs_r = [];
    end

%     if isfile([dir 'obstacles_3_' num2str(file_ctr) '.txt'])
%         ObstacleData = load([dir 'originalObs_1_1.txt']);
%     else
%         ObstacleData = load([dir 'originalObs_1_1.txt']);
%     end
       if isfile([dir 'obstacles_3_' num2str(file_ctr) '.txt'])
        ObstacleData = load([dir 'obstacles_1_' num2str(file_ctr) '.txt']);
    else
        ObstacleData = load([dir 'obstacles_1_' num2str(file_ctr) '.txt']);
    end
    if ~isempty(ObstacleData)
        obs_x = ObstacleData(:,1);
        obs_y = ObstacleData(:,2);
        obs_z = ObstacleData(:,3);
        obs_r = ObstacleData(:,4);
    else
        obs_x = [];
        obs_y = [];
        obs_z = [];
        obs_r = [];
    end
    
    if (isfile([dir 'robotMovePath_3_' num2str(file_ctr) '.txt']))
        MoveData = load([dir 'robotMovePath_1_' num2str(file_ctr) '.txt']);
    else
        MoveData = load([dir 'robotMovePath_1_' num2str(file_ctr) '.txt']);
    end
    move_x = MoveData(:,1);
    move_y = MoveData(:,2);
    move_z = MoveData(:,3);
    move_theta = 0; % MoveData(:,4);
        
    if (isfile([dir 'robotMovePath_2_' num2str(file_ctr) '.txt']))
        MoveData2 = load([dir 'robotMovePath_2_' num2str(file_ctr) '.txt']);
    else
        MoveData2 = load([dir 'robotMovePath_1_' num2str(file_ctr) '.txt']);
    end
    move_x2 = MoveData2(:,1);
    move_y2 = MoveData2(:,2);
    move_z2 = MoveData2(:,3);
    move_theta2 = 0; % MoveData(:,4);
    
    if (isfile([dir 'robotMovePath_3_' num2str(file_ctr) '.txt']))
        MoveData3 = load([dir 'robotMovePath_3_' num2str(file_ctr) '.txt']);
    else
        MoveData3 = load([dir 'robotMovePath_1_' num2str(file_ctr) '.txt']);
    end
    move_x3 = MoveData3(:,1);
    move_y3 = MoveData3(:,2);
    move_z3 = MoveData3(:,3);   
    move_theta3 = 0; % MoveData(:,4);
        
    if (isfile([dir 'robotMovePath_4_' num2str(file_ctr) '.txt']))
        MoveData4 = load([dir 'robotMovePath_4_' num2str(file_ctr) '.txt']);
    else
        MoveData4 = load([dir 'robotMovePath_1_' num2str(file_ctr) '.txt']);
    end
    move_x4 = MoveData4(:,1);
    move_y4 = MoveData4(:,2);
    move_z4 = MoveData4(:,3);
    move_theta4 = 0; % MoveData(:,4);
    
    if (isfile([dir 'path_3_' num2str(file_ctr) '.txt']))
        PathData = load([dir 'path_1_' num2str(file_ctr) '.txt']);
    else
        PathData = load([dir 'path_1_' num2str(file_ctr) '.txt']);
    end
    path_x = PathData(:,1);
    path_y = PathData(:,2);
    path_z = PathData(:,3);
    if (isfile([dir 'path_2_' num2str(file_ctr) '.txt']))
        PathData2 = load([dir 'path_2_' num2str(file_ctr) '.txt']);
    else
        PathData2 = load([dir 'path_1_' num2str(file_ctr) '.txt']);
    end
    path_x2 = PathData2(:,1);
    path_y2 = PathData2(:,2);
    path_z2 = PathData2(:,3);
    if (isfile([dir 'path_3_' num2str(file_ctr) '.txt']))
        PathData3 = load([dir 'path_3_' num2str(file_ctr) '.txt']);
    else
        PathData3 = load([dir 'path_1_' num2str(file_ctr) '.txt']);
    end
    path_x3 = PathData3(:,1);
    path_y3 = PathData3(:,2);
    path_z3 = PathData3(:,3);
    if (isfile([dir 'path_4_' num2str(file_ctr) '.txt']))
        PathData4 = load([dir 'path_4_' num2str(file_ctr) '.txt']);
    else
        PathData4 = load([dir 'path_1_' num2str(file_ctr) '.txt']);
    end
    path_x4 = PathData4(:,1);
    path_y4 = PathData4(:,2);
    path_z4 = PathData4(:,3);
    
    if (file_ctr > start_move_at_ctr)
        %Uncomment with 4
        MoveData2 = load([dir 'robotMovePath_2_' num2str(file_ctr) '.txt']);
        move_x2 = MoveData2(:,1);
        move_y2 = MoveData2(:,2);
        move_z2 = MoveData2(:,3);
    
        MoveData3 = load([dir 'robotMovePath_3_' num2str(file_ctr) '.txt']);
        move_x3 = MoveData3(:,1);
        move_y3 = MoveData3(:,2);
        move_z3 = MoveData3(:,3);
    
        MoveData4 = load([dir 'robotMovePath_4_' num2str(file_ctr) '.txt']);
        move_x4 = MoveData4(:,1);
        move_y4 = MoveData4(:,2);
        move_z4 = MoveData4(:,3);
    end
        
    
    if exist([dir 'obsNodes.txt'], 'file') 
      OBSNodesData = load([dir 'obsNodes.txt']);
      OBSNodes_x = OBSNodesData(:,1);
      OBSNodes_y = OBSNodesData(:,2);
      OBSNodes_z = OBSNodesData(:,3);
    else
      OBSNodes_x = [];
      OBSNodes_y = [];
      OBSNodes_z = [];
    end
    
    
   if exist([dir 'obsNodesNeighbors.txt'], 'file') 
      OBSNodesNData = load([dir 'obsNodesNeighbors.txt']);
      OBSNodesN_x = OBSNodesNData(:,1);
      OBSNodesN_y = OBSNodesNData(:,2);
      OBSNodesN_z = OBSNodesNData(:,3);
    else
      OBSNodesN_x = [];
      OBSNodesN_y = [];
      OBSNodesN_z = [];
    end
    
%         figure(3)
%         clf
%         plot(node_x, node_y, '.r')
%         hold on
%         plot(Cnode_x, Cnode_y, 'xk')
%         plot(x, y)
%         plot(obs_x,obs_y, 'g', 'LineWidth',2)
%         plot(path_x,path_y, 'c', 'LineWidth',2)
%         plot(move_x,move_y, 'k', 'LineWidth',2)
%         hold off
%         axis([-50 50 -50 50])
    
    
    
    % make a contour plot, need to grid up data into Z, an m x n matrix
    % we'll set the cost value of Z(i,j) to be the average value of all
    % nodes that exist within the transformed grid [i i+1] X [j j+1]
    % note that we need to define these well with respect to obstacles
    
    Xs = minXval:contourGranularity:(maxXval-1);
    Ys = minYval:contourGranularity:(maxYval-1);
    Zs = minZval:contourGranularity:(maxZval-1);
    
    
    Z = zeros(length(Ys), length(Xs), length(Zs));
    Zmin = inf(length(Ys), length(Xs), length(Zs));
    Counts = zeros(length(Ys), length(Xs), length(Zs)); % remember the number of points in each grid
    
    for v = 1:length(node_cost)
        j = max(min(floor((node_x(v) - minXval)/contourGranularity)+1, size(Z,2)),1);
        i = max(min(floor((node_y(v) - minYval)/contourGranularity)+1, size(Z,1)),1);
        k = max(min(floor((node_z(v) - minZval)/contourGranularity)+1, size(Z,3)),1);
        
        if isinf(node_cost(v))
            % nodes that are furhter than path-legnth from the goal
            % never have thier initial inf value removed (since
            % the wavefeont has not expanded to them yet)
          continue
        end
        
        Z(i,j,k) = Z(i,j,k) + node_cost(v);
        Counts(i,j,k) = Counts(i,j,k) + 1;
        
        Zmin(i,j,k) = min(Zmin(i,j,k), node_cost(v));
    end
    % now take the average
    Z(Z ~= 0) = Z(Z ~= 0)./Counts(Z ~= 0);
    
    % now, if there are cells do not have any values, we need to put
    % something in them, we'll average over non-zero neighbors
    % but do nothing if path to goal does not exist
    jjj = 1;
    while sum(Z(:)==0) > 0
     
        jjj = jjj + 1;
        
        if sum(Z(:)~=0) < 5 || jjj > 3 % max(size(Z))
            % whole thing is unpopulated
            Z(Z(:)==0) = countorLevels(end);
            Zmin(Z(:)==0) = countorLevels(end);
            break
        end
        
        [yZeroInds, xZeroInds] = find(Z == 0);
        zZeroInds = zeros(length(xZeroInds));
        for k = 1:length(xZeroInds)
            zZeroInds(k) = (floor(xZeroInds(k)/length(Xs)) + 1);
            xZeroInds(k) = mod(xZeroInds(k), length(Xs));
        if (xZeroInds(k) == 0)
            xZeroInds(k) = length(Xs);
            zZeroInds(k) = zZeroInds(k) - 1;
        end
        end
        
        dZ = zeros(size(Z));
        dZmin = zeros(size(Zmin));
        for k = 1:length(xZeroInds)
            xZind = xZeroInds(k);
            yZind = yZeroInds(k);
            zZind = zZeroInds(k);
            
            % find the "3 x 3" sub matrix around (yZind, xZind) in Z
            % note that it can be smaller if (yZind, xZind) is on the border
            minxZind = max(xZind - 1, 1);
            maxxZind = min(xZind + 1, size(Z, 2));
            minyZind = max(yZind - 1, 1);
            maxyZind = min(yZind + 1, size(Z, 1));
            minzZind = max(zZind - 1, 1);
            maxzZind = min(zZind + 1, size(Z, 3));
            
            if (videoFill == true) && (length(Ys)/2-2 < yZind) && (yZind < length(Ys)/2+2)  
                % use 2 X 3 (equal to and below) so that don't get cost
                % through wall obstacle due to contour interpolation
                
                minyZind = yZind;
            end
            
            subZ = Z(minyZind:maxyZind, minxZind:maxxZind, minzZind:maxzZind);
            testSize = size(dZ);
            if  sum(subZ(:) ~= 0) > 1
                dZ(yZind, xZind, zZind) = sum(subZ(:))/sum(subZ(:) ~= 0);
                dZmin(yZind, xZind, zZind) = min(subZ(subZ(:) ~= 0));
            end
        end
        Z = Z+dZ;
        Zmin(isinf(Zmin) & dZmin ~= 0) = dZmin(isinf(Zmin) & dZmin ~= 0);
    end
     
    % transform nodes so we can plot them on the contour plot
    c_node_x = (node_x-minXval)/contourGranularity+.5;
    c_node_y = (node_y-minYval)/contourGranularity+.5;
    c_node_z = (node_z-minZval)/contourGranularity+.5;
    
    % transform path so we can plot it on the contour plot
    c_path_x = (path_x-minXval)/contourGranularity+.5;
    c_path_y = (path_y-minYval)/contourGranularity+.5;
    c_path_z = (path_z-minZval)/contourGranularity+.5;
    
    if (~isnan(path_x4))
        c_path_x2 = (path_x2-minXval)/contourGranularity+.5;
        c_path_y2 = (path_y2-minYval)/contourGranularity+.5;
        c_path_z2 = (path_z2-minZval)/contourGranularity+.5;
        
        c_path_x3 = (path_x4-minXval)/contourGranularity+.5;
        c_path_y3 = (path_y4-minYval)/contourGranularity+.5;
        c_path_z3 = (path_z4-minZval)/contourGranularity+.5;
        
        c_path_x4 = (path_x4-minXval)/contourGranularity+.5;
        c_path_y4 = (path_y4-minYval)/contourGranularity+.5;
        c_path_z4 = (path_z4-minZval)/contourGranularity+.5;
    end
    
    % transform move path so we can plot it on the contour plot
    c_move_x = (move_x-minXval)/contourGranularity+.5;
    c_move_y = (move_y-minYval)/contourGranularity+.5;
    c_move_z = (move_z-minZval)/contourGranularity+.5;
    if (~isnan(move_x4))
    %uncomment with 4
    c_move_x2 = (move_x2-minXval)/contourGranularity+.5;
    c_move_y2 = (move_y2-minYval)/contourGranularity+.5;
    c_move_z2 = (move_z2-minZval)/contourGranularity+.5;
    
    c_move_x3 = (move_x3-minXval)/contourGranularity+.5;
    c_move_y3 = (move_y3-minYval)/contourGranularity+.5;
    c_move_z3 = (move_z3-minZval)/contourGranularity+.5;
    
    c_move_x4 = (move_x4-minXval)/contourGranularity+.5;
    c_move_y4 = (move_y4-minYval)/contourGranularity+.5;
    c_move_z4 = (move_z4-minZval)/contourGranularity+.5;
    end
    c_move_theta = move_theta - pi/2;
    
    % transform edges
    c_x = (x-minXval)/contourGranularity+.5;
    c_y = (y-minYval)/contourGranularity+.5;
    c_z = (z-minZval)/contourGranularity+.5;
    
    % transform obstacles
    c_obs_x = (obs_x-minXval)/contourGranularity+.5;
    c_obs_y = (obs_y-minYval)/contourGranularity+.5;
    c_obs_z = (obs_z-minZval)/contourGranularity+.5;
    
    % calculate ticks
    theXticks = -50:tickInterval:maxXval;
    theYticks = -50:tickInterval:maxXval;
    theZticks = -50:tickInterval:maxZval;
    c_theXticks = (theXticks-minXval)/contourGranularity+.5;
    c_theYticks = (theYticks-minYval)/contourGranularity+.5;
    c_theZticks = (theZticks-minZval)/contourGranularity+.5;
   
    c_OBSNodes_x = (OBSNodes_x-minXval)/contourGranularity+.5;
    c_OBSNodes_y = (OBSNodes_y-minXval)/contourGranularity+.5;
    c_OBSNodes_z = (OBSNodes_z-minZval)/contourGranularity+.5;
    
    c_OBSNodesN_x = (OBSNodesN_x-minXval)/contourGranularity+.5;
    c_OBSNodesN_y = (OBSNodesN_y-minXval)/contourGranularity+.5;
    c_OBSNodesN_z = (OBSNodesN_z-minZval)/contourGranularity+.5;
    
    % make path start at robot pose and not move goal 
    %c_path_x(1) = c_move_x(end);
    %c_path_y(1) = c_move_y(end);
    
     
    %     figure(4)
    %     clf
    %     contourf(Z,countorLevels)
    %     hold on
    %     plot(c_path_x, c_path_y, 'k', 'LineWidth',3)
    %     plot(c_path_x, c_path_y, 'w', 'LineWidth',1)
    %     hold off
    %
    %     axis([1, size(Z,1), 1, size(Z,2)-1])
    %     set(gca,'XTick',c_theXticks)
    %     set(gca,'XTickLabel', theXticks)
    %     set(gca,'YTick',c_theYticks)
    %     set(gca,'YTickLabel', theYticks)
    
   
    % extract polygon obstacle inds
    if ~isempty(c_obs_x)
        naninds = find(isnan(c_obs_x));
        polyStartInds = [1; (naninds(1:end-1))+1];
        polyEndEnds = naninds - 1;
    else
        polyStartInds = [];
        polyEndEnds = [];
    end
    
    
    % make costs above robot cost color grey
    robot_Zind_x = floor(c_move_x(end))+1;
    robot_Zind_y = floor(c_move_y(end))+1;
    robot_Zind_z = floor(c_move_z(end))+1;
    
    % but keep cells near robot original values
    robot_minxZind = max(robot_Zind_x - 1, 1);
    robot_maxxZind = min(robot_Zind_x + 1, size(Z, 2));
    robot_minyZind = max(robot_Zind_y - 1, 1);
    robot_maxyZind = min(robot_Zind_y + 1, size(Z, 1));
    robot_minzZind = max(robot_Zind_z - 1, 1);
    robot_maxzZind = min(robot_Zind_z + 1, size(Z, 3));
    robotSamplePatch = Zmin(robot_minyZind:robot_maxyZind, robot_minxZind:robot_maxxZind, robot_minzZind:robot_maxzZind);

    robotSampleVal = mean(robotSamplePatch(~isinf(robotSamplePatch)));
    %robotSampleVal = Zmin(robot_Zind_y, robot_Zind_x);
    
    plotContourValIndRobotVal = find(countorLevels >= robotSampleVal, 1,'first')+1; 
    if isempty(plotContourValIndRobotVal) || plotContourValIndRobotVal > length(countorLevels)
      plotContourValIndRobotVal = length(countorLevels);
    end
    plotContourValIndAboveRobotVal = min(plotContourValIndRobotVal + 1, length(countorLevels));
    %dividingCost = (countorLevels(plotContourValIndRobotVal) + countorLevels(plotContourValIndAboveRobotVal))/2;
    dividingCost = countorLevels(plotContourValIndRobotVal);
    
    Zmin(Zmin > dividingCost) = Inf;
    Zmin(robot_minyZind:robot_maxyZind, robot_minxZind:robot_maxxZind, robot_minzZind:robot_maxzZind) = robotSamplePatch;
    
    tempcolormap = mycolormap;
    tempcolormap(plotContourValIndAboveRobotVal:end,:) = 1;
    
    % this keeps the colormap constant vs time
    Zmin(1,length(Xs)+1) = countorLevels(end);
   
    
%     figure(3)
%     imagesc(Zmin)
%     axis xy

    
    clf
    subplot(2,2,1);
    
    %clf
    %colormap(tempcolormap)
    %contourf(Zmin,countorLevels, 'EdgeColor', 'none')
    hold on
    %view(135,90)
    %plot3(c_x, c_y, c_z, 'Color',[0.3,0.3,0.3],'LineWidth',.5)
    plot3(c_node_x, c_node_y, c_node_z, '.', 'Color',[0.3,0.3,0.3],'MarkerSize',.5)
    plot3(c_move_x, c_move_y, c_move_z, 'k', 'LineWidth',3)
    plot3(c_move_x, c_move_y, c_move_z, 'r', 'LineWidth',3)
    plot3(c_move_x2, c_move_y2, c_move_z2, 'k', 'LineWidth',3)
    plot3(c_move_x2, c_move_y2, c_move_z2, 'b', 'LineWidth',3)
    plot3(c_move_x3, c_move_y3, c_move_z3, 'k', 'LineWidth',3)
    %plot3(c_move_x3, c_move_y3, c_move_z3, 'b', 'LineWidth',3)
    plot3(c_move_x4, c_move_y4, c_move_z4, 'k', 'LineWidth',3)
    %plot3(c_move_x4, c_move_y4, c_move_z4, 'b', 'LineWidth',3)
    %plot(c_path_x, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x, c_path_y, 'w', 'LineWidth',1)
    %plot(c_path_x2, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x2, c_path_y, 'w', 'LineWidth',1)
    %plot(c_path_x3, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x3, c_path_y, 'w', 'LineWidth',1)
    %plot(c_path_x4, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x4, c_path_y, 'w', 'LineWidth',1)
    %contour(Zmin,countorLevels(1:(maxPlotcontourValindPlusOne-1)), 'k')
    %for p = 1:length(polyStartInds)
    %   poly_inds = polyStartInds(p):polyEndEnds(p);
    %   patch(c_obs_x(poly_inds), c_obs_y(poly_inds), 'k') 
    %end
    %plot3(c_obs_x,c_obs_y, c_obs_z, 'w', 'LineWidth',1)
    if start_move_at_ctr > file_ctr
      plot3(c_node_x(1), c_node_y(1), c_node_z(1), 'sw', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',8)
    else
      plot3(c_path_x(end), c_path_y(end), c_node_z(end), 'sw', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',8)
    end
    plot3(c_move_x(end), c_move_y(end), c_move_z(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    if (~isnan(move_x4))
    plot3(c_move_x2(end), c_move_y2(end), c_move_z2(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    plot3(c_move_x3(end), c_move_y3(end), c_move_z3(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    plot3(c_move_x4(end), c_move_y4(end), c_move_z4(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    end
    buildingOrigin9 = [((-17.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim = [1.8, 1.8, 3];
    plotcube(dim, buildingOrigin9, 1, [.3 .3 .3]);
    buildingOrigin10 = [((-17.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim2 = [1.8, 1.8, 7];
    plotcube(dim2, buildingOrigin10, 1, [.3 .3 .3]);
    buildingOrigin11 = [((-17.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim3 = [1.8, 1.8, 4];
    plotcube(dim3, buildingOrigin11, 1, [.3 .3 .3]);
    buildingOrigin12 = [((-17.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim4 = [1.8, 1.8, 3];
    plotcube(dim4, buildingOrigin12, 1, [.3 .3 .3]);
    buildingOrigin = [((-7.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim = [1.8, 1.8, 8];
    plotcube(dim, buildingOrigin, 1, [.3 .3 .3]);
    buildingOrigin2 = [((-7.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim2 = [1.8, 1.8, 12];
    plotcube(dim2, buildingOrigin2, 1, [.3 .3 .3]);
    buildingOrigin13 = [((-7.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim = [1.8, 1.8, 3];
    plotcube(dim, buildingOrigin13, 1, [.3 .3 .3]);
    buildingOrigin14 = [((-7.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim6 = [1.8, 1.8, 4];
    plotcube(dim6, buildingOrigin14, 1, [.3 .3 .3]);
    buildingOrigin3 = [((2.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim3 = [1.8, 1.8, 5];
    plotcube(dim3, buildingOrigin3, 1, [.3 .3 .3]);
    buildingOrigin4 = [((2.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim4 = [1.8, 1.8, 10];
    plotcube(dim4, buildingOrigin4, 1, [.3 .3 .3]);
    buildingOrigin5 = [((12.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim5 = [1.8, 1.8, 6];
    plotcube(dim5, buildingOrigin5, 1, [.3 .3 .3]);
    buildingOrigin6 = [((12.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim6 = [1.8, 1.8, 4];
    plotcube(dim6, buildingOrigin6, 1, [.3 .3 .3]);
    buildingOrigin7 = [((2.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim7 = [1.8, 1.8, 3];
    plotcube(dim7, buildingOrigin7, 1, [.3 .3 .3]);
    buildingOrigin8 = [((2.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim8 = [1.8, 1.8, 3];
    plotcube(dim8, buildingOrigin8, 1, [.3 .3 .3]);
    buildingOrigin15 = [((12.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim7 = [1.8, 1.8, 5];
    plotcube(dim7, buildingOrigin15, 1, [.3 .3 .3]);
    buildingOrigin16 = [((12.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim8 = [1.8, 1.8, 4];
    plotcube(dim8, buildingOrigin16, 1, [.3 .3 .3]);
    buildingFloor = [((-20.5 - minXval)/contourGranularity+.5), ((-20.5 - minYval)/contourGranularity+.5), ((-20- minZval)/contourGranularity+.5)];
    dim5 = [18, 18, .5];
    plotcube(dim5, buildingFloor, 1, [.3 1 .3]);
    
    if (size(ObstacleData,1) > 90)
        for k = 1:90
            obs_x = (ObstacleData(k,1)-minXval)/contourGranularity+.5;
            obs_y = (ObstacleData(k,2)-minYval)/contourGranularity+.5;
            obs_z = (ObstacleData(k,3)-minZval)/contourGranularity+.5;
            if (ObstacleData(k,4) > 0)
                [temp1, temp2, temp3] = sphere;
                temp1 = temp1 * (ObstacleData(k,4)/contourGranularity)*.55;
                temp2 = temp2 * (ObstacleData(k,4)/contourGranularity)*.55;
                temp3 = temp3 * (ObstacleData(k,4)/contourGranularity)*.55;
                if(ObstacleData(k,4) > 2.8)
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*1);
                else
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*.85);
                end
                %plot3(obs_x, obs_y, obs_z, 'o', 'LineWidth', 1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.7, 0.7, 0.7],'MarkerSize', 20*ObstacleData(k,4))
            end
        end
        k = size(ObstacleData,1);
        obs_x = (ObstacleData(k,1)-minXval)/contourGranularity+.5;
        obs_y = (ObstacleData(k,2)-minYval)/contourGranularity+.5;
        obs_z = (ObstacleData(k,3)-minZval)/contourGranularity+.5;
        if (ObstacleData(k,4) > 0)
            [temp1, temp2, temp3] = sphere;
            temp1 = temp1 * ObstacleData(k,4)/contourGranularity;
            temp2 = temp2 * ObstacleData(k,4)/contourGranularity;
            temp3 = temp3 * ObstacleData(k,4)/contourGranularity;
            %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*.85);
            %plot3(obs_x, obs_y, obs_z, 'o', 'LineWidth', 1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.7, 0.7, 0.7],'MarkerSize', 20*ObstacleData(k,4))
        end
    else
        for k = 1:size(ObstacleData,1)
            obs_x = (ObstacleData(k,1)-minXval)/contourGranularity+.45;
            obs_y = (ObstacleData(k,2)-minYval)/contourGranularity+.45;
            obs_z = (ObstacleData(k,3)-minZval)/contourGranularity+.45;
            if (ObstacleData(k,4) > 0)
                [temp1, temp2, temp3] = sphere;
                %if((ObstacleData(k,4) == 1.3) || (ObstacleData(k,4) == .78))
                %    ObstacleData(k,4) = ObstacleData(k,4)/.65;
                %end
                ObstacleData(k,4) = ObstacleData(k,4)/.75;
                temp1 = temp1 * ObstacleData(k,4)/contourGranularity*.8;
                temp2 = temp2 * ObstacleData(k,4)/contourGranularity*.8;
                temp3 = temp3 * ObstacleData(k,4)/contourGranularity*.8;
                if(ObstacleData(k,4) > 2.8)
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*1);
                else
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*1);
                end
                %plot3(obs_x, obs_y, obs_z, 'o', 'LineWidth', 1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.7, 0.7, 0.7],'MarkerSize',20*ObstacleData(k,4))
            end
        end
    end
%    plot(raw_x(2), raw_y(2), 'sw', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',8)
    plot3(c_OBSNodes_x, c_OBSNodes_y, c_OBSNodes_z,'.k')
    plot3(c_OBSNodesN_x, c_OBSNodesN_y, c_OBSNodesN_z, 'ok')
    %plotVehicleTheta(fig, [c_move_x(end) c_move_y(end)], c_move_theta(end), .5, 'k', sensor_radius, '--w')
    
    if length(c_path_x) < 2
        c_path_x = [c_path_x c_path_x];
        c_path_y = [c_path_y c_path_y];
        c_path_z = [c_path_z c_path_z];
    end
    
    if length(c_path_x2) < 2
        c_path_x = [c_path_x2 c_path_x2];
        c_path_y = [c_path_y2 c_path_y2];
        c_path_z = [c_path_z2 c_path_z2];
    end
    
    if length(c_path_x3) < 2
        c_path_x = [c_path_x3 c_path_x3];
        c_path_y = [c_path_y3 c_path_y3];
        c_path_z = [c_path_z3 c_path_z3];
    end
    
    if length(c_path_x4) < 2
        c_path_x = [c_path_x4 c_path_x4];
        c_path_y = [c_path_y4 c_path_y4];
        c_path_z = [c_path_z4 c_path_z4];
    end
    
    %plotVehicle(fig, [c_path_x(1) c_path_y(1)], [c_path_x(2) c_path_y(2)], .5, 'k', sensor_radius, '--w')

    
    hold off
    %axis([1, size(Z,1), 1, size(Z,2)-1])
    set(gca,'XTick',[])
    %set(gca,'XTickLabel', theXticks)
    set(gca,'YTick',[])
    %set(gca,'YTickLabel', theYticks)
    set(gca,'ZTick',[])
    %set(gca,'ZTickLabel', theZticks)
    %set(gca,'FontSize',20)
    %colorbar
    view(45,60)
    title('Side 1 View')
    
    %axis tight
    %set(gca,'nextplot','replacechildren');
    %set(gcf,'Renderer','zbuffer');
     subplot(2,2,2);
    
    %clf
    %colormap(tempcolormap)
    %contourf(Zmin,countorLevels, 'EdgeColor', 'none')
    hold on
    view(135,90)
    %plot3(c_x, c_y, c_z, 'Color',[0.3,0.3,0.3],'LineWidth',.5)
    plot3(c_node_x, c_node_y, c_node_z, '.', 'Color',[0.3,0.3,0.3],'MarkerSize',.5)
    plot3(c_move_x, c_move_y, c_move_z, 'k', 'LineWidth',3)
    plot3(c_move_x, c_move_y, c_move_z, 'r', 'LineWidth',3)
    plot3(c_move_x2, c_move_y2, c_move_z2, 'k', 'LineWidth',3)
    plot3(c_move_x2, c_move_y2, c_move_z2, 'b', 'LineWidth',3)
    plot3(c_move_x3, c_move_y3, c_move_z3, 'k', 'LineWidth',3)
    %plot3(c_move_x3, c_move_y3, c_move_z3, 'b', 'LineWidth',3)
    plot3(c_move_x4, c_move_y4, c_move_z4, 'k', 'LineWidth',3)
    %plot3(c_move_x4, c_move_y4, c_move_z4, 'b', 'LineWidth',3)
    %plot(c_path_x, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x, c_path_y, 'w', 'LineWidth',1)
    %plot(c_path_x2, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x2, c_path_y, 'w', 'LineWidth',1)
    %plot(c_path_x3, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x3, c_path_y, 'w', 'LineWidth',1)
    %plot(c_path_x4, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x4, c_path_y, 'w', 'LineWidth',1)
    %contour(Zmin,countorLevels(1:(maxPlotcontourValindPlusOne-1)), 'k')
    %for p = 1:length(polyStartInds)
    %   poly_inds = polyStartInds(p):polyEndEnds(p);
    %   patch(c_obs_x(poly_inds), c_obs_y(poly_inds), 'k') 
    %end
    %plot3(c_obs_x,c_obs_y, c_obs_z, 'w', 'LineWidth',1)
    if start_move_at_ctr > file_ctr
      plot3(c_node_x(1), c_node_y(1), c_node_z(1), 'sw', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',8)
    else
      plot3(c_path_x(end), c_path_y(end), c_node_z(end), 'sw', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',8)
    end
    plot3(c_move_x(end), c_move_y(end), c_move_z(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    if (~isnan(move_x4))
    plot3(c_move_x2(end), c_move_y2(end), c_move_z2(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    plot3(c_move_x3(end), c_move_y3(end), c_move_z3(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    plot3(c_move_x4(end), c_move_y4(end), c_move_z4(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    end
    buildingOrigin9 = [((-17.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim = [1.8, 1.8, 3];
    plotcube(dim, buildingOrigin9, 1, [.3 .3 .3]);
    buildingOrigin10 = [((-17.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim2 = [1.8, 1.8, 7];
    plotcube(dim2, buildingOrigin10, 1, [.3 .3 .3]);
    buildingOrigin11 = [((-17.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim3 = [1.8, 1.8, 4];
    plotcube(dim3, buildingOrigin11, 1, [.3 .3 .3]);
    buildingOrigin12 = [((-17.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim4 = [1.8, 1.8, 3];
    plotcube(dim4, buildingOrigin12, 1, [.3 .3 .3]);
    buildingOrigin = [((-7.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim = [1.8, 1.8, 8];
    plotcube(dim, buildingOrigin, 1, [.3 .3 .3]);
    buildingOrigin2 = [((-7.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim2 = [1.8, 1.8, 12];
    plotcube(dim2, buildingOrigin2, 1, [.3 .3 .3]);
    buildingOrigin13 = [((-7.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim = [1.8, 1.8, 3];
    plotcube(dim, buildingOrigin13, 1, [.3 .3 .3]);
    buildingOrigin14 = [((-7.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim6 = [1.8, 1.8, 4];
    plotcube(dim6, buildingOrigin14, 1, [.3 .3 .3]);
    buildingOrigin3 = [((2.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim3 = [1.8, 1.8, 5];
    plotcube(dim3, buildingOrigin3, 1, [.3 .3 .3]);
    buildingOrigin4 = [((2.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim4 = [1.8, 1.8, 10];
    plotcube(dim4, buildingOrigin4, 1, [.3 .3 .3]);
    buildingOrigin5 = [((12.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim5 = [1.8, 1.8, 6];
    plotcube(dim5, buildingOrigin5, 1, [.3 .3 .3]);
    buildingOrigin6 = [((12.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim6 = [1.8, 1.8, 4];
    plotcube(dim6, buildingOrigin6, 1, [.3 .3 .3]);
    buildingOrigin7 = [((2.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim7 = [1.8, 1.8, 3];
    plotcube(dim7, buildingOrigin7, 1, [.3 .3 .3]);
    buildingOrigin8 = [((2.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim8 = [1.8, 1.8, 3];
    plotcube(dim8, buildingOrigin8, 1, [.3 .3 .3]);
    buildingOrigin15 = [((12.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim7 = [1.8, 1.8, 5];
    plotcube(dim7, buildingOrigin15, 1, [.3 .3 .3]);
    buildingOrigin16 = [((12.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim8 = [1.8, 1.8, 4];
    plotcube(dim8, buildingOrigin16, 1, [.3 .3 .3]);
    buildingFloor = [((-20.5 - minXval)/contourGranularity+.5), ((-20.5 - minYval)/contourGranularity+.5), ((-20- minZval)/contourGranularity+.5)];
    dim5 = [18, 18, .5];
    plotcube(dim5, buildingFloor, 1, [.3 1 .3]);
    
    if (size(ObstacleData,1) > 90)
        for k = 1:90
            obs_x = (ObstacleData(k,1)-minXval)/contourGranularity+.5;
            obs_y = (ObstacleData(k,2)-minYval)/contourGranularity+.5;
            obs_z = (ObstacleData(k,3)-minZval)/contourGranularity+.5;
            if (ObstacleData(k,4) > 0)
                [temp1, temp2, temp3] = sphere;
                temp1 = temp1 * (ObstacleData(k,4)/contourGranularity)*.55;
                temp2 = temp2 * (ObstacleData(k,4)/contourGranularity)*.55;
                temp3 = temp3 * (ObstacleData(k,4)/contourGranularity)*.55;
                if(ObstacleData(k,4) > 2.8)
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*1);
                else
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*.85);
                end
                %plot3(obs_x, obs_y, obs_z, 'o', 'LineWidth', 1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.7, 0.7, 0.7],'MarkerSize', 20*ObstacleData(k,4))
            end
        end
        k = size(ObstacleData,1);
        obs_x = (ObstacleData(k,1)-minXval)/contourGranularity+.5;
        obs_y = (ObstacleData(k,2)-minYval)/contourGranularity+.5;
        obs_z = (ObstacleData(k,3)-minZval)/contourGranularity+.5;
        if (ObstacleData(k,4) > 0)
            [temp1, temp2, temp3] = sphere;
            temp1 = temp1 * ObstacleData(k,4)/contourGranularity;
            temp2 = temp2 * ObstacleData(k,4)/contourGranularity;
            temp3 = temp3 * ObstacleData(k,4)/contourGranularity;
            %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*.85);
            %plot3(obs_x, obs_y, obs_z, 'o', 'LineWidth', 1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.7, 0.7, 0.7],'MarkerSize', 20*ObstacleData(k,4))
        end
    else
        for k = 1:size(ObstacleData,1)
            obs_x = (ObstacleData(k,1)-minXval)/contourGranularity+.45;
            obs_y = (ObstacleData(k,2)-minYval)/contourGranularity+.45;
            obs_z = (ObstacleData(k,3)-minZval)/contourGranularity+.45;
            if (ObstacleData(k,4) > 0)
                [temp1, temp2, temp3] = sphere;
                %if((ObstacleData(k,4) == 1.3) || (ObstacleData(k,4) == .78))
                %    ObstacleData(k,4) = ObstacleData(k,4)/.65;
                %end
                ObstacleData(k,4) = ObstacleData(k,4)/.75;
                temp1 = temp1 * ObstacleData(k,4)/contourGranularity*.8;
                temp2 = temp2 * ObstacleData(k,4)/contourGranularity*.8;
                temp3 = temp3 * ObstacleData(k,4)/contourGranularity*.8;
                if(ObstacleData(k,4) > 2.8)
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*1);
                else
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*1);
                end
                %plot3(obs_x, obs_y, obs_z, 'o', 'LineWidth', 1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.7, 0.7, 0.7],'MarkerSize',20*ObstacleData(k,4))
            end
        end
    end
%    plot(raw_x(2), raw_y(2), 'sw', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',8)
    plot3(c_OBSNodes_x, c_OBSNodes_y, c_OBSNodes_z,'.k')
    plot3(c_OBSNodesN_x, c_OBSNodesN_y, c_OBSNodesN_z, 'ok')
    %plotVehicleTheta(fig, [c_move_x(end) c_move_y(end)], c_move_theta(end), .5, 'k', sensor_radius, '--w')
    
    if length(c_path_x) < 2
        c_path_x = [c_path_x c_path_x];
        c_path_y = [c_path_y c_path_y];
        c_path_z = [c_path_z c_path_z];
    end
    
    if length(c_path_x2) < 2
        c_path_x = [c_path_x2 c_path_x2];
        c_path_y = [c_path_y2 c_path_y2];
        c_path_z = [c_path_z2 c_path_z2];
    end
    
    if length(c_path_x3) < 2
        c_path_x = [c_path_x3 c_path_x3];
        c_path_y = [c_path_y3 c_path_y3];
        c_path_z = [c_path_z3 c_path_z3];
    end
    
    if length(c_path_x4) < 2
        c_path_x = [c_path_x4 c_path_x4];
        c_path_y = [c_path_y4 c_path_y4];
        c_path_z = [c_path_z4 c_path_z4];
    end
    
    %plotVehicle(fig, [c_path_x(1) c_path_y(1)], [c_path_x(2) c_path_y(2)], .5, 'k', sensor_radius, '--w')

    
    
    hold off
    %axis([1, size(Z,1), 1, size(Z,2)-1])
    set(gca,'XTick',[])
    %set(gca,'XTickLabel', theXticks)
    set(gca,'YTick',[])
    %set(gca,'YTickLabel', theYticks)
    set(gca,'ZTick',[])
    %set(gca,'ZTickLabel', theZticks)
    %set(gca,'FontSize',20)
    %colorbar
    view(0,0)
    title('Side 2 View')
    
    %axis tight
    %set(gca,'nextplot','replacechildren');
    %set(gcf,'Renderer','zbuffer');
    subplot(2,2,3);
    
    %clf
    %colormap(tempcolormap)
    %contourf(Zmin,countorLevels, 'EdgeColor', 'none')
    hold on
    %view(135,90)
    %plot3(c_x, c_y, c_z, 'Color',[0.3,0.3,0.3],'LineWidth',.5)
    plot3(c_node_x, c_node_y, c_node_z, '.', 'Color',[0.3,0.3,0.3],'MarkerSize',.5)
    plot3(c_move_x, c_move_y, c_move_z, 'k', 'LineWidth',3)
    plot3(c_move_x, c_move_y, c_move_z, 'r', 'LineWidth',3)
    plot3(c_move_x2, c_move_y2, c_move_z2, 'k', 'LineWidth',3)
    plot3(c_move_x2, c_move_y2, c_move_z2, 'b', 'LineWidth',3)
    plot3(c_move_x3, c_move_y3, c_move_z3, 'k', 'LineWidth',3)
    %plot3(c_move_x3, c_move_y3, c_move_z3, 'b', 'LineWidth',3)
    plot3(c_move_x4, c_move_y4, c_move_z4, 'k', 'LineWidth',3)
    %plot3(c_move_x4, c_move_y4, c_move_z4, 'b', 'LineWidth',3)
    %plot(c_path_x, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x, c_path_y, 'w', 'LineWidth',1)
    %plot(c_path_x2, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x2, c_path_y, 'w', 'LineWidth',1)
    %plot(c_path_x3, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x3, c_path_y, 'w', 'LineWidth',1)
    %plot(c_path_x4, c_path_y, 'k', 'LineWidth',3)
    %plot(c_path_x4, c_path_y, 'w', 'LineWidth',1)
    %contour(Zmin,countorLevels(1:(maxPlotcontourValindPlusOne-1)), 'k')
    %for p = 1:length(polyStartInds)
    %   poly_inds = polyStartInds(p):polyEndEnds(p);
    %   patch(c_obs_x(poly_inds), c_obs_y(poly_inds), 'k') 
    %end
    %plot3(c_obs_x,c_obs_y, c_obs_z, 'w', 'LineWidth',1)
    if start_move_at_ctr > file_ctr
      plot3(c_node_x(1), c_node_y(1), c_node_z(1), 'sw', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',8)
    else
      plot3(c_path_x(end), c_path_y(end), c_node_z(end), 'sw', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',8)
    end
    plot3(c_move_x(end), c_move_y(end), c_move_z(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    if (~isnan(move_x4))
    plot3(c_move_x2(end), c_move_y2(end), c_move_z2(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    plot3(c_move_x3(end), c_move_y3(end), c_move_z3(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    plot3(c_move_x4(end), c_move_y4(end), c_move_z4(end), 'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[1.000000000, 0.819607843, 0.137254902],'MarkerSize',12)
    end
    buildingOrigin9 = [((-17.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim = [1.8, 1.8, 3];
    plotcube(dim, buildingOrigin9, 1, [.3 .3 .3]);
    buildingOrigin10 = [((-17.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim2 = [1.8, 1.8, 7];
    plotcube(dim2, buildingOrigin10, 1, [.3 .3 .3]);
    buildingOrigin11 = [((-17.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim3 = [1.8, 1.8, 4];
    plotcube(dim3, buildingOrigin11, 1, [.3 .3 .3]);
    buildingOrigin12 = [((-17.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim4 = [1.8, 1.8, 3];
    plotcube(dim4, buildingOrigin12, 1, [.3 .3 .3]);
    buildingOrigin = [((-7.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim = [1.8, 1.8, 8];
    plotcube(dim, buildingOrigin, 1, [.3 .3 .3]);
    buildingOrigin2 = [((-7.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim2 = [1.8, 1.8, 12];
    plotcube(dim2, buildingOrigin2, 1, [.3 .3 .3]);
    buildingOrigin13 = [((-7.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim = [1.8, 1.8, 3];
    plotcube(dim, buildingOrigin13, 1, [.3 .3 .3]);
    buildingOrigin14 = [((-7.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim6 = [1.8, 1.8, 4];
    plotcube(dim6, buildingOrigin14, 1, [.3 .3 .3]);
    buildingOrigin3 = [((2.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim3 = [1.8, 1.8, 5];
    plotcube(dim3, buildingOrigin3, 1, [.3 .3 .3]);
    buildingOrigin4 = [((2.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim4 = [1.8, 1.8, 10];
    plotcube(dim4, buildingOrigin4, 1, [.3 .3 .3]);
    buildingOrigin5 = [((12.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim5 = [1.8, 1.8, 6];
    plotcube(dim5, buildingOrigin5, 1, [.3 .3 .3]);
    buildingOrigin6 = [((12.5 - minXval)/contourGranularity+.5), ((2.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim6 = [1.8, 1.8, 4];
    plotcube(dim6, buildingOrigin6, 1, [.3 .3 .3]);
    buildingOrigin7 = [((2.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim7 = [1.8, 1.8, 3];
    plotcube(dim7, buildingOrigin7, 1, [.3 .3 .3]);
    buildingOrigin8 = [((2.5 - minXval)/contourGranularity+.5), ((12.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim8 = [1.8, 1.8, 3];
    plotcube(dim8, buildingOrigin8, 1, [.3 .3 .3]);
    buildingOrigin15 = [((12.5 - minXval)/contourGranularity+.5), ((-17.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim7 = [1.8, 1.8, 5];
    plotcube(dim7, buildingOrigin15, 1, [.3 .3 .3]);
    buildingOrigin16 = [((12.5 - minXval)/contourGranularity+.5), ((-7.5 - minYval)/contourGranularity+.5), ((-18- minZval)/contourGranularity+.5)];
    dim8 = [1.8, 1.8, 4];
    plotcube(dim8, buildingOrigin16, 1, [.3 .3 .3]);
    buildingFloor = [((-20.5 - minXval)/contourGranularity+.5), ((-20.5 - minYval)/contourGranularity+.5), ((-20- minZval)/contourGranularity+.5)];
    dim5 = [18, 18, .5];
    plotcube(dim5, buildingFloor, 1, [.3 1 .3]);
    
    if (size(ObstacleData,1) > 90)
        for k = 1:90
            obs_x = (ObstacleData(k,1)-minXval)/contourGranularity+.5;
            obs_y = (ObstacleData(k,2)-minYval)/contourGranularity+.5;
            obs_z = (ObstacleData(k,3)-minZval)/contourGranularity+.5;
            if (ObstacleData(k,4) > 0)
                [temp1, temp2, temp3] = sphere;
                temp1 = temp1 * (ObstacleData(k,4)/contourGranularity)*.55;
                temp2 = temp2 * (ObstacleData(k,4)/contourGranularity)*.55;
                temp3 = temp3 * (ObstacleData(k,4)/contourGranularity)*.55;
                if(ObstacleData(k,4) > 2.8)
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*1);
                else
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*.85);
                end
                %plot3(obs_x, obs_y, obs_z, 'o', 'LineWidth', 1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.7, 0.7, 0.7],'MarkerSize', 20*ObstacleData(k,4))
            end
        end
        k = size(ObstacleData,1);
        obs_x = (ObstacleData(k,1)-minXval)/contourGranularity+.5;
        obs_y = (ObstacleData(k,2)-minYval)/contourGranularity+.5;
        obs_z = (ObstacleData(k,3)-minZval)/contourGranularity+.5;
        if (ObstacleData(k,4) > 0)
            [temp1, temp2, temp3] = sphere;
            temp1 = temp1 * ObstacleData(k,4)/contourGranularity;
            temp2 = temp2 * ObstacleData(k,4)/contourGranularity;
            temp3 = temp3 * ObstacleData(k,4)/contourGranularity;
            %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*.85);
            %plot3(obs_x, obs_y, obs_z, 'o', 'LineWidth', 1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.7, 0.7, 0.7],'MarkerSize', 20*ObstacleData(k,4))
        end
    else
        for k = 1:size(ObstacleData,1)
            obs_x = (ObstacleData(k,1)-minXval)/contourGranularity+.45;
            obs_y = (ObstacleData(k,2)-minYval)/contourGranularity+.45;
            obs_z = (ObstacleData(k,3)-minZval)/contourGranularity+.45;
            if (ObstacleData(k,4) > 0)
                [temp1, temp2, temp3] = sphere;
                %if((ObstacleData(k,4) == 1.3) || (ObstacleData(k,4) == .78))
                %    ObstacleData(k,4) = ObstacleData(k,4)/.65;
                %end
                ObstacleData(k,4) = ObstacleData(k,4)/.75;
                temp1 = temp1 * ObstacleData(k,4)/contourGranularity*.8;
                temp2 = temp2 * ObstacleData(k,4)/contourGranularity*.8;
                temp3 = temp3 * ObstacleData(k,4)/contourGranularity*.8;
                if(ObstacleData(k,4) > 2.8)
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*1);
                else
                    %surf((temp1 + obs_x), (temp2 + obs_y), (temp3 + obs_z), ones(size(temp1))*1);
                end
                %plot3(obs_x, obs_y, obs_z, 'o', 'LineWidth', 1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.7, 0.7, 0.7],'MarkerSize',20*ObstacleData(k,4))
            end
        end
    end
%    plot(raw_x(2), raw_y(2), 'sw', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',8)
    plot3(c_OBSNodes_x, c_OBSNodes_y, c_OBSNodes_z,'.k')
    plot3(c_OBSNodesN_x, c_OBSNodesN_y, c_OBSNodesN_z, 'ok')
    %plotVehicleTheta(fig, [c_move_x(end) c_move_y(end)], c_move_theta(end), .5, 'k', sensor_radius, '--w')
    
    if length(c_path_x) < 2
        c_path_x = [c_path_x c_path_x];
        c_path_y = [c_path_y c_path_y];
        c_path_z = [c_path_z c_path_z];
    end
    
    if length(c_path_x2) < 2
        c_path_x = [c_path_x2 c_path_x2];
        c_path_y = [c_path_y2 c_path_y2];
        c_path_z = [c_path_z2 c_path_z2];
    end
    
    if length(c_path_x3) < 2
        c_path_x = [c_path_x3 c_path_x3];
        c_path_y = [c_path_y3 c_path_y3];
        c_path_z = [c_path_z3 c_path_z3];
    end
    
    if length(c_path_x4) < 2
        c_path_x = [c_path_x4 c_path_x4];
        c_path_y = [c_path_y4 c_path_y4];
        c_path_z = [c_path_z4 c_path_z4];
    end
    
    %plotVehicle(fig, [c_path_x(1) c_path_y(1)], [c_path_x(2) c_path_y(2)], .5, 'k', sensor_radius, '--w')

    
    
    hold off
    %axis([1, size(Z,1), 1, size(Z,2)-1])
    set(gca,'XTick',[])
    %set(gca,'XTickLabel', theXticks)
    set(gca,'YTick',[])
    %set(gca,'YTickLabel', theYticks)
    set(gca,'ZTick',[])
    %set(gca,'ZTickLabel', theZticks)
    %set(gca,'FontSize',20)
    %colorbar
    view(90,90)
    title('Top View')
    
    %axis tight
    %set(gca,'nextplot','replacechildren');
    %set(gcf,'Renderer','zbuffer');

    F = getframe(fig);
    writeVideo(writerObj,F);
    
    file_ctr = file_ctr + 1; 
    %if file_ctr >= 20
    %    pause()
    %end
end
% close video file
close(writerObj);