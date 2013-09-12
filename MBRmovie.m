%% function to create a movie simulation from the time and state of MBR
%% from running MBR_gillespie_func.m

function [] = MBRmovie(timeVec, state)
    
    % video options
    timestep = 1/30;% 30fps
    
    % unpack state
    x = state.posn(:,1);
    y = state.posn(:,2);
    th = state.posn(:,3);
    thunwrap = rad2deg(unwrap(deg2rad(th)));
    th = thunwrap;
    numcell = size([state.cellposn],1);
       
    xLimits = [min(x) max(x)];
    xLimits = round(xLimits/10)*10;
    yLimits = [min(y) max(y)];
    yLimits = round(yLimits/10)*10;
    
    F = state.F;
    
    % create object for movie
    mbrViewObj = VideoWriter('mbrGillespieSim1.avi');
    open(mbrViewObj);
    
    % create movie figure
    mov1 = figure;
    axis([-50 100 -100 50])

    %buffer = 50;
    %axis([xLimits(1) - buffer, xLimits(2) + buffer,yLimits(1) - buffer, yLimits(2) + buffer])
    axis equal
    axis manual
    % plot first frame
    % plot MBR
    diagl = sqrt(20^2 + 20^2);
    diagl = sqrt(30^2 + 30^2);
    cornerTh = [45:90:360]'+th(1);
   
    
    
    mbrCorners = diagl*[cosd(cornerTh) sind(cornerTh)] + repmat([x(1),y(1)],[4,1]);
    mbrCorners = [mbrCorners;mbrCorners(1,:)];
    
    % convert from um to m:
    %mbrCorners = mbrCorners * 10^-6;
    celllength = 3;
    
    mbrsq = line(mbrCorners(:,1),mbrCorners(:,2));
    set(mbrsq,'color',[0 0 0]);
        
    rotationMatrix = [cosd(th(1)) -sind(th(1));sind(th(1)) cosd(th(1))];
    
    bacHead = rotationMatrix * state.cellposn(:,1:2)';
    bacHead = bacHead';
    cellAngle = state.cellAngle + repmat(th,[1,numcell]); % in world frame
    dbac = celllength*[cosd(cellAngle(1,:))' sind(cellAngle(1,:))'];
    bacTail = bacHead + dbac;
        
    % plot cells
    for j = 1:size(state.cellposn,1)
        bacX = [bacHead(j,1);bacTail(j,1)];
        bacY = [bacHead(j,2);bacTail(j,2)];
        hold on
        cellplot(j) = plot(bacX,bacY,'-b','LineWidth',2);
    end
    
    %htime =  text(xLimits(2),yLimits(2),sprintf('Time = %.2f s',timeVec(1)),'horizontalAlignment','center');
    htime =  text(40,40,sprintf('Time = %.2f s',timeVec(1)),'horizontalAlignment','center');
    
    hold on 
    % plot trajectory
    plot(x,y,'-g')
    xlabel('X-coordinate')
    ylabel('Y-coordinate')
    title('Gillespie Simulation of MBR')
    
    % save frame to video
    frame = getframe(gcf);
    open(mbrViewObj)
    writeVideo(mbrViewObj,frame);
    % close(mbrViewObj)
    
    %for i = 2:length(timeVec)
    i = 1;
    timeCurrent = 0;
    plot(x(1),y(1),'o','MarkerSize',5,'LineWidth',3,'Color','g')
    plot(x(end),y(end),'o','MarkerSize',5,'LineWidth',3,'Color','r')
    while timeCurrent + timestep < timeVec(end)
        
                
        
        % update timeCurrent
        timeCurrent = timeCurrent + timestep;
        
        % update the index i
        i = find(timeVec<timeCurrent,1,'last');

        % update time vector
        set(htime,'string',sprintf('Time = %.2f s',timeCurrent))
        
        % compute angle of MBR
        thRate = (th(i+1) - th(i))/(timeVec(i+1)-timeVec(i));
        thNow = th(i) + thRate*(timeCurrent-timeVec(i));
        
        cornerTh = [45:90:360]'+ thNow;       
        % update MBR position
        vel = [diff([x(i),x(i+1)]);diff([y(i),y(i+1)])]/(timeVec(i+1) - timeVec(i));
        xNow = x(i) + vel(1)*(timeCurrent - timeVec(i));
        yNow = y(i) + vel(2)*(timeCurrent - timeVec(i));
        
        mbrCorners = diagl*[cosd(cornerTh) sind(cornerTh)] + repmat([xNow,yNow],[4,1]);
        mbrCorners = [mbrCorners;mbrCorners(1,:)];
        
        set(mbrsq,'XData',mbrCorners(:,1),'YData',mbrCorners(:,2));
                
        rotationMatrix = [cosd(thNow) -sind(thNow);sind(thNow) cosd(thNow)];
        
        bacHead = rotationMatrix * state.cellposn(:,1:2)';
        bacHead = bacHead' + repmat([xNow,yNow],[numcell,1]);
        
        dbac = celllength*[cosd(cellAngle(i,:)+thRate*(timeCurrent-timeVec(i)))' sind(cellAngle(i,:)+thRate*(timeCurrent-timeVec(i)))'];

        
        bacTail = bacHead + dbac;
        
        % plot cells
        for j = 1:size(state.cellposn,1)
            if F(i,j) == 0
                
                set(cellplot(j),'XData',bacHead(j,1),'YData',bacHead(j,2))
                set(cellplot(j),'Color',[1 0 0],'MarkerSize',3,'Marker','o')
                
            else
                bacX = [bacHead(j,1);bacTail(j,1)];
                bacY = [bacHead(j,2);bacTail(j,2)];
                set(cellplot(j),'XData',bacX,'YData',bacY)
                set(cellplot(j),'Color','b','MarkerSize',3,'Marker','none')
            end
        end
        axis([-50 100 -100 50])
        %axis equal

        drawnow
        
        frame = getframe(gcf);
        open(mbrViewObj)
        writeVideo(mbrViewObj,frame);
%         close(mbrViewObj)
        
        tfactor = 0.1;
        if i == length(timeVec)
            disp('End of Simulation')
            close(mbrViewObj)
        else
            pause((timeVec(i+1) - timeVec(i))*tfactor)
        end
        
    end