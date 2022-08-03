function animate(sim_name,... % e.g. 'R_equals_1L_sphere6_f_220'
    start_fraction,... % number between 0 and 1 (inclusive)
    end_fraction,... % number between start_fraction and 1 (inclusive)
    show_figure,... % true or false
    view_vec,... % 2, 3 or a vector of the form [az el]
    show_time,... % true or false
    num_shown_frames,... % an integer >= 1
    body_fixed_view,... % true or false
    type,... % 'filament', 'tip_disp_vector_field' or 'prescribed'
    show_rotation_axis,... % true or false
    show_body_fixed_frame... % true or false
    )

% N.B. Only works for NBOD=1 because I don't currently run multi-body simulations.

% N.B. If the simulation was run in the normal git directory on a machine known as "remote" in your ssh
% config file, all strictly necessary data files can be downloaded by running the following bash command:
%
% for type in par reference state; do scp remote:~/phd_codes/tethered_filaments/sim_name*$type* .; done
%
% Some options, such as show_rotation_axis = true, may require additional files.

% -------------------------------------------------------------------- %
% Open the data files
% -------------------------------------------------------------------- %

fid_body = fopen([sim_name, '_body_states.dat']);
fid_seg = fopen([sim_name, '_seg_states.dat']);

parameters = load([sim_name, '.par']);
NFIL = parameters(1);
NSEG = parameters(2);
RSEG = parameters(6);

L = 2.2*RSEG*NSEG;

num_frames = 1 + floor((parameters(11)-1)/parameters(13));

if exist([sim_name, '_blob_references.dat'], 'file') == 2

    is_planar_sim = false;

    if show_rotation_axis

        W = load([sim_name '_body_vels.dat']);
        W = W(:,5:7);

    end

else

    is_planar_sim = true;

    show_rotation_axis = false;

end

fil_references = reshape(load([sim_name, '_fil_references.dat']), 3, NFIL);

start_fraction = min(1, start_fraction);
start_fraction = max(0, start_fraction);
end_fraction = min(1, end_fraction);
end_fraction = max(end_fraction, start_fraction);

try

    % We may have not bothered to download the backup file, but if we
    % happen to have it then we can use it to get a tight upper bound on
    % how far the simulation reached.
    temp = load([sim_name '.backup']);
    goal_steps = parameters(11);
    upper_bound_on_step_reached = min([temp(1) + parameters(13), goal_steps]);

    if upper_bound_on_step_reached <= start_fraction*goal_steps

        % We expect to reach the end of the data before we actually start
        % plotting.
        fprintf('\n Backup file suggests the simulation ran for at most %g%% of the maximum steps.\n The choice start_fraction = %g is not expected to produce results. \n\n', 100*upper_bound_on_step_reached/parameters(11), start_fraction);
        return;
        
    end
    
    % If we're asking our video to finish after the end of the file, reduce
    % the requested video length.
    end_fraction = min([end_fraction, upper_bound_on_step_reached/parameters(11)]);

end

% -------------------------------------------------------------------- %
% Prepare the video file
% -------------------------------------------------------------------- %

if body_fixed_view

    video_name = 'body_frame';

else

    video_name = 'lab_frame';

end

if strcmp(type, 'prescribed')
    
    video_name = [sim_name '_' video_name '_' type '_animation_from_' num2str(start_fraction*parameters(11)/parameters(14)) 'T_to_' num2str(end_fraction*parameters(11)/parameters(14)) 'T'];
    
else
    
     video_name = [sim_name '_' video_name '_' type '_animation_from_' num2str(start_fraction*parameters(11)*parameters(12)/36.3107) 'T_to_' num2str(end_fraction*parameters(11)*parameters(12)/36.3107) 'T'];
    
end

if show_rotation_axis

    video_name = [video_name '_with_rotation_axis'];

end

if show_body_fixed_frame

    video_name = [video_name '_with_body_fixed_frame'];

end

try

    % .mp4 is better than .avi in terms of space, but some remote machines
    % refuse to do .mp4
    video = VideoWriter(video_name, 'MPEG-4');

catch

    video = VideoWriter(video_name);

end

video.Quality = 100;
video.FrameRate = 30;
open(video);

% -------------------------------------------------------------------- %
% Initialise arrays etc.
% -------------------------------------------------------------------- %

read_line = @(fid, line_length) textscan(fid, repmat('%f', [1 line_length]), 1, 'CommentStyle', '%', 'Delimiter', ' ');

seg_pos = NaN(3, NSEG*NFIL + NFIL - 1, num_shown_frames);

if strcmp(type, 'tip_disp_vector_field')

    tip_disp = NaN(3, NFIL, num_shown_frames);
    base_pos = NaN(3, NFIL, num_shown_frames);

end

if is_planar_sim

    xmin = min(fil_references(1,:)) - L;
    xmax = max(fil_references(1,:)) + L;
    ymin = min(fil_references(2,:)) - L;
    ymax = max(fil_references(2,:)) + L;
    zmin = 0;
    zmax = 1.1*L;

else

    % TO-DO: How could I generalise this to non-spherical bodies?

    [s1, s2, s3] = sphere(100);

    sphere_radius = norm(fil_references(:,1));

    xmin = -sphere_radius - L;
    xmax = sphere_radius + L;
    ymin = -sphere_radius - L;
    ymax = sphere_radius + L;
    zmin = -sphere_radius - L;
    zmax = sphere_radius + L;

    sphere_radius = sphere_radius - 1.1*RSEG;
    s1 = sphere_radius*s1;
    s2 = sphere_radius*s2;
    s3 = sphere_radius*s3;

end

% -------------------------------------------------------------------- %
% Ensure the figure borders are minimal
% -------------------------------------------------------------------- %
% Note: the quality of this 'border removal' is dependent on the choice of
% view(...) and consequently how well the axis box fits the enclosed image.

% TO-DO: Is there a less 'hacky' way of doing this?

if show_figure

    h = figure;

else

    h = figure('visible', 'off');

end

h.CurrentAxes = axes;
ax = h.CurrentAxes;
axis(ax, 'equal');
axis(ax, [xmin xmax ymin ymax zmin zmax]);
view(ax, view_vec);
box(ax, 'on');
xticks(ax, []);
yticks(ax, []);
zticks(ax, []);
M = print(h, '-RGBImage', '-r1000');
M1 = double(M(:,:,1));
[row, col] = find(M1 - 255); % uint8 runs from 0 to 255
c = zeros(4,1);
r = zeros(4,1);
[c(1), id] = min(col); r(1) = row(id);
[c(2), id] = max(col); r(2) = row(id);
[r(3), id] = min(row); c(3) = col(id);
[r(4), id] = max(row); c(4) = col(id);
ax.Units = 'pixels';
h.Units = 'pixels';
height = max(r) - min(r);
width = max(c) - min(c);
if (height > width) % keep shortest dimension same as standard figure window
    ax.Position(4) = ax.Position(3)*height/width;
    h.Position(4) = h.Position(3)*height/width;
else
    ax.Position(3) = ax.Position(4)*width/height;
    h.Position(3) = h.Position(4)*width/height;
end
h.Position(1:2) = [100 100];
set(ax, 'Units', 'normalized', 'Position', [0 0 1 1]);
set(h, 'color', 'w'); % white background for better contrast

% -------------------------------------------------------------------- %
% Attempt to make the filament width accurate
% -------------------------------------------------------------------- %

[seg_s1, seg_s2, seg_s3] = sphere(100);
hold on;
seg_sphere_handle = surf(ax, RSEG*seg_s1 + 0.5*(xmin + xmax), RSEG*seg_s2 + 0.5*(ymin + ymax), RSEG*seg_s3 + 0.5*(zmin + zmax), 'EdgeColor', 'none', 'FaceColor', [0 0 0]);
axis off;
M = print(h, '-RGBImage', '-r1000');
M1 = double(M(:,:,1));
[row, col] = find(M1 - 255); % uint8 runs from 0 to 255
c = zeros(4,1);
r = zeros(4,1);
[c(1), id] = min(col); r(1) = row(id);
[c(2), id] = max(col); r(2) = row(id);
[r(3), id] = min(row); c(3) = col(id);
[r(4), id] = max(row); c(4) = col(id);
seg_height = max(r) - min(r);
seg_width = max(c) - min(c);
seg_size = 0.5*(seg_width + seg_height);
set(ax, 'Units', 'points');
line_width = 0.5*seg_size*(ax.Position(3)/size(M1,2) + ax.Position(4)/size(M1,1)); % take the average of the four possible widths we could calculate.
set(ax, 'Units', 'normalized', 'Position', [0 0 1 1]);
delete(seg_sphere_handle);

% -------------------------------------------------------------------- %
% Create the animation
% -------------------------------------------------------------------- %

hold on;

if ~is_planar_sim
    surf(ax, s1, s2, s3, 'Edgecolor', 'none', 'FaceColor', [0.8 0.8 0.8]);
end

if show_time
     h_string = text(ax, 0, 1, '$t/T = 0$', 'Interpreter', 'latex', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 24);
end

for n = 1:num_frames

    try % fails if we reach the end of the file early; i.e. the simulation didn't run to the final time-step.

        D_body = read_line(fid_body, 8);
        
        if strcmp(type, 'prescribed')
            
            D_seg = read_line(fid_seg, 1 + NFIL*NSEG*3);
            
        else
            
            D_seg = read_line(fid_seg, 1 + NFIL*NSEG*4);
        
        end

        if ((n > start_fraction*num_frames) && (n <= end_fraction*num_frames))

            if strcmp(type, 'prescribed')
                
                curr_time  = D_body{1}/parameters(14);
                
            else
                
                curr_time  = D_body{1}*parameters(12)/36.3107;
                
            end

            D_body = cell2mat(D_body(2:end));

            if is_planar_sim

                R = 1;

            else

                q = D_body(4:7);
                qsq = q.^2;
                R = [1 - 2*(qsq(3) + qsq(4)), 2*(q(2)*q(3) - q(1)*q(4)), 2*(q(2)*q(4) + q(1)*q(3));...
                    2*(q(2)*q(3) + q(1)*q(4)), 1 - 2*(qsq(2) + qsq(4)), 2*(q(3)*q(4) - q(1)*q(2));...
                    2*(q(2)*q(4) - q(1)*q(3)), 2*(q(3)*q(4) + q(1)*q(2)), 1 - 2*(qsq(2) + qsq(3))];

            end
           
            D_seg = cell2mat(D_seg(2:end));

            try

                delete(handles);

                if show_rotation_axis

                    delete(rot_axis_handle);

                end
                
                if show_body_fixed_frame
                    
                    delete(frame_handle);
                    
                end

            end

            for m=1:num_shown_frames-1

                seg_pos(:,:,m) = seg_pos(:,:,m+1);

                if strcmp(type, 'tip_disp_vector_field')

                    tip_disp(:,:,m) = tip_disp(:,:,m+1);
                    base_pos(:,:,m) = base_pos(:,:,m+1);

                end

            end

            for j=NFIL:-1:1
                
                if strcmp(type, 'prescribed')
                    
                    Dfil = D_seg(1 + 3*(j-1)*NSEG: 3*j*NSEG);
                    
                else

                    Dfil = D_seg(1 + 4*(j-1)*NSEG: 4*j*NSEG);
                
                end

                p = j + (j-1)*NSEG;

                if body_fixed_view

                    seg_pos(:,p,end) = fil_references(:,j);

                else

                    seg_pos(:,p,end) = R*fil_references(:,j);

                end

                if strcmp(type, 'tip_disp_vector_field')

                    base_pos(:,j,end) = seg_pos(:,p,end);

                end

                for k=2:NSEG
                    
                    if strcmp(type, 'prescribed')
                        
                        id = 3*(k-1) + 1;
                        p = p + 1;
                        
                        if body_fixed_view
                            
                            seg_pos(:,p,end) = R'*(Dfil(id:id+2)' - D_body(1:3)');
                            
                        else
                            
                            seg_pos(:,p,end) = Dfil(id:id+2)' - D_body(1:3)';
                            
                        end
                        
                    else

                        id = 4*(k-2) + 1;

                        q = Dfil(id:id+3);
                        qsq = q.^2;
                        t = [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];

                        id = id + 4;

                        q = Dfil(id:id+3);
                        qsq = q.^2;
                        t = t + [1 - 2*(qsq(3) + qsq(4)); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))];

                        t = 1.1*RSEG*t;

                        p = p + 1;

                        if body_fixed_view

                            seg_pos(:,p,end) = seg_pos(:,p-1,end) + R'*t;

                        else

                            seg_pos(:,p,end) = seg_pos(:,p-1,end) + t;

                        end
                    
                    end

                end

                if strcmp(type, 'tip_disp_vector_field')

                    tip_disp(:,j,end) = seg_pos(:,p,end) - base_pos(:,j,end);

                    if ~is_planar_sim % project into tangent plane if on a spherical surface

                        tip_disp(:,j,end) = tip_disp(:,j,end) - dot(tip_disp(:,j,end), base_pos(:,j,end))*base_pos(:,j,end)/dot(base_pos(:,j,end), base_pos(:,j,end));

                    end

                end

            end

            for m=num_shown_frames:-1:1

                if (strcmp(type, 'filament') || strcmp(type, 'prescribed'))

                    handles(m) = plot3(ax, seg_pos(1,:,m), seg_pos(2,:,m), seg_pos(3,:,m), '-', 'LineWidth', line_width, 'Color', [1 1 1]*(1 - m/num_shown_frames)); % darker curves for more recent data

                elseif strcmp(type, 'tip_disp_vector_field')

                    if is_planar_sim

                        handles(m) = quiver(ax, base_pos(1,:,m), base_pos(2,:,m), tip_disp(1,:,m), tip_disp(2,:,m), 'Color', [1 1 1]*(1 - m/num_shown_frames));

                    else

                        handles(m) = quiver3(ax, base_pos(1,:,m), base_pos(2,:,m), base_pos(3,:,m), tip_disp(1,:,m), tip_disp(2,:,m), tip_disp(3,:,m), 'Color', [1 1 1]*(1 - m/num_shown_frames));

                    end

                end

            end

            if show_rotation_axis

                Wvec = W(n,:)';
                Wvec = (L + sphere_radius)*Wvec/norm(Wvec);

                if body_fixed_view

                    Wvec = R' * Wvec;

                end

                rot_axis_handle = quiver3(ax, [0 0], [0 0], [0 0], Wvec(1)*[1 -1], Wvec(2)*[1 -1], Wvec(3)*[1 -1], 0, 'r-');

            end
            
            if show_body_fixed_frame

                if body_fixed_view

                    frame_handle = quiver3(ax, [0 0 0], [0 0 0], [0 0 0], (L + sphere_radius)*[1 0 0], (L + sphere_radius)*[0 1 0], (L + sphere_radius)*[0 0 1], 0, 'b-');
                    
                else
                    
                    frame_handle = quiver3(ax, [0 0 0], [0 0 0], [0 0 0], (L + sphere_radius)*R(1,:), (L + sphere_radius)*R(2,:), (L + sphere_radius)*R(3,:), 0, 'b-');

                end

            end

            if show_time
                h_string.String = sprintf('$t/T = %g$', curr_time);
            end

            if show_figure

                frame = getframe(h);

            else

                % Bypass the need to use getframe(), which sometimes fails
                % when 'visible' is set to 'off'.
                print(h, 'temp_frame_for_animation', '-djpeg', '-r300');
                frame = imread('temp_frame_for_animation.jpg');

            end

            writeVideo(video, frame);
            fprintf('Written frame %i.\n', n);
            
        elseif (n > end_fraction*num_frames)
            
            break;

        end

    catch my_error
        
        disp(my_error);

        close(video);

        if ~show_figure
            close(h);
        end

        fclose all;
        return;

    end

end

close(video);

if ~show_figure
    close(h);
end

fclose all;

end
