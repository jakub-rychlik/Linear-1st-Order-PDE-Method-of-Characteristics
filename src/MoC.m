function MoC(U,vtk,tol,plot_index)

    % ---INITIALIZATION---
    % Interval I is parametrizing curve which creates entry boundary. 
    global b1 b2 c f g set set_bdd Act_en_bdd S;
    format long
    eps = 250;
    I = U(1):0.15:U(2);  % partition of I
    I = [I,U(2)];

    % Discretizing entry boundary and saving some indexes from the creation
    % which will be needed for creating even more points afterwards in 'Adding characteristics' step.
    disp('Discretizing entry boundary...');
    [disc_en_bdd, bdd_steps] = Disc_entry_bdd(I,tol,eps);


    % ---CREATING CHARACTERISTICS---
    disp('Creating characteristic curves...');
    h=0.05;  % iteration step for Runge-Kutta method
    char_approx = cell(1, length(disc_en_bdd(1,:)));
    parametrization_point_dist = cell(1, length(disc_en_bdd(1,:)));

    for K=1:length(disc_en_bdd(1,:))
        [char_approx{K}, parametrization_point_dist{K}] = Create_characteristic(disc_en_bdd(:,K),h,tol);
    end
    

    % ---ADDING CHARACTERISTICS---
    % Adding characteristics, so they are not far away from each other.
    % It will also create new pointss on disc_en_bdd 
    disp('Adding more curves...');
    [char_approx, disc_en_bdd,parametrization_point_dist] = Add_characteristics(char_approx,disc_en_bdd,tol,bdd_steps,I,eps,parametrization_point_dist);


    % ---REMOVING CHARACTERISTICS---
    % Removing unnecessary characteristics, because it is easier than adding specific quantity of new ones.
    disp('Filtration of all points...');
    check = true;
    while check
        k = 1;
        K = length(char_approx);
        check = false;

        while k+2 < K
            exit_dist = round(norm(char_approx{k}(:,end) - char_approx{k+2}(:,end)),8);
            entry_dist = round(norm(char_approx{k}(:,1) - char_approx{k+2}(:,1)),8);

            if ~(exit_dist > tol) && ~(exit_dist < entry_dist)
                char_approx = {char_approx{1:k},char_approx{k+2:end}};
                disc_en_bdd = [disc_en_bdd(:,1:k),disc_en_bdd(:,k+2:end)];
                parametrization_point_dist = {parametrization_point_dist{1:k},parametrization_point_dist{k+2:end}};
                check = true;
                K = K-1;
            end

            k = k+1;
            if k+2 > K
                break
            end
        end
    end


    % ---REMOVING POINTS---
    % Removing unnecessary points on characteristics.
    Total_points = 0;
    for K=1:length(char_approx)
        check = true;
        line = char_approx{K};
        point_dist = parametrization_point_dist{K};

        while check
            l = 1;
            L = length(line(1,:));
            check = false;

            while l+2 < L
                d = round(norm(line(:,l) - line(:,l+2)),8)+10^(-8);
                if ~(d > tol)
                    line = [line(:,1:l),line(:,l+2:end)];
                    point_dist(l) = point_dist(l) + point_dist(l+1);
                    point_dist = [point_dist(1:l),point_dist(l+2:end)];
                    check = true;
                    L = L-1;
                end

                l = l+1;
                if l+2 > L
                    break
                end
            end
        end

        parametrization_point_dist{K} = point_dist;
        char_approx{K} = line;
        Total_points = Total_points + length(char_approx{K}(1,:));
    end


    % ---DISCRETIZING IVC---
    disc_ivc = [zeros(1,length(disc_en_bdd(1,:)))];
    for k=1:length(disc_en_bdd(1,:))
        disc_ivc(k) = Ivc(disc_en_bdd(1,k),disc_en_bdd(2,k));
    end
    

    % ---CREATING TRIANGLE NET---
    % Creating triangle net by assigning vertices of triangles to every point of all characteristic (using structures).
    disp('Creating triangle net...');
    [Chars_w_vertices,vtkPoints,vtkCells] = Create_net(char_approx,Total_points);


    % ---SOLUTION---
    % Finding approximate solution in vertices.
    disp('Approximating solution in vertices of triangles...');
    Chars_w_vertices = Approximate_solution_vertices(Chars_w_vertices,parametrization_point_dist,disc_ivc);
   

    % ---CREATING VTK FILE---
    % If vtk = true, it creates .vtk file for Paraview
    if vtk
        disp('Creating VTK file...');
        createVTK(Chars_w_vertices,vtkPoints,vtkCells);
    end
    

    % ---CREATING BRIEF BOUNDARIES FOR OMEGA---
    xL=Inf;
    xR=-Inf;
    yD=Inf;
    yU=-Inf;
    for line=1:length(Chars_w_vertices)
        for point=1:length(Chars_w_vertices{line})
            if Chars_w_vertices{line}(point).coor(1) < xL
                xL = Chars_w_vertices{line}(point).coor(1);
            end

            if Chars_w_vertices{line}(point).coor(1) > xR
               xR  = Chars_w_vertices{line}(point).coor(1);
            end
            
            if Chars_w_vertices{line}(point).coor(2) < yD
                yD = Chars_w_vertices{line}(point).coor(2);
            end
            
            if Chars_w_vertices{line}(point).coor(2) > yU
                yU = Chars_w_vertices{line}(point).coor(2);
            end
        end
    end
    
    
    % 3D plots 12 lines, which were created by cutting the surface of approximated solution
    if plot_index == 1
        figure;
    end
    i=1;
    disp('Approximating solution inside triangles...');
    for y = yD:((yU-yD)/12):yU   
        Sol_approx = [];
        x_list = [];
        for x = xL:0.05:xR
            u = Approximate_solution_inside(Chars_w_vertices,[x;y]);
            if ~isnan(u)
                Sol_approx = [Sol_approx, u];
                x_list = [x_list, x];
            end
        end
        y_list = [zeros(1,length(x_list))]+y;
        if i > 1
            plot3(x_list,y_list,Sol_approx,'Color','black');
            hold on;
        end
        i = i+1;
    end
    hold on;

    % plots characteristics on xy-plane
    for l=1:(length(char_approx)) 
        plot(char_approx{l}(1,:),char_approx{l}(2,:));
        hold on;
    end

    % plots boundary of Omega on xy-plane
    for i = 1:length(S)
        border = S{i};
        x_list = [];
        y_list = [];
        interval = border.int;
        for step = interval(1):tol/3:interval(2)+0.1
            x_list = [x_list, border.xfunc(step)];
            y_list = [y_list, border.yfunc(step)];
        end
        plot(x_list,y_list,'Color','black','LineWidth',1);
        hold on;
    end

end