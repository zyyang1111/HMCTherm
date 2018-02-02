function [ best_topology, best_topology_delay ] = identify_best_topology( topologies, num_MC, core_area, MC_area, single_cycle_WL )

    max_area_per_layer=400; %maximum footprint area
    show_plots=0;

    num_cores=topologies(1,1)*topologies(1,2)*topologies(1,3);
    
    MC_per_core=num_MC/num_cores;
    
    total_area=num_cores*core_area+num_MC*MC_area;
    min_layers=ceil(total_area/max_area_per_layer);

    num_layers=topologies(:,3);
    num_cols=topologies(:,1);
    num_rows=topologies(:,2);
    
    for i=1:size(topologies,1)
        AR(i,1)=topologies(i,2)/topologies(i,1);
        y(i,1)=sqrt((core_area+MC_per_core*MC_area)/AR(i));
        x(i,1)=(core_area+MC_per_core*MC_area)/y(i);
    end
    
    router_dist=max(x,y);
    router_delay=ceil(router_dist/single_cycle_WL);
    
    norm_num_layers=num_layers/min(num_layers);
    norm_cost=router_delay+norm_num_layers;
    norm_cost(find(topologies(:,3)<min_layers))=Inf;
%     norm_cost=router_delay;
    
    best_topology_idx=find(norm_cost==min(norm_cost),1);
    best_topology=topologies(best_topology_idx,:);
    best_topology_delay=router_delay(best_topology_idx);
    
%%%%%%%%%%%%%%PLOTS

    if show_plots
    
        horz_core_w=x(best_topology_idx);
        horz_core_h=core_area/x(best_topology_idx);
        vert_core_w=core_area/y(best_topology_idx);
        vert_core_h=y(best_topology_idx);

        horz_MC_w=horz_core_w*num_cols(best_topology_idx);
        horz_MC_h=MC_area/horz_MC_w;
        vert_MC_h=vert_core_h*num_rows(best_topology_idx);
        vert_MC_w=MC_area/vert_MC_h;

        MC_per_row=num_MC/num_layers(best_topology_idx)/num_rows(best_topology_idx);
        MC_per_col=num_MC/num_layers(best_topology_idx)/num_cols(best_topology_idx);

        figure
        MC_so_far=0.5;
        MC_prev=0;
        x_cursor=0;
        y_cursor=0;
        for i=1:best_topology(2)
            for j=1:best_topology(1)
                hold on
                rectangle('Position',[x_cursor,y_cursor,horz_core_w,horz_core_h],'FaceColor','g')
                x_cursor=x_cursor+horz_core_w;
            end
            x_cursor=0;
            y_cursor=y_cursor+horz_core_h;
            for k=1:fix(MC_so_far+MC_per_row)-MC_prev
                rectangle('Position',[x_cursor,y_cursor,horz_MC_w,horz_MC_h],'FaceColor','r')
                y_cursor=y_cursor+horz_MC_h;
            end
            MC_prev=fix(MC_so_far+MC_per_row);
            MC_so_far=MC_so_far+MC_per_row;
        end

        figure
        MC_so_far=0.5;
        MC_prev=0;
        x_cursor=0;
        y_cursor=0;
        for i=1:best_topology(1)
            for j=1:best_topology(2)
                hold on
                rectangle('Position',[x_cursor,y_cursor,vert_core_w,vert_core_h],'FaceColor','g')
                y_cursor=y_cursor+vert_core_h;
            end
            y_cursor=0;
            x_cursor=x_cursor+vert_core_w;
            for k=1:fix(MC_so_far+MC_per_col)-MC_prev
                rectangle('Position',[x_cursor,y_cursor,vert_MC_w,vert_MC_h],'FaceColor','r')
                x_cursor=x_cursor+vert_MC_w;
            end
            MC_prev=fix(MC_so_far+MC_per_col);
            MC_so_far=MC_so_far+MC_per_col;
        end    
    end

end

