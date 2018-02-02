function [ out_res, in_cap, par_cap, wire_res, wire_cap ] = calc_cacti_parameters( tech, vdd_real, verbose )

    %%%%%%%UNITS%%%%
    %A
    %V
    %um
    %F

%     tech=45;
%     vdd_real=1;

    F_sz_um=tech/1000;
    Cpolywire=0;
    c_junc_sidewall=0.25e-15; %F/micron
    BULK_CU_RESISTIVITY = 0.018; %ohm-micron
    CU_RESISTIVITY = 0.022; %ohm-micron
    PERMITTIVITY_FREE_SPACE = 8.854e-18; %F/micron
    fringe_cap=0.115*1e-15;
    vert_dielectric_constant=3.9;
    miller_value=1.5;
    alpha_scatter = 1;
    aspect_ratio=[2.0 2.0 2.2];
    wire_pitch=[2.5 4 8]*F_sz_um;

    if tech==32
        %32 nm tech
        % F_sz_um=0.032; %micron
        n_to_p_eff_curr_drv_ratio=2.41;
        nmos_effective_resistance_multiplier=1.49;
        vdd=0.9; %V
        v_th=0.21835;
        I_on_n=2211.7*1e-6*(vdd_real-v_th)/(vdd-v_th); %A/micron
        c_g_ideal=5.34e-16; %F/micron
        c_fringe=0.04e-15; %F/micron
        Lphy=0.013; %micron
        c_junc_area=1e-15; %F/micron^2


        % wire_pitch=[2.5 4 8]*F_sz_um;
        % aspect_ratio=[2.0 2.0 2.2];
        barrier_thickness = 0.003;
        % alpha_scatter = 1;
        % wire_width=wire_pitch./2;
        % wire_thickness=aspect_ratio.*wire_width;
        % wire_spacing=wire_pitch-wire_width;
        % dishing_thickness = [0 0 0.1].*wire_thickness;
        ild_thickness=[0.21 0.21 0.385];
        % miller_value=[1.5 1.5 1.5];
        horiz_dielectric_constant=2.214;
        % vert_dielectric_constant=[3.9 3.9 3.9];
        % fringe_cap=0.115*1e-15;
    elseif tech==45
        %45 nm tech
        % F_sz_um=0.045; %micron
        n_to_p_eff_curr_drv_ratio=2.41;
        nmos_effective_resistance_multiplier=1.51;
        vdd=1; %V
        v_th=0.18035;
        I_on_n=2046.6*1e-6*(vdd_real-v_th)/(vdd-v_th); %A/micron
        c_g_ideal=6.78e-16; %F/micron
        c_fringe=0.05e-15; %F/micron
        Lphy=0.018; %micron
        c_junc_area=1e-15; %F/micron^2


        % wire_pitch=[2.5 4 8]*F_sz_um;
        % aspect_ratio=[2.0 2.0 2.2];
        barrier_thickness = 0.004;
        % alpha_scatter = 1;
        % wire_width=wire_pitch./2;
        % wire_thickness=aspect_ratio.*wire_width;
        % wire_spacing=wire_pitch-wire_width;
        % dishing_thickness = [0 0 0.1].*wire_thickness;
        ild_thickness=[0.315 0.315 0.55];
        % miller_value=[1.5 1.5 1.5];
        horiz_dielectric_constant=2.46;
        % vert_dielectric_constant=[3.9 3.9 3.9];
        % fringe_cap=0.115*1e-15;
    elseif tech==65
        %65 nm tech
        % F_sz_um=0.065;
        n_to_p_eff_curr_drv_ratio=2.41;
        nmos_effective_resistance_multiplier=1.50;
        vdd=1.1;
        v_th=0.19491;
        I_on_n=1197.2*1e-6*(vdd_real-v_th)/(vdd-v_th);
        c_g_ideal=4.69e-16;
        c_fringe=0.077e-15;
        Lphy=0.025;
        c_junc_area=1e-15;

        % wire_pitch=[2.5 4 8]*F_sz_um;
        % aspect_ratio=[2.0 2.0 2.2];
        barrier_thickness = 0.006;
        % alpha_scatter = 1;
        % wire_width=wire_pitch./2;
        % wire_thickness=aspect_ratio.*wire_width;
        % wire_spacing=wire_pitch-wire_width;
        % dishing_thickness = [0 0 0.1].*wire_thickness;
        ild_thickness=[0.405 0.405 0.77];
        % miller_value=[1.5 1.5 1.5];
        horiz_dielectric_constant=2.734;
        % vert_dielectric_constant=[3.9 3.9 3.9];
        % fringe_cap=0.115*1e-15;
    elseif tech==90
        %90 nm tech
        % F_sz_um=0.090;
        n_to_p_eff_curr_drv_ratio=2.45;
        nmos_effective_resistance_multiplier=1.54;
        vdd=1.2;
        v_th=0.23707;
        I_on_n=1076.9*1e-6*(vdd_real-v_th)/(vdd-v_th);
        c_g_ideal=6.64e-16;
        c_fringe=0.08e-15;
        Lphy=0.037;
        c_junc_area=1e-15;

        % wire_pitch=[2.5 4 8]*F_sz_um;
        % aspect_ratio=[2.0 2.0 2.2];
        barrier_thickness = 0.008;
        % alpha_scatter = 1;
        % wire_width=wire_pitch./2;
        % wire_thickness=aspect_ratio.*wire_width;
        % wire_spacing=wire_pitch-wire_width;
        % dishing_thickness = [0 0 0.1].*wire_thickness;
        ild_thickness=[0.48 0.48 1.1];
        % miller_value=[1.5 1.5 1.5];
        horiz_dielectric_constant=3.038;
        % vert_dielectric_constant=[3.9 3.9 3.9];
        % fringe_cap=0.115*1e-15;
    else
        display('Invalid Technology Node')
        return
    end

    % %180 nm tech
    % F_sz_um=0.180;
    % n_to_p_eff_curr_drv_ratio=2.45;
    % nmos_effective_resistance_multiplier=1.54;
    % vdd=1.5;
    % I_on_n=750*1e-6;
    % c_g_ideal=2*6.64e-16;
    % c_fringe=2*0.08e-15;
    % Lphy=0.12;
    % c_junc_area=2*1e-15;

    %output resistance

    min_w_nmos=(3/2)*F_sz_um;
    min_w_pmos=min_w_nmos*n_to_p_eff_curr_drv_ratio;

    restrans_n=nmos_effective_resistance_multiplier*vdd_real/I_on_n;
    restrans_p=restrans_n*n_to_p_eff_curr_drv_ratio;

    tr_R_on_n=restrans_n/min_w_nmos;
    tr_R_on_p=restrans_p/min_w_pmos;

    out_res=(tr_R_on_n+tr_R_on_p)/2;

    %input capacitance

    c_overlap=0.2*c_g_ideal;

    in_cap=c_g_ideal + c_overlap + 3*c_fringe*(min_w_nmos+min_w_pmos) + Lphy*Cpolywire;

    % parasitic capacitance

    cell_h_def=50*F_sz_um;
    HPOWERRAIL=2*F_sz_um;
    MIN_GAP_BET_P_AND_N_DIFFS = 5 * F_sz_um;
    w_poly_contact = F_sz_um;
    spacing_poly_to_contact = F_sz_um;

    drain_C_metal_connecting_folded_tr=0;
    h_tr_region=cell_h_def-2*HPOWERRAIL;
    ratio_p_to_n = 2/(2+1);
    w_folded_tr_n = (1 - ratio_p_to_n) * (h_tr_region - MIN_GAP_BET_P_AND_N_DIFFS);
    w_folded_tr_p = (ratio_p_to_n) * (h_tr_region - MIN_GAP_BET_P_AND_N_DIFFS);

    num_folded_tr_n = ceil(min_w_nmos/w_folded_tr_n);
    num_folded_tr_p = ceil(min_w_pmos/w_folded_tr_p);

    if num_folded_tr_n < 2
        w_folded_tr_n = min_w_nmos;
    end

    if num_folded_tr_p < 2
        w_folded_tr_p = min_w_pmos;
    end

    total_drain_w = w_poly_contact + 2*spacing_poly_to_contact;

    drain_h_for_sidewall_n = w_folded_tr_n;
    drain_h_for_sidewall_p = w_folded_tr_p;

    total_drain_height_for_cap_wrt_gate_n = w_folded_tr_n;
    total_drain_height_for_cap_wrt_gate_p = w_folded_tr_p;

    drain_C_area_n = c_junc_area*total_drain_w*w_folded_tr_n;
    drain_C_area_p = c_junc_area*total_drain_w*w_folded_tr_p;

    drain_C_sidewall_n = c_junc_sidewall * (drain_h_for_sidewall_n + 2 * total_drain_w);
    drain_C_sidewall_p = c_junc_sidewall * (drain_h_for_sidewall_p + 2 * total_drain_w);

    drain_C_wrt_gate_n = (2*c_fringe + 2*c_overlap)*total_drain_height_for_cap_wrt_gate_n;
    drain_C_wrt_gate_p = (2*c_fringe + 2*c_overlap)*total_drain_height_for_cap_wrt_gate_p;

    par_cap = drain_C_area_n + drain_C_area_p + drain_C_sidewall_n + drain_C_sidewall_p +drain_C_wrt_gate_n + drain_C_wrt_gate_p + 2*drain_C_metal_connecting_folded_tr;

    %wire_resistance

    wire_width=wire_pitch./2;
    wire_thickness=aspect_ratio.*wire_width;
    wire_spacing=wire_pitch-wire_width;
    dishing_thickness = [0 0 0.1].*wire_thickness;

    wire_res=alpha_scatter*CU_RESISTIVITY./((wire_thickness-barrier_thickness-dishing_thickness).*(wire_width-2*barrier_thickness));

    %wire_capacitance

    vertical_cap = 2 * PERMITTIVITY_FREE_SPACE .* vert_dielectric_constant .* wire_width ./ ild_thickness;
    sidewall_cap = 2 * PERMITTIVITY_FREE_SPACE .* miller_value .* horiz_dielectric_constant .* wire_thickness ./ wire_spacing;
    wire_cap = vertical_cap + sidewall_cap + fringe_cap;

    %print results
    if verbose
        display(['output resistance = ' num2str(out_res) ' ohms'])
        display(['input capacitance = ' num2str(in_cap*1e15) ' fF'])
        display(['paracitic capacitance = ' num2str(par_cap*1e15) ' fF'])
        display(['wire resistance per unit length = ' num2str(wire_res) ' ohm/micron'])
        display(['wire capacitance per unit length = ' num2str(wire_cap*1e15) ' fF/micron'])
    end
end