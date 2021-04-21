  YAML::Node default_config;
  default_config["species"]["name"] = "species";
  default_config["species"]["num"] = "0";
  default_config["species"]["diameter"] = "1";
  default_config["species"]["length"] = "0";
  default_config["species"]["insertion_type"] = "random";
  default_config["species"]["insert_file"] = "none";
  default_config["species"]["overlap"] = "false";
  default_config["species"]["draw_type"] = "orientation";
  default_config["species"]["color"] = "0";
  default_config["species"]["posit_flag"] = "false";
  default_config["species"]["spec_flag"] = "false";
  default_config["species"]["n_posit"] = "100";
  default_config["species"]["n_spec"] = "100";
  default_config["rigid_filament"]["max_length"] = "500";
  default_config["rigid_filament"]["min_length"] = "5";
  default_config["rigid_filament"]["stationary_flag"] = "false";
  default_config["rigid_filament"]["constrain_motion_flag"] = "false";
  default_config["rigid_filament"]["packing_fraction"] = "-1";
  default_config["rigid_filament"]["n_equil"] = "0";
  default_config["filament"]["packing_fraction"] = "-1";
  default_config["filament"]["persistence_length"] = "400";
  default_config["filament"]["perlen_ratio"] = "-1";
  default_config["filament"]["polydispersity_flag"] = "false";
  default_config["filament"]["max_length"] = "500";
  default_config["filament"]["min_length"] = "5";
  default_config["filament"]["min_bond_length"] = "1.5";
  default_config["filament"]["driving_factor"] = "0";
  default_config["filament"]["n_equil"] = "0";
  default_config["filament"]["nematic_driving"] = "false";
  default_config["filament"]["nematic_driving_freq"] = "0";
  default_config["filament"]["peclet_number"] = "-1";
  default_config["filament"]["flexure_number"] = "-1";
  default_config["filament"]["radius_of_curvature"] = "-1";
  default_config["filament"]["intrinsic_curvature"] = "0";
  default_config["filament"]["intrinsic_curvature_sig"] = "0";
  default_config["filament"]["randomize_intrinsic_curvature_handedness"] = "false";
  default_config["filament"]["intrinsic_curvature_min"] = "0";
  default_config["filament"]["highlight_handedness"] = "false";
  default_config["filament"]["highlight_curvature"] = "false";
  default_config["filament"]["draw_center_of_curvature"] = "false";
  default_config["filament"]["error_analysis"] = "false";
  default_config["filament"]["theta_analysis"] = "false";
  default_config["filament"]["lp_analysis"] = "false";
  default_config["filament"]["global_order_analysis"] = "false";
  default_config["filament"]["polar_order_analysis"] = "false";
  default_config["filament"]["polar_order_n_bins"] = "100";
  default_config["filament"]["polar_order_contact_cutoff"] = "3";
  default_config["filament"]["polar_order_width"] = "10";
  default_config["filament"]["msd_analysis"] = "false";
  default_config["filament"]["curvature_cluster_analysis"] = "false";
  default_config["filament"]["cluster_lifetime_min"] = "100";
  default_config["filament"]["cluster_by_handedness"] = "false";
  default_config["filament"]["spiral_init_flag"] = "false";
  default_config["filament"]["spiral_analysis"] = "false";
  default_config["filament"]["spiral_number_fail_condition"] = "0";
  default_config["filament"]["orientation_corr_analysis"] = "false";
  default_config["filament"]["orientation_corr_n_steps"] = "1000";
  default_config["filament"]["crossing_analysis"] = "false";
  default_config["filament"]["flocking_analysis"] = "false";
  default_config["filament"]["flock_polar_min"] = "0.5";
  default_config["filament"]["flock_contact_min"] = "0.5";
  default_config["filament"]["highlight_flock"] = "false";
  default_config["filament"]["flock_color_int"] = "1.57";
  default_config["filament"]["flock_color_ext"] = "4.71";
  default_config["filament"]["number_fluctuation_analysis"] = "false";
  default_config["filament"]["number_fluctuation_boxes"] = "6";
  default_config["filament"]["number_fluctuation_centers"] = "10";
  default_config["filament"]["in_out_analysis"] = "false";
  default_config["filament"]["drive_from_bond_center"] = "true";
  default_config["filament"]["flagella_flag"] = "false";
  default_config["filament"]["flagella_freq"] = "1";
  default_config["filament"]["flagella_period"] = "2";
  default_config["filament"]["flagella_amplitude"] = "1";
  default_config["filament"]["friction_ratio"] = "2";
  default_config["filament"]["dynamic_instability_flag"] = "false";
  default_config["filament"]["force_induced_catastrophe_flag"] = "false";
  default_config["filament"]["optical_trap_flag"] = "false";
  default_config["filament"]["optical_trap_spring"] = "20";
  default_config["filament"]["optical_trap_fixed"] = "false";
  default_config["filament"]["cilia_trap_flag"] = "false";
  default_config["filament"]["fic_factor"] = "0.828";
  default_config["filament"]["f_shrink_to_grow"] = "0.017";
  default_config["filament"]["f_shrink_to_pause"] = "0.0";
  default_config["filament"]["f_pause_to_grow"] = "0.0";
  default_config["filament"]["f_pause_to_shrink"] = "0.0";
  default_config["filament"]["f_grow_to_pause"] = "0.0";
  default_config["filament"]["f_grow_to_shrink"] = "0.00554";
  default_config["filament"]["v_poly"] = "0.44";
  default_config["filament"]["v_depoly"] = "0.793";
  default_config["filament"]["custom_set_tail"] = "false";
  default_config["filament"]["reference_frame_flag"] = "false";
  default_config["br_bead"]["driving_factor"] = "0";
  default_config["br_bead"]["driving_torque"] = "0";
  default_config["br_bead"]["density"] = "-1";
  default_config["br_bead"]["chiral_handedness"] = "0";
  default_config["br_bead"]["rotational_noise"] = "1";
  default_config["br_bead"]["translational_noise"] = "1";
  default_config["br_bead"]["randomize_handedness"] = "false";
  default_config["br_bead"]["highlight_handedness"] = "false";
  default_config["br_bead"]["alignment_interaction"] = "false";
  default_config["br_bead"]["alignment_torque"] = "0";
  default_config["br_bead"]["packing_fraction"] = "-1";
  default_config["br_bead"]["draw_shape"] = "sphere";
  default_config["spherocylinder"]["diffusion_analysis"] = "false";
  default_config["spherocylinder"]["n_diffusion_samples"] = "1";
  default_config["spherocylinder"]["midstep"] = "false";
  default_config["spindle"]["n_filaments_bud"] = "0";
  default_config["spindle"]["n_filaments_mother"] = "0";
  default_config["spindle"]["alignment_potential"] = "false";
  default_config["spindle"]["k_spring"] = "1000";
  default_config["spindle"]["k_align"] = "0";
  default_config["spindle"]["spring_length"] = "0";
  default_config["spindle"]["spb_diameter"] = "5";
  default_config["spindle"]["nuc_site_insertion"] = "random";
  default_config["crosslink"]["anchors"][0]["velocity_s"] = "0";
  default_config["crosslink"]["anchors"][1]["velocity_s"] = "0";
  default_config["crosslink"]["anchors"][0]["velocity_d"] = "0";
  default_config["crosslink"]["anchors"][1]["velocity_d"] = "0";
  default_config["crosslink"]["anchors"][0]["color"] = "0";
  default_config["crosslink"]["anchors"][1]["color"] = "0";
  default_config["crosslink"]["anchors"][0]["bind_file"] = "none";
  default_config["crosslink"]["anchors"][1]["bind_file"] = "none";
  default_config["crosslink"]["anchors"][0]["k_on_s"] = "10";
  default_config["crosslink"]["anchors"][1]["k_on_s"] = "10";
  default_config["crosslink"]["anchors"][0]["k_off_s"] = "2";
  default_config["crosslink"]["anchors"][1]["k_off_s"] = "2";
  default_config["crosslink"]["anchors"][0]["k_on_d"] = "10";
  default_config["crosslink"]["anchors"][1]["k_on_d"] = "10";
  default_config["crosslink"]["anchors"][0]["k_off_d"] = "2";
  default_config["crosslink"]["anchors"][1]["k_off_d"] = "2";
  default_config["crosslink"]["concentration"] = "0";
  default_config["crosslink"]["begin_with_bound_crosslinks"] = "0";
  default_config["crosslink"]["use_binding_volume"] = "true";
  default_config["crosslink"]["infinite_reservoir_flag"] = "false";
  default_config["crosslink"]["bind_site_density"] = "1";
  default_config["crosslink"]["static_flag"] = "false";
  default_config["crosslink"]["diffusion_s"] = "0";
  default_config["crosslink"]["diffusion_d"] = "0";
  default_config["crosslink"]["k_on_s"] = "10";
  default_config["crosslink"]["k_off_s"] = "2";
  default_config["crosslink"]["k_on_d"] = "10";
  default_config["crosslink"]["k_off_d"] = "2";
  default_config["crosslink"]["energy_dep_factor"] = "0";
  default_config["crosslink"]["force_dep_length"] = "0";
  default_config["crosslink"]["polar_affinity"] = "1";
  default_config["crosslink"]["k_spring"] = "10";
  default_config["crosslink"]["k_spring_compress"] = "-1.";
  default_config["crosslink"]["f_stall"] = "100";
  default_config["crosslink"]["force_dep_vel_flag"] = "true";
  default_config["crosslink"]["k_align"] = "0";
  default_config["crosslink"]["rest_length"] = "0";
  default_config["crosslink"]["step_direction"] = "0";
  default_config["crosslink"]["tether_draw_type"] = "orientation";
  default_config["crosslink"]["tether_diameter"] = "0.5";
  default_config["crosslink"]["tether_color"] = "3.1416";
  default_config["crosslink"]["minus_end_pausing"] = "false";
  default_config["crosslink"]["plus_end_pausing"] = "false";
  default_config["crosslink"]["r_capture"] = "5";
  default_config["crosslink"]["lut_grid_num"] = "256";
  default_config["receptor"]["component"] = "cortex";
  default_config["receptor"]["concentration"] = "-1";
  default_config["receptor"]["induce_catastrophe"] = "false";
  default_config["seed"] = "7859459105545";
  default_config["n_runs"] = "1";
  default_config["n_random"] = "1";
  default_config["run_name"] = "sc";
  default_config["n_dim"] = "3";
  default_config["n_periodic"] = "0";
  default_config["boundary"] = "0";
  default_config["system_radius"] = "100";
  default_config["n_steps"] = "1000000";
  default_config["i_step"] = "0";
  default_config["prev_step"] = "0";
  default_config["delta"] = "0.001";
  default_config["dynamic_timestep"] = "false";
  default_config["dynamic_timestep_ramp"] = "0.001";
  default_config["graph_flag"] = "false";
  default_config["n_graph"] = "1000";
  default_config["graph_diameter"] = "0";
  default_config["invert_background"] = "false";
  default_config["draw_boundary"] = "true";
  default_config["load_checkpoint"] = "false";
  default_config["checkpoint_run_name"] = "sc";
  default_config["n_load"] = "0";
  default_config["movie_flag"] = "false";
  default_config["movie_directory"] = "frames";
  default_config["time_analysis"] = "false";
  default_config["bud_height"] = "680";
  default_config["bud_radius"] = "300";
  default_config["lj_epsilon"] = "1";
  default_config["wca_eps"] = "1";
  default_config["wca_sig"] = "1";
  default_config["ss_a"] = "1";
  default_config["ss_rs"] = "1.5";
  default_config["ss_eps"] = "1";
  default_config["f_cutoff"] = "2000";
  default_config["constant_pressure"] = "false";
  default_config["constant_volume"] = "false";
  default_config["target_pressure"] = "0";
  default_config["target_radius"] = "100";
  default_config["pressure_time"] = "100";
  default_config["compressibility"] = "1";
  default_config["zero_temperature"] = "false";
  default_config["thermo_flag"] = "false";
  default_config["n_thermo"] = "1000";
  default_config["insert_radius"] = "-1";
  default_config["interaction_flag"] = "true";
  default_config["remove_duplicate_interactions"] = "false";
  default_config["coarse_grained_mesh_interactions"] = "false";
  default_config["mesh_coarsening"] = "2";
  default_config["species_insertion_failure_threshold"] = "10000";
  default_config["species_insertion_reattempt_threshold"] = "10";
  default_config["uniform_crystal"] = "false";
  default_config["n_steps_equil"] = "0";
  default_config["n_steps_target"] = "100000";
  default_config["static_particle_number"] = "false";
  default_config["checkpoint_from_spec"] = "false";
  default_config["potential"] = "wca";
  default_config["soft_potential_mag"] = "10";
  default_config["soft_potential_mag_target"] = "-1";
  default_config["like_like_interactions"] = "true";
  default_config["auto_graph"] = "false";
  default_config["local_order_analysis"] = "false";
  default_config["local_order_width"] = "50";
  default_config["local_order_bin_width"] = "0.5";
  default_config["local_order_n_analysis"] = "100";
  default_config["density_analysis"] = "0";
  default_config["density_bin_width"] = "0.1";
  default_config["density_com_only"] = "false";
  default_config["overlap_analysis"] = "false";
  default_config["highlight_overlaps"] = "false";
  default_config["reduced"] = "false";
  default_config["reload_reduce_switch"] = "false";
  default_config["checkpoint_flag"] = "false";
  default_config["n_checkpoint"] = "10000";
  default_config["knockout_xlink"] = "false";
  default_config["no_midstep"] = "false";
