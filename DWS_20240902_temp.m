[X,Y,T,U,V,Xdrift,Ydrift] = readPIV_bin("G:\Combined_Contractility_123\output\cal520_bluebeads_PT_20240618_10s_20240618165855_plateB_Well__G_008_r_0004_c_0004_BLUE_BEADS\piv\deformations_pass1.bin");

%'cal520_bluebeads_PT_20240618_10s_20240618165855_plateB_Well__B_006_r_0005_c_0005_BLUE_BEADS_tfm_trace_peaks'
%cal520_bluebeads_PT_20240618_10s_20240618165855_plateB_Well__G_008_r_0004_c_0004_BLUE_BEADS_tfm_trace_peaks
%max(max(U))

[MU,IU] = max(abs(U),[],"all");

[MV,IV] = max(abs(V),[],"all");