
%runPilotsConDuration
%
%     author: steeve laquitaine
%       date: 151107
%    purpose: Wrapper file to run pilot psyphy.
%             Test different contrasts and stim durations
%
%
%1) [stim: 8Hz flk,0.3s]
%
%s025: det thrsh : 0.0009 (0.09%)
%
%priors 80 and 40 deg

%test con [0.0009 0.0018 0.0036] 
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd040_mean225_con000090000200003_loc36_t100_075_032perCon_151107','responseObject=dot','responseDevice=0')
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd080_mean225_con000010000200003_loc36_t100_075_034perCon_151107','responseObject=dot','responseDevice=0')
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd080_mean225_con000010000200003_loc36_t100_075_034perCon_151107','responseObject=dot','responseDevice=0')
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd040_mean225_con000090000200003_loc36_t100_075_032perCon_151107','responseObject=dot','responseDevice=0')
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd080_mean225_con000010000200003_loc36_t100_075_034perCon_151107','responseObject=dot','responseDevice=0')

%test con [0.0014 0.0027 0.0054]
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd040_mean225_con00014_00027_00054_loc36_t100_075_032perCon_151107','responseObject=dot','responseDevice=0')
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd080_mean225_con00014_00027_00054_loc36_t100_075_034perCon_151107','responseObject=dot','responseDevice=0')
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd040_mean225_con00014_00027_00054_loc36_t100_075_032perCon_151107','responseObject=dot','responseDevice=0')
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd080_mean225_con00014_00027_00054_loc36_t100_075_034perCon_151107','responseObject=dot','responseDevice=0')

%test con 1
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd080_mean225_con1_loc36_t100_075_034perCon_151107','responseObject=dot','responseDevice=0')
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd040_mean225_con1_loc36_t100_075_032perCon_151107','responseObject=dot','responseDevice=0')



%another subject XX

%test con [0.0009 0.0018 0.0036] 
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd040_mean225_con00009_00018_00036_loc36_t100_075_032perCon_151107','responseObject=dot','responseDevice=0')
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd080_mean225_con00009_00018_00036_loc36_t100_075_034perCon_151107','responseObject=dot','responseDevice=0')
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd080_mean225_con00009_00018_00036_loc36_t100_075_034perCon_151107','responseObject=dot','responseDevice=0')
%slTaskLocEstPilotConDuration('psychophysics','displayName=testPsyphy','params=steeve_exp04_metho_Pstd080_mean225_con00009_00018_00036_loc36_t100_075_034perCon_151107','responseObject=dot','responseDevice=0')


