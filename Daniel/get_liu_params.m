function [params] = get_liu_params(tmdc,approx)
%GET_LIU_PARAMS liefert eine materialspezifische Parametrisierung des
%Tight-Binding-Modells nach Liu et al. mit Spin-Bahn-Wechselwirkung 
% (Look-Up-Table).
%
%   Variable    Typ         Bedeutung
%   -----------------------------------------------------------------------
%   tmdc        string      Name des Materials. Folgende Möglichkeiten:
%                            - MoS2
%                            - WS2
%                            - MoSe2
%                            - WSe2
%                            - MoTe2
%                            - WTe2
%   approx      string      Parametersatz
%                            - LDA
%                            - GGA



% params = [E1 E2 t0 t1 t2 t11 t12 t22 r0 r1 r2 r11 r12 u0 u1 u2 u11 u12 ...
%               u22 lambda E12 Em10 E10];

switch approx
    case {'LDA', 'lda'}
        switch tmdc
            case 'MoS2'
                params = [0.820 1.931 -0.176 -0.101 0.531 0.084 0.169 0.070 ...
                    0.070 -0.252 0.084 0.019 0.093 -0.043 0.047 0.005 0.304 ...
                    -0.192 -0.162 0.073 4.840 1.395 3.176];
            case 'WS2'
                params = [0.905 2.167 -0.175 -0.090 0.611 0.043 0.181 0.008 ...
                    0.075 -0.282 0.127 0.001 0.114 -0.063 0.047 0.004 0.374 ...
                    -0.224 -0.177 0.211 5.473 1.526 3.667];
            case 'MoSe2'
                params = [0.715 1.687 -0.154 -0.134 0.437 0.124 0.119 0.072 ...
                    0.048 -0.248 0.090 0.066 0.045 -0.067 0.041 0.005 0.327 ...
                    -0.194 -0.151 0.091 4.296 1.128 2.862];
            case 'WSe2'
                params = [0.860 1.892 -0.152 -0.125 0.508 0.094 0.129 0.009 ...
                    0.044 -0.278 0.129 0.059 0.058 -0.090 0.039 0.001 0.392 ...
                    -0.224 -0.165 0.228 4.815 1.267 3.275];
            case 'MoTe2'
                params = [0.574 1.410 -0.148 -0.173 0.333 0.203 0.186 0.127 ...
                    0.007 -0.280 0.067 0.073 0.081 -0.054 0.008 0.037 0.145 ...
                    -0.078 0.035 0.107 3.991 0.798 2.918];
            case 'WTe2'
                params = [0.685 1.489 -0.124 -0.159 0.362 0.196 0.101 0.044 ... 
                    -0.009 -0.250 0.129 0.131 - 0.007 -0.086 0.012 -0.020   ...
                    0.361 -0.193 -0.129 0.237 4.412 1.004 3.347];
        end        
    case {'GGA', 'gga'}
        switch tmdc
            case 'MoS2'
                params = [0.683 1.707 -0.146 -0.114 0.506 0.085 0.162 0.073 ... 
                    0.060 -0.236 0.067 0.016 0.087 -0.038 0.046 0.001 0.266 ...
                    -0.176 -0.150 0.073 4.840 1.395 3.176];
            case 'WS2'
                params = [0.717 1.916 -0.152 -0.097 0.590 0.047 0.178 0.016 ...
                    0.069 -0.261 0.107 -0.003 0.109 -0.054 0.045 0.002 0.325...
                    -0.206 -0.163 0.211 5.473 1.526 3.667];                
            case 'MoSe2'
                params = [0.684 1.546 -0.146 -0.130 0.432 0.144 0.117 0.075 ... 
                    0.039 -0.209 0.069 0.052 0.060 -0.042 0.035 0.008 0.272 ...
                    -0.172 -0.150 0.091 4.296 1.128 2.862];
            case 'WSe2'
                params = [0.728 1.655 -0.146 -0.124 0.507 0.117 0.127 0.015 ...
                    0.036 -0.234 0.107 0.044 0.075 -0.051 0.032 0.007 0.329 ...
                    -0.202 -0.164 0.228 4.815 1.267 3.275];
            case 'MoTe2'
                params = [0.588 1.303 -0.226 -0.234 0.036 0.400 0.098 0.017 ...
                    0.003 -0.025 -0.169 0.082 0.051 0.057 0.103 0.187 -0.045...
                    -0.141 0.087 0.107 3.991 0.798 2.918];
            case 'WTe2'
                params = [0.697 1.380 -0.109 -0.164 0.368 0.204 0.093 0.038 ...
                    -0.015 -0.209 0.107 0.115 0.009 -0.066 0.011 -0.013     ...
                    0.312 -0.177 -0.132 0.237 4.412 1.004 3.347];
        end
end