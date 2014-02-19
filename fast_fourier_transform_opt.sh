#!/bin/bash
# ./bin/test_fourier_transform hpce.fs1910.fast_fourier_transform_opt
./bin/time_fourier_transform hpce.fs1910.fast_fourier_transform_opt 1 | grep 881743 | cut -d\  -f2- 
# ./bin/time_fourier_transform hpce.fs1910.fast_fourier_transform_opt 2 | grep 881743 | cut -d\  -f2-
# ./bin/time_fourier_transform hpce.fs1910.fast_fourier_transform_opt 3 | grep 881743 | cut -d\  -f2-
# ./bin/time_fourier_transform hpce.fs1910.fast_fourier_transform_opt 4 | grep 881743 | cut -d\  -f2-


