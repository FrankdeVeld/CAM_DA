clc
!wsl g++ -o DATTS DATimeToSafety.cpp -I /usr/local/include -L /usr/local/lib -ldace
!wsl ./DATTS
disp("Done")
