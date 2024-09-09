clc
!wsl g++ -o DAItVal DAIterativeValidation.cpp -I /usr/local/include -L /usr/local/lib -ldace
!wsl ./DAItVal
disp("Done")
