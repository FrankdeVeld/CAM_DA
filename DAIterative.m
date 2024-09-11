%clc
!wsl g++ -o DAIterative DAGreedyTestRK78Iterative.cpp -I /usr/local/include -L /usr/local/lib -ldace
!wsl ./DAIterative
disp("Done")
