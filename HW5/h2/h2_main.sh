g++ h2_main_x.cpp ../util.cpp -o h2_main_x -O2 -larmadillo -llapack -lblas
./h2_main_x H1 H2

g++ h2_main_y.cpp ../util.cpp -o h2_main_y -O2 -larmadillo -llapack -lblas
./h2_main_y H1 H2

g++ h2_main_z.cpp ../util.cpp -o h2_main_z -O2 -larmadillo -llapack -lblas
./h2_main_z H1 H2
