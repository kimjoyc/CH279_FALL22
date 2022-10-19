g++ part1.cpp -o part1 util.cpp -O2 -larmadillo -llapack -lblas

# row 1
./part1 H1_H1_1 H1_H1_2 H1_H1_3
./part1 H1_H2_1 H1_H2_2 H1_H2_3

#row 2
./part1 H2_H1_1 H2_H1_2 H2_H1_3
./part1 H2_H2_1 H2_H2_2 H2_H2_3


g++ part2.cpp -o part2 util.cpp -O2 -larmadillo -llapack -lblas
./part2 


