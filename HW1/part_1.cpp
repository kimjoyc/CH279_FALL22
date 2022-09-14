#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <armadillo>



using std::string;
using std::vector;

class Atom {
private:
	int atomic_number;
	double x;
	double y;
    double z;



public:
	Atom(int input_atomic_number, double input_x, double input_y,double input_z) {
		atomic_number = input_atomic_number;
		x = input_x;
		y = input_y;
        z = input_z;


	}

	void print_info() {
		std::cout << "Atomic Number : " << atomic_number <<
		", x coordinate : " << x << ", y coordinate : " << y <<  ", z coordinate : " << z
		<< std::endl;
	}


	double lennard_jones(double x_j, double y_j,double z_j, double gold_ep_ij,double gold_sig_ij) {
        double x_ij= pow((x-x_j),2);
        double y_ij= pow((y-y_j),2);
        double z_ij= pow((z-z_j),2);
        double r_ij= sqrt(x_ij+y_ij+z_ij);

		
		return gold_ep_ij*(pow(gold_sig_ij/r_ij,12) - 2*pow(gold_sig_ij/r_ij,6));	
	}


};

void readfile(vector<Atom> &atoms, string &filename) {
	std::ifstream infile(filename);
	if (infile.is_open()) {

        int atomic_number;
        double x;
        double y;
        double z;

		while (infile >> atomic_number >> x >> y >> z) {




			Atom new_atom(atomic_number, x, y, z);
			atoms.push_back(new_atom);
			
		}
		infile.close();

	} else {
		throw std::invalid_argument("Can't open file to read.");
	}
}





int main(int argc, char* argv[]) {

	vector<Atom> atoms;
	string filename = argv[1];



	try {
		readfile(atoms, filename);

	}
	catch (std::invalid_argument &e) {
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		std::cout << "Something wrong happened!" << std::endl;
		return EXIT_FAILURE;
	}

	
	for (auto x : atoms) {
		x.print_info();
		arma::mat input_1 = {x.lennard_jones(0,0,1,2.951,5.29)};
		input_1.print();

	}
	return EXIT_SUCCESS;
}
