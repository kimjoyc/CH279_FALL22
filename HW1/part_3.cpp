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

    //finite diff central diff approx for xyz
    double lennard_jones_center_diff(double x_j, double y_j,double z_j,double gold_ep_ij,double gold_sig_ij,double h) {
        double x_ij= pow((x-x_j),2);
        double y_ij= pow((y-y_j),2);
        double z_ij= pow((z-z_j),2);

        //x shift
        double x_ij_der_left= pow((x-x_j-h),2);
        double x_ij_der_right= pow((x-x_j+h),2);
        //x shift calc
        double r_ij_x_left= sqrt(x_ij_der_left+y_ij+z_ij);
        double r_ij_x_right= sqrt(x_ij_der_right+y_ij+z_ij);
        //energy calc w/shift for x
		double lj_energy_left_x = gold_ep_ij*(pow(gold_sig_ij/r_ij_x_left,12) - 2*pow(gold_sig_ij/r_ij_x_left,6));
        double lj_energy_right_x= gold_ep_ij*(pow(gold_sig_ij/r_ij_x_right,12) - 2*pow(gold_sig_ij/r_ij_x_right,6));	
        //cent x 
        double central_diff_x= -(lj_energy_right_x-lj_energy_left_x)/(2*h);


        //y shift
        double y_ij_der_left= pow((y-y_j-h),2);
        double y_ij_der_right= pow((y-y_j+h),2);
        //yshift cacl
        double r_ij_y_left= sqrt(y_ij_der_left+x_ij+z_ij);
        double r_ij_y_right= sqrt(y_ij_der_right+x_ij+z_ij);
        //energy calc w/shift for y
		double lj_energy_left_y = gold_ep_ij*(pow(gold_sig_ij/r_ij_y_left,12) - 2*pow(gold_sig_ij/r_ij_y_left,6));
        double lj_energy_right_y= gold_ep_ij*(pow(gold_sig_ij/r_ij_y_right,12) - 2*pow(gold_sig_ij/r_ij_y_right,6));	
        //ycent calc
        double central_diff_y= -(lj_energy_right_y-lj_energy_left_y)/(2*h);
        
        
        //z shift
        double z_ij_der_left= pow((z-z_j-h),2);
        double z_ij_der_right= pow((z-z_j+h),2);
        //zshift cacl
        double r_ij_z_left= sqrt(z_ij_der_left+x_ij+y_ij);
        double r_ij_z_right= sqrt(z_ij_der_right+x_ij+y_ij);
        //energy calc w/shift for z
		double lj_energy_left_z = gold_ep_ij*(pow(gold_sig_ij/r_ij_z_left,12) - 2*pow(gold_sig_ij/r_ij_z_left,6));
        double lj_energy_right_z= gold_ep_ij*(pow(gold_sig_ij/r_ij_z_right,12) - 2*pow(gold_sig_ij/r_ij_z_right,6));	
        //zcent calc
        double central_diff_z= -(lj_energy_right_z-lj_energy_left_z)/(2*h);

        double central_diff_xyz=sqrt(pow(central_diff_x,2)+pow(central_diff_y,2)+pow(central_diff_z,2));

        return central_diff_xyz; 
	}


    //finite diff forward diff approx

    double lennard_jones_forward_diff(double x_j, double y_j,double z_j, double gold_ep_ij,double gold_sig_ij,double h) {
        double x_ij= pow((x-x_j),2);
        double y_ij= pow((y-y_j),2);
        double z_ij= pow((z-z_j),2);
        
        double r_ij= sqrt(x_ij+y_ij+z_ij);
        double lj_energy = gold_ep_ij*(pow(gold_sig_ij/r_ij,12) - 2*pow(gold_sig_ij/r_ij,6));

        //x shift
        double x_ij_der_right_= pow((x-x_j+h),2);
        //x shift calc
        double r_ij_x_right_= sqrt(x_ij_der_right_+y_ij+z_ij);
        //energy calc w/shift for x
        double lj_energy_right_x_= gold_ep_ij*(pow(gold_sig_ij/r_ij_x_right_,12) - 2*pow(gold_sig_ij/r_ij_x_right_,6));	
        //finite forward x 
        double forward_diff_x= -(lj_energy_right_x_-lj_energy)/(h);


        //y shift
        double y_ij_der_right_= pow((y-y_j+h),2);
        //yshift calc
        double r_ij_y_right_= sqrt(y_ij_der_right_+y_ij+z_ij);
        //energy calc w/shift for y
        double lj_energy_right_y_= gold_ep_ij*(pow(gold_sig_ij/r_ij_y_right_,12) - 2*pow(gold_sig_ij/r_ij_y_right_,6));	
        //finite forward y
        double forward_diff_y= -(lj_energy_right_y_-lj_energy)/(h);


        //y shift
        double z_ij_der_right_= pow((z-z_j+h),2);
        //yshift calc
        double r_ij_z_right_= sqrt(z_ij_der_right_+y_ij+z_ij);
        //energy calc w/shift for z
        double lj_energy_right_z_= gold_ep_ij*(pow(gold_sig_ij/r_ij_z_right_,12) - 2*pow(gold_sig_ij/r_ij_z_right_,6));	
        //finite forward y
        double forward_diff_z= -(lj_energy_right_z_-lj_energy)/(h);


        double forward_diff_xyz=sqrt(pow(forward_diff_x,2)+pow(forward_diff_y,2)+pow(forward_diff_z,2));




        return forward_diff_xyz; 
	}

   
    //analytical forces

    double lennard_jones_analytical_forces(double x_j, double y_j,double z_j, double gold_ep_ij,double gold_sig_ij,double h) {
  
        //x partial deriv
        double x_ij= 2*(x-x_j);
        double r_ij_x= 1/2*sqrt(x_ij);
        double lj_energy_x = gold_ep_ij*((pow(gold_sig_ij,12)/pow(r_ij_x,13)) - 12*(pow(gold_sig_ij,6)/pow(r_ij_x,7)));

        //y partial deriv
        double y_ij= 2*(y-y_j);
        double r_ij_y= 1/2*sqrt(y_ij);
        double lj_energy_y = gold_ep_ij*((pow(gold_sig_ij,12)/pow(r_ij_y,13)) - 12*(pow(gold_sig_ij,6)/pow(r_ij_y,7)));

        //z partial deriv
        double z_ij= 2*(z-z_j);
        double r_ij_z= 1/2*sqrt(z_ij);
        double lj_energy_z = gold_ep_ij*((pow(gold_sig_ij,12)/pow(r_ij_z,13)) - 12*(pow(gold_sig_ij,6)/pow(r_ij_z,7)));


        double lj_energy_xyz=sqrt(pow(lj_energy_x,2)+pow(lj_energy_y,2)+pow(lj_energy_z,2));


        return lj_energy_xyz;


	}


    double steepest_descent_central(double starting_point_x, double starting_point_y,double starting_point_z,double stepsize, int total)
    {
        double gold_ep_ij=2.951;
        double gold_sig_ij=5.29;
        double deriv=lennard_jones_center_diff(starting_point_x, starting_point_y,starting_point_z,gold_ep_ij,gold_sig_ij,stepsize);
        int count=0;

        while(deriv > total && count < 1e-6)
        {
            double new_point_x = starting_point_x - stepsize * deriv;
            double new_point_y = starting_point_y - stepsize * deriv;
            double new_point_z = starting_point_z - stepsize * deriv;

            if (lennard_jones(new_point_x,new_point_y,new_point_z,gold_ep_ij,gold_sig_ij)<lennard_jones(starting_point_x, starting_point_y,starting_point_z,gold_ep_ij,gold_sig_ij))
            {
                starting_point_x=new_point_x;
                starting_point_y=new_point_y;
                starting_point_z=new_point_z;
                deriv=lennard_jones_center_diff(starting_point_x,starting_point_y,starting_point_z,gold_ep_ij,gold_sig_ij,stepsize);
                stepsize*=1.2;


            }

            else{
                stepsize/=2;
            
            }

            count+=1;
            return lennard_jones(starting_point_x,starting_point_y,starting_point_z,gold_ep_ij,gold_sig_ij);

        }

        
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

		arma::mat input_1 = {x.steepest_descent_central(0,0,1,0.01,80)};
		input_1.print("center diff steep descent");


	}
	return EXIT_SUCCESS;
}



