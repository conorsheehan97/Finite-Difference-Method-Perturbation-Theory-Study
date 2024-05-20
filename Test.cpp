#include "FDM.h"

using namespace std; 

int main()
{
	// Using Scoping and Smarrt Pointers to Assign/Delete memory at runtime for efficiency
	/* {
		// First off the Double Well implementation 
		std::unique_ptr<FDM> obj = std::make_unique<FDM>();
		obj->discretize_x();
		obj->calculate_k();
		obj->calculate_v();
		obj->calculate_hamiltonian();
		obj->calculateEigenValuesAndVectors(obj->vector_to_eigen_matrix(obj->hamiltonian), "Double_Well_Potential.txt");
	}*/
	/* {
		// Secondly the Harmonic Oscillator Potential
		std::unique_ptr<FDM> obj_qho = std::make_unique<FDM>();
		obj_qho->discretize_x();
		obj_qho->calculate_k();
		obj_qho->calculate_v();
		obj_qho->calculate_hamiltonian();
		obj_qho->calculateEigenValuesAndVectors(obj_qho->vector_to_eigen_matrix(obj_qho->hamiltonian), "Quantum_Harmonic_Oscillator.txt");
	}*/
	{
		// Now for the QHO Perturbation Example
		for (int i = 0; i < 1000; ++i)
		{
			double lambda = 0.001 + i * .001;
			std::cout << "Presently at lambda = " << lambda << "\n";
			std::string filename = "Perturbation/" + to_string(lambda) + ".txt";
			std::unique_ptr<FDM> obj_per = std::make_unique<FDM>(200);
			obj_per->discretize_x();
			obj_per->calculate_k();
			obj_per->new_potential(lambda);
			obj_per->calculate_hamiltonian();
			obj_per->eigen_info(obj_per->vector_to_eigen_matrix(obj_per->hamiltonian), 2, filename);
		}
		
	}
	return 0;
}
