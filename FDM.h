#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include <memory>
#include <fstream>
#include <string>
#include <Eigen/Dense>

class FDM
{
private:
	int N; // Hamiltonian & Wavefunction dimensions
	double xmin, xmax; // Max & Min X-Values
	double hbar, m; // Physical constants
	double delta_x; // Our step size
	std::vector <double> x_vals; // Vector to contain our discretized x-axis
	std::vector<std::vector <double>> kinetic; // Our Kinetic matrix
	std::vector<std::vector <double>> potential; // Method that calculates the Potential segment of our Hamiltonian
	std::vector<std::vector<double>> eigenvectors; // We do need somewhere to store the entire point of this project to be fair 
	std::vector<std::vector<double>> perturbed_potential; //We need this to get H = T + (V_NP + V_P)
public:
	//Destructor and Constructor with default values
	FDM() : hbar(1), m(1), N(200), xmin(-5), xmax(5), delta_x(0.001)
	{
		int N;
		std::cout << "Enter the number of intervals required :" << std::endl;
		std::cin >> N;
		this->N = N;
		// Resizing our kinetic, potential and hamiltonian 'matrices'
		kinetic.resize(N, std::vector<double>(N, 0));
		potential.resize(N, std::vector<double>(N, 0));
		hamiltonian.resize(N, std::vector<double>(N, 0));
		perturbed_potential.resize(N, std::vector<double>(N, 0));
	};
	//Destructor and Constructor with default values
	FDM(int n) : hbar(1), m(1), N(200), xmin(-4), xmax(4), delta_x(0.001)
	{
		this->N = n;
		// Resizing our kinetic, potential and hamiltonian 'matrices'
		kinetic.resize(N, std::vector<double>(N, 0));
		potential.resize(N, std::vector<double>(N, 0));
		hamiltonian.resize(N, std::vector<double>(N, 0));
		perturbed_potential.resize(N, std::vector<double>(N, 0));
	};
	~FDM() {};
	// A few getter/setter methods
	void set_h();
	void set_m(); 
	void set_lims();
	double get_m();
	double get_h();
	std::vector<std::vector <double>> hamiltonian; // Our Hamiltonian, which we need the Eigenvalues of
	double double_well_potential(double a, double b,double c, double x); // A function that implements the double well Ax**4 + Bx**2 for a 1D case
	double harmonic_oscillator_potential(double m, double omega, double x);// This simply implements the QHO Potential with default parameters
	double perturbation(double lambda, double x); //This allows the perturbation calculation
	void discretize_x(); // Discretizing our x-axis for the FDM
	void calculate_k(); // Method that calculates the Kinetic segment of our Hamiltonian
	void calculate_v(); // Method that calculates the Potential segment of our Hamiltonian
	void calculate_hamiltonian(); // Method that creates our Hamiltonian
	Eigen::MatrixXd vector_to_eigen_matrix(const std::vector<std::vector <double>> &vec); // A helper method to convert our vectors to the Eigen library Dynamic Matrix type
	void calculateEigenValuesAndVectors(const Eigen::MatrixXd& matrix, const std::string& filename); //Method to write the first 10 Eigenvectots to a text file using the Eigen library
	bool isSymmetric(const Eigen::MatrixXd& matrix); // This was a helper method, to make sure our Hamiltonian was symmetric, to allow us to use the SelfAdjoint Eigensolver function from the Eigen Library
	void new_potential(double lambda); // This method will be our Perturbation check for the QHO
	void eigen_info(const Eigen::MatrixXd& matrix, int n, const std::string & filename); //Slightly tailored Eigensolver, to only get th first 2 states, and write to a subdirectory
};

