#include "FDM.h"

void FDM::set_h()
{
	// This allows us to either set hbar to its default value of 1, or the literature value
	char choice;
	do {
		std::cout << "\n-------------------------------------------------------------------";
		std::cout << "\nWould you like to set hbar to 1, or use the literature value?\n";
		std::cout << "Press 'l' for literature value, or 'd' for the default value of 1 \n";
		std::cin >> choice;
		switch (choice)
		{
		case 'l':
		case 'L': 
			std::cout << "\nYou chose the literature value, how educated of you\n"; 
			std::cout << "\n-------------------------------------------------------------------";
			hbar = 1.054E-34; break;
		case 'd':
		case 'D':
			std::cout << "\n You chose the default value. How uneducated of you\n";
			std::cout << "\n-------------------------------------------------------------------";
			hbar = 1; break;
		}
	} while (choice != 'l' && choice != 'L' && choice != 'd' && choice != 'D');
}
void FDM::set_m()
{
	// Nearly identical to the other method, this just allows for integer calues of m. No need for doubles here in my opinion
	// Also doubles wouldn't let us use a switch statement
	std::cout << "\n-------------------------------------------------------------------";
	bool condition = false;
	int m;
	do {
		std::cout << "\nPlease enter a numeric value for m. I recommend 1\n";
		std::cin >> m;
		switch(m)
		{
		case 1: 
			std::cout << "\nYou played it safe. Good idea\n";
			std::cout << "\n-------------------------------------------------------------------";
			m = m; break;
			condition = true;
		default:
			std::cout << "\nYou crazy catfish, fair enough\n";
			std::cout << "\n-------------------------------------------------------------------";
			condition = true;
		}
	} while (condition = false);
}
void FDM::set_lims()
{
	// This is a setter method for x-axis limits. Simple enough, but includes a quick check of the value sizes
	double a, b;
	std::cout << "\n-------------------------------------------------------------------\n";
	std::cout << "Please enter values for X-MIN and X_MAX\n"; 
	std::cin >> a;
	std::cin >> b;
	if (b < a)
	{
		double temp = b;
		b = a; 
		a = temp;
		xmin = a;
		xmax = b;
		std::cout << "\n-------------------------------------------------------------------\n";
	}
	else
	{
		xmin = a;
		xmax = b;
		std::cout << "\n-------------------------------------------------------------------\n";
	}
}
double FDM::get_h()
{
	return hbar;
}
double FDM::get_m()
{
	return m;
}
void FDM::discretize_x()
{
	double x_n = xmin;
	delta_x = (xmax - xmin) / (N-1);
	for (int i = 0; i < N; i++)
	{
		x_n = xmin + i * delta_x;
		x_vals.push_back(x_n);
	}
}
void FDM::calculate_k()
{
	double prefix = -1*(hbar * hbar) / (2 * m * pow(this->delta_x,2));
	for (int i = 0; i < kinetic.size(); ++i)
	{
		for (int j = 0; j < kinetic.size(); ++j)
		{
			if (i == j)
			{
				kinetic[i][j] = -2*prefix;

			}
			else if (abs(i - j) == 1)
			{
				kinetic[i][j] = 1*prefix;
			}
			else
			{
				kinetic[i][j] = 0;
			}
		}
		
	}
}
double FDM::double_well_potential(double a, double b, double c, double x)
{
	return (a * pow(x, 4) - b * pow(x, 2) + c);
}
double FDM::harmonic_oscillator_potential(double m, double omega, double x)
{
	return (0.5 * m * omega * omega * x * x);
}
double FDM::perturbation(double lambda, double x)
{
	return (lambda * pow(x, 4));
}
void FDM::calculate_v() 
{
	// This method allows the user to pick a potential to investigate, and calculate eigenvalues for. 
	std::cout << "\n-------------------------------------------------------------------";
	bool condition = false;
	char potential_choice;
	do
	{
		std::cout << "\nPlease choose a Potential Energy to investigate:\n";
		std::cout << "Double Well Potential - Please enter 'd'\n";
		std::cout << "Harmonic Oscillator - Please enter 'h'\n";
		std::cin >> potential_choice;
		switch (potential_choice)
		{
		case 'd':
		case 'D':
			for (int i = 0; i < potential.size(); ++i)
			{
				for (int j = 0; j < potential.size(); ++j)
				{
					if (i == j)
					{
						potential[i][j] = double_well_potential(1, 3, 0, x_vals[i]);
					}
					else
					{
						potential[i][j] = 0;
					}
				}
			}

			std::cout << "\n-------------------------------------------------------------------\n";
			condition = true;
			break;
		case 'h':
		case 'H':
			for (int i = 0; i < potential.size(); ++i)
			{
				for (int j = 0; j < potential.size(); ++j)
				{
					if (i == j)
					{
						potential[i][j] = harmonic_oscillator_potential(1,1,x_vals[i]);
					}
					else
					{
						potential[i][j] = 0;
					}
				}
			}
			std::cout << "\n-------------------------------------------------------------------\n";
			condition = true;
			break;
		}


		} while (condition = false);
	}
void FDM::calculate_hamiltonian()
{
	// This simply adds the Kinetic and Potential Matrices to give the Hamiltonian Matrix
	for (int i = 0; i < hamiltonian.size(); ++i)
	{
		for (int j = 0; j < hamiltonian.size(); ++j)
		{
			hamiltonian[i][j] = kinetic[i][j] + potential[i][j];
		}
	}
}
Eigen::MatrixXd FDM::vector_to_eigen_matrix(const std::vector<std::vector<double>> &vec)
{
	//This method allows us to convert our vector to the appropriate data type for use with the Eigen library 
	int rows = vec.size();
	int cols = vec[0].size();
	Eigen::MatrixXd matrix(rows, cols);

	for (int i = 0; i < rows; ++i) 
	{
		for (int j = 0; j < cols; ++j) 
		{
				matrix(i, j) = vec[i][j];
		}
		}

	return matrix;
}
void FDM::calculateEigenValuesAndVectors(const Eigen::MatrixXd& matrix, const std::string& filename) {
	// Check if the matrix is square
	if (matrix.rows() != matrix.cols()) {
		std::cerr << "Matrix must be square!" << std::endl;
		return;
	}

	// Calculate eigenvalues and eigenvectors
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);

	// Extract the eigenvalues and eigenvectors
	Eigen::VectorXd eigenvalues = solver.eigenvalues().real();
	Eigen::MatrixXd eigenvectors = solver.eigenvectors().real();

	// Open a file to write the results
	std::ofstream outFile(filename);
	if (!outFile.is_open()) 
	{
		std::cerr << "Unable to open file: " << filename << std::endl;
		return;
	}
	int numEigenvalues = std::min(5, static_cast<int>(eigenvalues.size()));

	// Write eigenvalues
	outFile << "Eigenvalues:\n";
	for (int i = 0; i < numEigenvalues; ++i) {
		outFile << eigenvalues[i] << "\t";
	}

	// Write eigenvectors
	outFile << "\nEigenvectors:\n";
	for (int i = 0; i < eigenvectors.rows(); ++i)
	{
		outFile << x_vals[i] << "\t";
		for (int j = 0; j < numEigenvalues; ++j) 
		{
			outFile << eigenvectors(i, j) << "\t";
		}
		outFile << "\n";
	}

	outFile.close();
	std::cout << "First 5 eigenvalues and eigenvectors written to " << filename << std::endl;
}
bool FDM::isSymmetric(const Eigen::MatrixXd& matrix)
{
	// Compare the matrix with its transpose
	return matrix.isApprox(matrix.transpose());
}
void FDM::new_potential(double lambda)
{
	//First we need to construct our new Potential for our Hamiltonian. We'll assume V_NP is already set
	for (int i = 0; i < potential.size(); ++i)
	{
		for (int j = 0; j < potential.size(); ++j)
		{
			if (i == j)
			{
				potential[i][j] = harmonic_oscillator_potential(1,1,x_vals[i]) + perturbation(lambda, x_vals[i]);
			}
			else
			{
				potential[i][j] = 0;
			}
		}
	}
};
void FDM::eigen_info(const Eigen::MatrixXd &matrix, int n, const std::string& filename)
{
	// Check if the matrix is square
	if (matrix.rows() != matrix.cols()) {
		std::cerr << "Matrix must be square!" << std::endl;
		return;
	}

	// Calculate eigenvalues and eigenvectors
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);

	// Extract the eigenvalues and eigenvectors
	Eigen::VectorXd eigenvalues = solver.eigenvalues().real();
	Eigen::MatrixXd eigenvectors = solver.eigenvectors().real();

	// Open a file to write the results
	std::ofstream outFile(filename);
	if (!outFile.is_open())
	{
		std::cerr << "Unable to open file: " << filename << std::endl;
		return;
	}
	int numEigenvalues = n;

	// Write eigenvalues
	outFile << "Eigenvalues:\n";
	for (int i = 0; i < numEigenvalues; ++i) {
		outFile << eigenvalues[i] << "\t";
	}

	// Write eigenvectors
	outFile << "\nEigenvectors:\n";
	for (int i = 0; i < eigenvectors.rows(); ++i)
	{
		outFile << x_vals[i] << "\t";
		for (int j = 0; j < numEigenvalues; ++j)
		{
			outFile << eigenvectors(i, j) << "\t";
		}
		outFile << "\n";
	}

	outFile.close();
}
