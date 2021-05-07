#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>

float rotateX (float angle_radians)
{
	float sineVariable;
	sineVariable = sinf (angle_radians);
	return sineVariable;
}

float rotateZ (float angle_radians)
{
	float cosineVariable;
	cosineVariable = cosf (angle_radians);
	return cosineVariable;
}

float degreeToRadians(float angle_degrees)
{
	float angle_radians;
	angle_radians = angle_degrees/57.2958;
	return angle_radians;
}

struct information
{
	// Atom ID is stored
	int CH3, CH2, CH, chain_id_number;
};

int main(int argc, char const *argv[])
{
	if(argc > 1)
	{
		// Processing argc and argv
		if(strstr(argv[1], "--help"))
		{
			printf("\n[*] Input:\n    ~~~~~~\n\n\t\"dump_equilibrate600.melt\"\n\n[*] Output:\n    ~~~~~~~\n\n\t\"NOT YET IMPLEMENTED\"\n\n[*] Description:\n    ~~~~~~~~~~~~\n\nThis program calculates order parameters from input dump file. This version works only for single chain isotactic polypropylene. Some of the order parameters calculated were\n\n\t[*] Local order paramter (vs Distance, vs time),\n\t[*] Intermolecular order parameter (vs radial distance),\n\t[*] Chain orientation order parameter (along X, Y and Z axes; vs time)\n\nChain orientation order parameter can be calculated based on any custom vector as a basis and it's not restricted to X/Y/Z axes.\n\n");
			exit(1);
		}
		else
		{
			printf("\nError:\n~~~~~~\n\n\tUnknown command passed. Type '--help' for more information or run the executable directly\n\n");
			exit(1);
		}
	}

	// Input file
	FILE *read;
	read = fopen("../dump_equilibrate600.melt", "r");

	// Output file
	FILE *LOPdump, *avgLOP, *COOP, *IOP;
	LOPdump = fopen("dump_lop.csv", "w");
	avgLOP = fopen("avg_lop.csv", "w");
	COOP = fopen("coop.csv", "w");
	IOP = fopen("iop.csv", "w");


	// Adding a header to *avgLOP file
	fprintf(avgLOP, "timeframe,average local_order_parameter\n");

	// To count the scanning timeframe
	int currentTimeframe = 0;

	// currentLine = to scan through and identify line by number from dump file
	// natoms = to store total number of atoms, from dump file
	// nmobileatoms = to store the total number of mobile atoms (aka UAs in chains)
	int currentLine = 0, natoms, nmobileatoms = 0, isFirstTimeframe = 1;
	char lineString[3000];

	struct information *monomer;

	// Declaring variables for loops
	int a;

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// ~~~~~~~~~~~~~~~~~~~ VARIABLE DECLARATIONS ~~~~~~~~~~~~~~~~~~~
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// ~~~~~~~~~~~~~~~~~~~ VARIABLES TO STORE DUMP INPUT ~~~~~~~~~~~~~~~~~~~
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		int *id_dump, *molType_dump, *atomType_dump, *ix_dump, *iy_dump, *iz_dump, currentAtom = 0, currentMonomer = 0, currentChain = 1;
		float *x_dump, *y_dump, *z_dump, xlo, xhi, ylo, yhi, zlo, zhi;

		// ~~~~~~~~~~~~~~~~~~~ VARIABLES TO COMPUTE LOCAL ORDER PARAMETER ~~~~~~~~~~~~~~~~~~~
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		float x1, y1, z1, x2, y2, z2, vx1, vx2, vy1, vy2, vz1, vz2, dotproduct, magnitude1, magnitude2, theta, costheta, *LOCAL_ORDER_PARAMETER, AVERAGE_LOCAL_ORDER_PARAMETER, sum_local_order_parameter = 0, *CENTER_OF_LOCAL_ORDER_PARAMETER_X, *CENTER_OF_LOCAL_ORDER_PARAMETER_Y;
		int i, j, k, l, NUMBER_OF_LOPs = 0, degree_of_polymerization, *LOP_ATOM1, *LOP_ATOM2, *LOP_ATOM3, *LOP_ATOM4; // Four atoms involved in calculating vectors and local order parameter

		// ~~~~~~~~~~~~~~~~~~~ VARIABLES TO COMPUTE LOP AS A FUNCTION OF Y-AXIS ~~~~~~~~~~~~~~~~~~~
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		int nBins;
		float y_lower_end, y_higher_end, y_distance, y_bin_length, *y_bin_boundary, *Y_vs_LOCAL_ORDER_PARAMETER_AVERAGE, y_local_order_parameter_sum = 0, y_local_order_parameter_N = 0;
		// To compute LOP as a function of distance from nucleating agent.
		nBins = 2;
		// printf("How many bins are required for computing LOP vs Y-axis\n");
		// scanf("%d", &nBins);
		// Splitting the simulation box into various bins
		y_bin_boundary = (float *) calloc ((nBins + 1), sizeof (float)); // To check if CoM is inbetween the bin boundaries
		Y_vs_LOCAL_ORDER_PARAMETER_AVERAGE = (float *) calloc (nBins, sizeof (float)); // To plot LOP as a function of Y-axis


		// ~~~~~~~~~~~~~~~~~~~ VARIABLES TO COMPUTE INTERMOLECULAR ORDER PARAMETER ~~~~~~~~~~~~~~~~~~~
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		int NUMBER_OF_IOPs = 0, index1_iop, index2_iop, index3_iop, index4_iop, nbins_iop, n_iop_within_bin, NUMBER_OF_IOPs_actual = 0;
		float *x1_iop, *y1_iop, *z1_iop, *x2_iop, *y2_iop, *z2_iop, *x_com_iop, *y_com_iop, *z_com_iop, distance_iop;
		int r, dr; // 'r' starts from 0, the center. 'dr' is user input variable; r and dr are int because float cannot be used as array index.
		int *denominator_iop, iop_array_variable = 0;
		float INTERMOLECULAR_ORDER_PARAMETER, *sum_iop, *AVERAGE_INTERMOLECULAR_ORDER_PARAMETER, *SD_INTERMOLECULAR_ORDER_PARAMETER; // Intermolecular order parameter is the order parameter between pairs of vectors; average intermolecular order parameter is the average of all IOPs within a radial shell; standard deviation of intermolacular order parameter is the SD corresponding to that average; avg and sd are arrays, storing values for radial shells in a particular timeframe.

		// Get 'dr' from user
		dr = 1;
		// printf ("Enter dr value (for intermolecular order parameter): \n");
		// scanf ("%d", &dr);

		// ~~~~~~~~~~~~~~~~~~~ VARIABLES TO COMPUTE CHAIN-CRYSTAL ORDER PARAMETER ~~~~~~~~~~~~~~~~~~~
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		int NUMBER_OF_COOPs;
		float angle_degrees, angle_radians;
		float X1_COOP, X2_COOP, Y1_COOP, Y2_COOP, Z1_COOP, Z2_COOP;
		// First array var is for angle and second array var is for IOPs scanned
		// Example: CHAIN_ORIENTATION_ORDER_PARAMETER[angle_radians][NUMBER_OF_IOPs]
		float **CHAIN_ORIENTATION_ORDER_PARAMETER;
		float *CHAIN_ORIENTATION_ORDER_PARAMETER_AVERAGE, *CHAIN_ORIENTATION_ORDER_PARAMETER_SUM;
		int CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT = 0;

		// Memory allocation for CHAIN_ORIENTATION_ORDER_PARAMETER
		CHAIN_ORIENTATION_ORDER_PARAMETER_AVERAGE = (float *) calloc (360, sizeof(float));
		CHAIN_ORIENTATION_ORDER_PARAMETER_SUM = (float *) calloc (360, sizeof(float));

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// ~~~~~~~~~~~~~~~~~~~ VARIABLE DECLARATIONS END ~~~~~~~~~~~~~~~~~~~
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	while(fgets(lineString, 3000, read) != NULL)
	{
		currentLine++;

		// Get the total number of atoms and declare all variables
		if (currentLine == 4 && isFirstTimeframe == 1)
		{
			sscanf(lineString, "%d", &natoms);

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// MEM ALLOCATION - DUMP PARAMETERS
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Maximum number of chains possible: (natoms / 3); since there are atleast 3 UAs in every chain. If we account for end groups, it'll be even fewer.
			monomer = (struct information *) calloc ((natoms / 3), sizeof (struct information));
			id_dump = (int *) malloc (natoms * sizeof (int));
			molType_dump = (int *) malloc (natoms * sizeof (int));
			atomType_dump = (int *) malloc (natoms * sizeof (int));
			ix_dump = (int *) malloc (natoms * sizeof (int));
			iy_dump = (int *) malloc (natoms * sizeof (int));
			iz_dump = (int *) malloc (natoms * sizeof (int));
			x_dump = (float *) malloc (natoms * sizeof (float));
			y_dump = (float *) malloc (natoms * sizeof (float));
			z_dump = (float *) malloc (natoms * sizeof (float));

			degree_of_polymerization = (natoms / 3);

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// MEM ALLOCATION - LOCAL ORDER PARAMETER
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// This is gross over-estimate of total number of LOPs present.
			// Used for allocating memory
			NUMBER_OF_LOPs = (degree_of_polymerization * 2);
			LOCAL_ORDER_PARAMETER = (float *) calloc (NUMBER_OF_LOPs, sizeof (float));
			LOP_ATOM1 = (int *) calloc (NUMBER_OF_LOPs, sizeof (int));
			LOP_ATOM2 = (int *) calloc (NUMBER_OF_LOPs, sizeof (int));
			LOP_ATOM3 = (int *) calloc (NUMBER_OF_LOPs, sizeof (int));
			LOP_ATOM4 = (int *) calloc (NUMBER_OF_LOPs, sizeof (int));
			CENTER_OF_LOCAL_ORDER_PARAMETER_X = (float *) calloc (NUMBER_OF_LOPs, sizeof (float));
			CENTER_OF_LOCAL_ORDER_PARAMETER_Y = (float *) calloc (NUMBER_OF_LOPs, sizeof (float));

			// Atom type 1 signifies backbone atoms and 0 signifies pendant groups/end groups.
			// By using calloc, all elements are auto assigned to zero.
			// In switch statement, first UA, which is the end group of the chain, is ignored, which should be zero.
			// Using Calloc means, it's zero by default.
			// atomType = (int *) calloc (natoms, sizeof(float));

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// MEM ALLOCATION - INTERMOLECULAR ORDER PARAMETER
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			NUMBER_OF_IOPs = degree_of_polymerization + 1;
			x1_iop = (float *) calloc (NUMBER_OF_IOPs, sizeof(float));
			y1_iop = (float *) calloc (NUMBER_OF_IOPs, sizeof(float));
			z1_iop = (float *) calloc (NUMBER_OF_IOPs, sizeof(float));
			x2_iop = (float *) calloc (NUMBER_OF_IOPs, sizeof(float));
			y2_iop = (float *) calloc (NUMBER_OF_IOPs, sizeof(float));
			z2_iop = (float *) calloc (NUMBER_OF_IOPs, sizeof(float));
			x_com_iop = (float *) calloc (NUMBER_OF_IOPs, sizeof(float));
			y_com_iop = (float *) calloc (NUMBER_OF_IOPs, sizeof(float));
			z_com_iop = (float *) calloc (NUMBER_OF_IOPs, sizeof(float));

			// Intermolecular order parameter, variable allocation
			// Number of bins in intermolecular order parameter calculations (radial bins)
			// Change this number based on system size and 'dr' distance
			// nbins_iop is directly proportional to system size and inversely proportional to 'dr' distance
			nbins_iop = 200;
			AVERAGE_INTERMOLECULAR_ORDER_PARAMETER = (float *) calloc (nbins_iop, sizeof(float));
			SD_INTERMOLECULAR_ORDER_PARAMETER = (float *) calloc (nbins_iop, sizeof(float));
			sum_iop = (float *) calloc (nbins_iop, sizeof(float));
			denominator_iop = (int *) calloc (nbins_iop, sizeof(float));

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// MEM ALLOCATION - CHAIN ORIENTATION ORDER PARAMETER
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			NUMBER_OF_COOPs = (degree_of_polymerization - 3) * 2;
			// CHAIN_ORIENTATION_ORDER_PARAMETER_X = (float *) calloc (NUMBER_OF_COOPs, sizeof(float));
			// CHAIN_ORIENTATION_ORDER_PARAMETER_Y = (float *) calloc (NUMBER_OF_COOPs, sizeof(float));
			// CHAIN_ORIENTATION_ORDER_PARAMETER_Z = (float *) calloc (NUMBER_OF_COOPs, sizeof(float));

			CHAIN_ORIENTATION_ORDER_PARAMETER = (float **) calloc (720, sizeof(float));
			for (int arrayVar = 0; arrayVar < 720; ++arrayVar)
			{
				CHAIN_ORIENTATION_ORDER_PARAMETER[arrayVar] = (float *) calloc (NUMBER_OF_IOPs, sizeof(float));
			}
		}

		// Gather information about simulation box dimension
		// sim box dimension along x-axis
		if (currentLine == 6 && isFirstTimeframe == 1)
		{
			sscanf(lineString, "%f %f", &xlo, &xhi);
			// printf("Box dimensions:\n\n\txlo: %.3f; xhi: %.3f;\n", xlo, xhi);
		}
		// sim box dimension along y-axis
		if (currentLine == 7 && isFirstTimeframe == 1)
		{
			sscanf(lineString, "%f %f", &ylo, &yhi);
			// printf("\tylo: %.3f; yhi: %.3f;\n", ylo, yhi);
		}
		// sim box dimension along z-axis
		if (currentLine == 8 && isFirstTimeframe == 1)
		{
			sscanf(lineString, "%f %f", &zlo, &zhi);
			// printf("\tzlo: %.3f; zhi: %.3f;\n", zlo, zhi);
		}

		// Block to read atom entries
		// The atom coordinates are then corrected based on image and then unwraped
		// Coordinates are stored in *x_dump, *y_dump, *z_dump
		if (currentLine >= 10 && currentLine <= (natoms + 9))
		{
			sscanf(lineString, "%d %d %d %f %f %f %*f %*f %*f %d %d %d", 
				&id_dump[currentAtom], 
				&molType_dump[currentAtom], 
				&atomType_dump[currentAtom], 
				&x_dump[currentAtom], 
				&y_dump[currentAtom], 
				&z_dump[currentAtom], 
				&ix_dump[currentAtom], 
				&iy_dump[currentAtom], 
				&iz_dump[currentAtom]);

			if (molType_dump[currentAtom] == 1)
			{
				nmobileatoms++;
			}

			// checking the uncorrected coordinates
			// if (ix_dump[currentAtom] != 0 || iy_dump[currentAtom] != 0 || iz_dump[currentAtom] != 0)
			// {
			// 	printf("original coords: x: %.3f; y: %.3f; z: %.3f;\n", x_dump[currentAtom], y_dump[currentAtom], z_dump[currentAtom]);
			// }

			// Coordinate correction based on ix, iy and iz
			if (ix_dump[currentAtom] > 0)
			{
				x_dump[currentAtom] = xhi + ((abs(ix_dump[currentAtom]) - 1) * (xhi - xlo)) + (x_dump[currentAtom] - xlo);
			}
			else if (ix_dump[currentAtom] < 0)
			{
				x_dump[currentAtom] = xlo - ((abs(ix_dump[currentAtom]) - 1) * (xhi - xlo)) - (xhi - x_dump[currentAtom]);
			}
			if (iy_dump[currentAtom] > 0)
			{
				x_dump[currentAtom] = xhi + ((abs(ix_dump[currentAtom]) - 1) * (xhi - xlo)) + (x_dump[currentAtom] - xlo);
			}
			else if (iy_dump[currentAtom] < 0)
			{
				x_dump[currentAtom] = xlo - ((abs(ix_dump[currentAtom]) - 1) * (xhi - xlo)) - (xhi - x_dump[currentAtom]);
			}
			if (iz_dump[currentAtom] > 0)
			{
				x_dump[currentAtom] = xhi + ((abs(ix_dump[currentAtom]) - 1) * (xhi - xlo)) + (x_dump[currentAtom] - xlo);
			}
			else if (iz_dump[currentAtom] < 0)
			{
				x_dump[currentAtom] = xlo - ((abs(ix_dump[currentAtom]) - 1) * (xhi - xlo)) - (xhi - x_dump[currentAtom]);
			}
			currentAtom++;

			// Checking the corrected coordinates
			// if (ix_dump[currentAtom] != 0 || iy_dump[currentAtom] != 0 || iz_dump[currentAtom] != 0)
			// {
			// 	printf("corrected coords: x: %.3f; y: %.3f; z: %.3f;\n", x_dump[currentAtom], y_dump[currentAtom], z_dump[currentAtom]);
			// 	printf("image information: ix_dump: %d; iy_dump: %d; iz_dump: %d;\n", ix_dump[currentAtom], iy_dump[currentAtom], iz_dump[currentAtom]);
			// 	sleep(2);
			// }

			// The repeating units along the chain length based on atom type: [1, 3, 2]
		} 

		// The code block is executed at the end of every timeframe
		// Under this block, all calculations are to be implemented
		if (currentLine == (natoms + 9))
		{
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// BLOCK TO IDENTIFY UNITED ATOMS
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Re-allocating variables at the end of every scan (and before starting the computations)
			currentLine = 0;
			currentAtom = 0;
			currentMonomer = 0;
			currentChain = 1;
			isFirstTimeframe = 0;
			iop_array_variable = 0;

			// Cycle through all the stored atom entries to identify them.
			for (i = 0; i < natoms; ++i)
			{
				// Checks if the mol_type is equal to 1, which corresponds to chain. Mol_type == 2 corresponds to crystal substrate
				if ((i + 2) < natoms && molType_dump[i] == 1)
				{
					// Checks if the current atom is a part of a monomer
					if (atomType_dump[i] == 1 && atomType_dump[i + 1] == 3 && atomType_dump[i + 2] == 2)
					{
						monomer[currentMonomer].CH = id_dump[i];
						monomer[currentMonomer].CH3 = id_dump[i + 1];
						monomer[currentMonomer].CH2 = id_dump[i + 2];
						monomer[currentMonomer].chain_id_number = currentChain;

						// Checking the above struct array
						if (currentChain > 0)
						{
							// printf("CH: %d; CH2: %d; CH3: %d; currentChain: %d;\n", 
							// 	monomer[currentMonomer].CH,
							// 	monomer[currentMonomer].CH2,
							// 	monomer[currentMonomer].CH3,
							// 	monomer[currentMonomer].chain_id_number);
						}

						i = i + 2; // Skip the other two atoms, since they are already taken into account above
						currentMonomer++;
					}
					// Checks for end of chain
					// If 'i' corresponds to chain end, then increment 'currentChain'
					else if (atomType_dump[i] == 1 && atomType_dump[i + 1] == 3 && atomType_dump[i + 2] == 3)
					{
						currentChain++;
					}
				}
			}
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// BLOCK TO IDENTIFY UNITED ATOMS ENDS
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// LOCAL ORDER PARAMETER CALCULATIONS - STARTS
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			// Reset 'NUMBER_OF_LOPs' and store the actual number of LOPs present.
			// This value can be used later to cycle through all the available LOP values.
			NUMBER_OF_LOPs = 0;	
			// for-loop-1 for 1st combination of local order parameter
			// Cycle through all monomers (currentMonomer was last used to increment struct array)
			for (a = 0; a < currentMonomer; a++)
			{
				// i = CH3(currentmonomer); j = CH3(nextmonomer);
				// k = CH2(currentmonomer); l = CH2(previousmonomer)
				// k = CH2(nextmonomer); l = CH2(nextmonomer+1)
				if (a >= 3 && (a+3) < currentMonomer)
				{
					i = monomer[a].CH3;
					j = monomer[a+3].CH3;
					k = monomer[a].CH2;
					l = monomer[a-3].CH2;

					// Checking if all monomers under consideration are from the same chain
					if ((i <= natoms) && (j <= natoms) && (k <= natoms) && (k > 1) && (l <= natoms) && (l > 1) && monomer[a].chain_id_number == monomer[a+3].chain_id_number && monomer[a].chain_id_number == monomer[a-3].chain_id_number)
					{
						// printf("loop-1 ==> i: %d; j: %d; k: %d; l: %d\n", i, j, k, l);	
						// printf("Atom1 = %d; Atom2 = %d; <===> Atom3 = %d; Atom4 = %d\n", i, j, k, l);
						// printf("Atom1 = (%.3f, %.3f, %.3f); Atom2 = (%.3f, %.3f, %.3f); <===> Atom3 = (%.3f, %.3f, %.3f); Atom4 = (%.3f, %.3f, %.3f)\n\n", x[i-1], y[i-1], z[i-1], x[j-1], y[j-1], z[j-1], x[k-1], y[k-1], z[k-1], x[l-1], y[l-1], z[l-1]);

						// Local order parameter calculations
						vx1 = x_dump[i-1] - x_dump[j-1];
						vy1 = y_dump[i-1] - y_dump[j-1];
						vz1 = z_dump[i-1] - z_dump[j-1];
						vx2 = x_dump[k-1] - x_dump[l-1];
						vy2 = y_dump[k-1] - y_dump[l-1];
						vz2 = z_dump[k-1] - z_dump[l-1];

						dotproduct = (vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2);
						magnitude1 = (vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1);
						magnitude2 = (vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2);
						costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
						theta = acosf (costheta);
						LOCAL_ORDER_PARAMETER[NUMBER_OF_LOPs] = ((3 * costheta * costheta) - 1) / 2;
						LOP_ATOM1[NUMBER_OF_LOPs] = i;
						LOP_ATOM2[NUMBER_OF_LOPs] = j;
						LOP_ATOM3[NUMBER_OF_LOPs] = k;
						LOP_ATOM4[NUMBER_OF_LOPs] = l;
						CENTER_OF_LOCAL_ORDER_PARAMETER_X[NUMBER_OF_LOPs] = (x_dump[i-1] + x_dump[j-1] + x_dump[k-1] + x_dump[l-1]) / 4;
						CENTER_OF_LOCAL_ORDER_PARAMETER_Y[NUMBER_OF_LOPs] = (y_dump[i-1] + y_dump[j-1] + y_dump[k-1] + y_dump[l-1]) / 4;

						// printf("LOCAL_ORDER_PARAMETER: %f\n", LOCAL_ORDER_PARAMETER[NUMBER_OF_LOPs]);
						NUMBER_OF_LOPs++;
					}
				}
				if ((a+5) < currentMonomer)
				{
					i = monomer[a].CH3;
					j = monomer[a+3].CH3;
					k = monomer[a+2].CH2;
					l = monomer[a+5].CH2;

					// Checking if all monomers under consideration are from the same chain
					if ((i <= natoms) && (j <= natoms) && (k <= natoms) && (k > 1) && (l <= natoms) && (l > 1) && monomer[a].chain_id_number == monomer[a+3].chain_id_number && monomer[a].chain_id_number == monomer[a+5].chain_id_number && monomer[a].chain_id_number == monomer[a+2].chain_id_number)
					{
						// printf("loop-2 ==> i: %d; j: %d; k: %d; l: %d\n", i, j, k, l);
						// printf("Atom1 = %d; Atom2 = %d; <===> Atom3 = %d; Atom4 = %d\n", i, j, k, l);
						// printf("Atom1 = (%.3f, %.3f, %.3f); Atom2 = (%.3f, %.3f, %.3f); <===> Atom3 = (%.3f, %.3f, %.3f); Atom4 = (%.3f, %.3f, %.3f)\n\n", x[i-1], y[i-1], z[i-1], x[j-1], y[j-1], z[j-1], x[k-1], y[k-1], z[k-1], x[l-1], y[l-1], z[l-1]);

						// Local order parameter calculations
						vx1 = x_dump[i-1] - x_dump[j-1];
						vy1 = y_dump[i-1] - y_dump[j-1];
						vz1 = z_dump[i-1] - z_dump[j-1];
						vx2 = x_dump[k-1] - x_dump[l-1];
						vy2 = y_dump[k-1] - y_dump[l-1];
						vz2 = z_dump[k-1] - z_dump[l-1];

						dotproduct = (vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2);
						magnitude1 = (vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1);
						magnitude2 = (vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2);
						costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
						theta = acosf (costheta);
						LOCAL_ORDER_PARAMETER[NUMBER_OF_LOPs] = ((3 * costheta * costheta) - 1) / 2;
						LOP_ATOM1[NUMBER_OF_LOPs] = i;
						LOP_ATOM2[NUMBER_OF_LOPs] = j;
						LOP_ATOM3[NUMBER_OF_LOPs] = k;
						LOP_ATOM4[NUMBER_OF_LOPs] = l;
						CENTER_OF_LOCAL_ORDER_PARAMETER_X[NUMBER_OF_LOPs] = (x_dump[i-1] + x_dump[j-1] + x_dump[k-1] + x_dump[l-1]) / 4;
						CENTER_OF_LOCAL_ORDER_PARAMETER_Y[NUMBER_OF_LOPs] = (y_dump[i-1] + y_dump[j-1] + y_dump[k-1] + y_dump[l-1]) / 4;

						// printf("LOCAL_ORDER_PARAMETER: %f\n", LOCAL_ORDER_PARAMETER[NUMBER_OF_LOPs]);
						NUMBER_OF_LOPs++;
					}
				}
			}

			// Incrementing and tracking the current timeframe
			currentTimeframe++;

			// Printing current scanning timeframe
			printf("%d\n", currentTimeframe);

			// Printing the local order parameters to output file
			// Header format: (timeframe starts,NUMBER_OF_LOPs: %d)
			fprintf(LOPdump, "timeframe,%d,NUMBER_OF_LOPs:,%d\n", currentTimeframe, NUMBER_OF_LOPs);
			// Output format: (si.no,lop_value,lop_atom1,lop_atom2,lop_atom3,lop_atom4)
			for (i = 0; i < NUMBER_OF_LOPs; ++i)
			{
				fprintf(LOPdump, "%d,%.4f,%d,%d,%d,%d\n", i+1, LOCAL_ORDER_PARAMETER[i], LOP_ATOM1[i], LOP_ATOM2[i], LOP_ATOM3[i], LOP_ATOM4[i]);
			}

			// Calculate average_local_order_parameter for every timeframes
			for (i = 0; i < NUMBER_OF_LOPs; ++i)
			{
				sum_local_order_parameter += LOCAL_ORDER_PARAMETER[i];
			}

			// It is not possible to store average in an array because the total timeframes is not known
			// DON'T FORGET TO SAVE THIS TO A FILE AND THEN USE IT LATER ON
			AVERAGE_LOCAL_ORDER_PARAMETER = (sum_local_order_parameter / NUMBER_OF_LOPs);
			// printf("AVERAGE_LOCAL_ORDER_PARAMETER: %.3f\n", AVERAGE_LOCAL_ORDER_PARAMETER);
			fprintf(avgLOP, "%d\t%f\n", currentTimeframe, AVERAGE_LOCAL_ORDER_PARAMETER);

			// Re-initialize 'sum_local_order_parameter' for next timeframe
			sum_local_order_parameter = 0;

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// LOCAL ORDER PARAMETER CALCULATIONS - ENDS
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// LOCAL ORDER PARAMETER AS A FUNCTION OF Y-AXIS - START
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			y_distance = (yhi - ylo);
			y_bin_length = (y_distance / nBins);
			// printf("y_distance: %.3f\n", y_distance);

			// Assigning y_bin_boundary values
			y_bin_boundary[0] = ylo;
			// printf("y_bin_boundary[0]: %.3f\n", y_bin_boundary[0]);
			// printf("i: %d; y_bin_boundary[0]: %.3f\n", i, y_bin_boundary[0]);
			for (i = 1; i < nBins; ++i)
			{
				y_bin_boundary[i] = y_bin_boundary[i-1] + y_bin_length;
				// printf("y_bin_boundary[%d]: %.3f\n", i, y_bin_boundary[i]);
				// printf("i: %d; y_bin_boundary[%d]: %.3f\n", i, i, y_bin_boundary[i]);
			}
			y_bin_boundary[nBins] = yhi;
			// printf("y_bin_boundary[%d]: %.3f\n", nBins, y_bin_boundary[nBins]);

			// Calculate LOP as a function of Y-axis
			// Loop through all bins
			// printf("\nLOP vs Y-axis:\n");
			for (j = 0; j < nBins; ++j)
			{
				// Initializing sum and denominator, for calculating average
				y_local_order_parameter_sum = 0;
				y_local_order_parameter_N = 0;
				// Loop through all LOP values in this timeframe, then repeated for every other timeframes
				for (i = 0; i < NUMBER_OF_LOPs; ++i)
				{
					// Check if LOP CoM is between the LOP boundaries
					if (CENTER_OF_LOCAL_ORDER_PARAMETER_Y[i] > y_bin_boundary[j] && CENTER_OF_LOCAL_ORDER_PARAMETER_Y[i] < y_bin_boundary[j+1])
					{
						y_local_order_parameter_sum += LOCAL_ORDER_PARAMETER[i];
						y_local_order_parameter_N++;
					}
				}
				Y_vs_LOCAL_ORDER_PARAMETER_AVERAGE[j] = (y_local_order_parameter_sum / y_local_order_parameter_N);
				// printf("\n\tBin-%d:\n\t\tY_vs_LOCAL_ORDER_PARAMETER_AVERAGE[%d]: %f\n", j, j, Y_vs_LOCAL_ORDER_PARAMETER_AVERAGE[j]);
			}

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// LOCAL ORDER PARAMETER AS A FUNCTION OF Y-AXIS - END
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// INTERMOLECULAR ORDER PARAMETER - START
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			// Assigning the center of IOPs
			for (int a = 0; a < currentMonomer; ++a)
			{
				index1_iop = monomer[a].CH2;
				index2_iop = monomer[a].CH;
				index3_iop = monomer[a+3].CH2;
				index4_iop = monomer[a+3].CH;

				// Checking if monomer 'a' and monomer 'a+3' are part of the same chain
				if (index1_iop <= natoms && index2_iop <= natoms && index3_iop <= natoms && index4_iop <= natoms && monomer[a].chain_id_number == monomer[a+3].chain_id_number)
				{
					// Finding the center point of CH---CH2 bond in two monomers
					x1_iop[a] = (x_dump[index1_iop] + x_dump[index2_iop]) / 2;
					y1_iop[a] = (y_dump[index1_iop] + y_dump[index2_iop]) / 2;
					z1_iop[a] = (z_dump[index1_iop] + z_dump[index2_iop]) / 2;
					x2_iop[a] = (x_dump[index3_iop] + x_dump[index4_iop]) / 2;
					y2_iop[a] = (y_dump[index3_iop] + y_dump[index4_iop]) / 2;
					z2_iop[a] = (z_dump[index3_iop] + z_dump[index4_iop]) / 2;
					NUMBER_OF_IOPs_actual++;

					// Center for IOPs (i.e) center of the above calculated two points
					x_com_iop[a] = (x1_iop[a] + x2_iop[a]) / 2;
					y_com_iop[a] = (y1_iop[a] + y2_iop[a]) / 2;
					z_com_iop[a] = (z1_iop[a] + z2_iop[a]) / 2;
				}
			}
			// iop_array_variable = 0;
			// printf("iop_array_variable: %d\n", iop_array_variable);

			for (a = 0; a < NUMBER_OF_IOPs; ++a)
			{
				if (x1_iop[a] + y1_iop[a] + z1_iop[a] + x2_iop[a] + y2_iop[a] + z2_iop[a] != 0) 
				{
					for (r = 0, iop_array_variable = 0; (r+dr) < nbins_iop; iop_array_variable++) // max distance to check = nbins_iop (equals 200 currently. Change it if the system size is different)
					{
						for (j = 0; j < NUMBER_OF_IOPs; ++j) // Loop over other IOP values, calculate the distance
						{
							distance_iop = sqrt(
								((x_com_iop[j] - x_com_iop[i]) * (x_com_iop[j] - x_com_iop[i])) + 
								((y_com_iop[j] - y_com_iop[i]) * (y_com_iop[j] - y_com_iop[i])) + 
								((z_com_iop[j] - z_com_iop[i]) * (z_com_iop[j] - z_com_iop[i]))
								);
							
							if (distance_iop > 0 && distance_iop > r && distance_iop < (r+dr) && (x1_iop[a] + y1_iop[a] + z1_iop[a]) != 0 && (x2_iop[a] + y2_iop[a] + z2_iop[a]) != 0 && (x1_iop[j] + y1_iop[j] + z1_iop[j]) != 0 && (x2_iop[j] + y2_iop[j] + z2_iop[j]) != 0) // greater than r and less than r+dr; sum of coordinates of each points not equal to zero
							{
								// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
								// Intermolecular order parameter calculations
								vx1 = x2_iop[a] - x1_iop[a];
								vy1 = y2_iop[a] - y1_iop[a];
								vz1 = z2_iop[a] - z1_iop[a];
								vx2 = x2_iop[j] - x1_iop[j];
								vy2 = y2_iop[j] - y1_iop[j];
								vz2 = z2_iop[j] - z1_iop[j];

								dotproduct = (vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2);
								magnitude1 = (vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1);
								magnitude2 = (vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2);
								costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
								theta = acosf (costheta);
								INTERMOLECULAR_ORDER_PARAMETER = ((3 * costheta * costheta) - 1) / 2;
								sum_iop[iop_array_variable] += INTERMOLECULAR_ORDER_PARAMETER;
								denominator_iop[iop_array_variable]++;
								// fprintf(stdout, "sum_iop[%d]: %f; denominator_iop[%d]: %d; INTERMOLECULAR_ORDER_PARAMETER: %f\n", iop_array_variable, sum_iop[iop_array_variable], iop_array_variable, denominator_iop[iop_array_variable], INTERMOLECULAR_ORDER_PARAMETER);
								// fprintf(stdout, "a: %d; r -- (r+dr): (%d -- %d); j: %d; iop_array_variable: %d;\n", a, r, r+dr, j, iop_array_variable);
								// fflush (stdout);
								// sleep(1);
							}
						}
						r = r+dr;
					}
				}
			}

			// Checking IOP
			// for (int r = 0; r < 65; ++r)
			// {
				// fprintf(stdout, "r: %d --> iop: %f\n", iop, sum_iop[iop]);
				// printf("sum_iop[%d]: %f; denominator_iop[%d]: %d; avg_iop[%d]: %f\n", r, sum_iop[r], r, denominator_iop[r], r, sum_iop[r]/denominator_iop[r]);
			// }
			// sleep(1);

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// INTERMOLECULAR ORDER PARAMETER - END
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// ORIENTATIONAL ORDER PARAMETER - START
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			X1_COOP = 0;
			Y1_COOP = 0;
			Z1_COOP = 0;
			Y2_COOP = 0;

			// printf("NUMBER_OF_IOPs_actual: %d\nDOP: %d", NUMBER_OF_IOPs_actual, degree_of_polymerization);

			for (int Ltheta = 0; Ltheta < 359; ++Ltheta)
			{
				angle_degrees = Ltheta;
				angle_radians = degreeToRadians(angle_degrees);

				X2_COOP = rotateX (angle_radians);
				Z2_COOP = rotateZ (angle_radians);
				// printf("(%.5f, %.5f)\n", X2_COOP, Z2_COOP);

				vx2 = X2_COOP - X1_COOP;
				vy2 = Y2_COOP - Y1_COOP;
				vz2 = Z2_COOP - Z1_COOP;

				for (a = 0; a < NUMBER_OF_IOPs_actual; ++a)
				{
					if (x1_iop[a] + y1_iop[a] + z1_iop[a] + x2_iop[a] + y2_iop[a] + z2_iop[a] != 0) 
					{
						vx1 = x2_iop[a] - x1_iop[a];
						vy1 = y2_iop[a] - y1_iop[a];
						vz1 = z2_iop[a] - z1_iop[a];

						dotproduct = (vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2);

						magnitude1 = (vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1);
						magnitude2 = (vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2);

						costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
						theta = acosf (costheta);

						// fprintf(stdout, "check1\n");
						// fflush(stdout);

						CHAIN_ORIENTATION_ORDER_PARAMETER[Ltheta][a] = ((3 * costheta * costheta) - 1) / 2;

						// fprintf(stdout, "check2\n");
						// fflush(stdout);

						// fprintf(stdout, "CHAIN_ORIENTATION_ORDER_PARAMETER[%.0f][%d] = %.4f\n", angle_degrees, a, CHAIN_ORIENTATION_ORDER_PARAMETER[Ltheta][a]);
						// fflush(stdout);
						// sleep(1);
					}
				}
				// sleep(5);
			}

			// Finding the sum of CHAIN_ORIENTATION_ORDER_PARAMETER and tracking denominator for average
			for (int Ltheta = 0; Ltheta < 359; ++Ltheta)
			{
				// Loop over all NUMBER_OF_IOPs
				for (int a = 0; a < NUMBER_OF_IOPs_actual; ++a)
				{
					CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta] += CHAIN_ORIENTATION_ORDER_PARAMETER[Ltheta][a];
					CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT++;
					// fprintf(stdout, "CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT: %d\n", CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT);
					// fprintf(stdout, "CHAIN_ORIENTATION_ORDER_PARAMETER[%d][%d] = %.4f\n", Ltheta, a, CHAIN_ORIENTATION_ORDER_PARAMETER[Ltheta][a]);
					// fprintf(stdout, "CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[%d] = %.4f\n", Ltheta, CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta]);
				}
				// fprintf(stdout, "CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[%d] = %.4f\n", Ltheta, CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta]);
				// sleep(3);
			}
			// CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT++;
			// fprintf(stdout, "CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT: %d\n", CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT);
			// sleep(1);

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// ORIENTATIONAL ORDER PARAMETER - END
			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			// To stop the program after first time frame, then use exit(1) in next line
			// exit(1);
			// To check output by using printf statements, use the following sleep()
			// sleep(1);
			currentLine = 0;
			NUMBER_OF_IOPs_actual = 0;
		}
	}

	// Calculating the normalized value (feature scaling) of chain orientation order parameter
	// Finding the max and minimum value in CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[]
	float coop_min, coop_max, coop_normalized;
	for (int Ltheta = 0; Ltheta < 359; ++Ltheta)
	{
		if (Ltheta == 0)
		{
			coop_min = CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta];
			coop_max = coop_min;
		}
		else
		{
			if (CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta] <= coop_min)
			{
				coop_min = CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta];
			}
			else if (CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta] >= coop_max)
			{
				coop_max = CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta];
			}
		}
	}
	printf("coop_min: %f; coop_max: %f\n", coop_min, coop_max);
	// Printing normalized values
	for (int Ltheta = 0; Ltheta < 359; ++Ltheta)
	{
		coop_normalized = CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta] - coop_min;
		coop_normalized = coop_normalized / (coop_max - coop_min);
		fprintf(COOP, "%d\t%f\n", Ltheta, coop_normalized);
	}
	// sleep(10);

	// fprintf(stdout, "CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT: %d\n", CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT);
	// Calculating average chain orientation order parameter
	// for (int Ltheta = 0; Ltheta < 359; ++Ltheta)
	// {
		// fprintf(stdout, "CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[%d] = %.4f\n", Ltheta, CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta]);
		// CHAIN_ORIENTATION_ORDER_PARAMETER_AVERAGE[Ltheta] = CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta] / CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT;
		// fprintf(stdout, "CHAIN_ORIENTATION_ORDER_PARAMETER_AVERAGE[%d] = %f (count: %d)\n", Ltheta, CHAIN_ORIENTATION_ORDER_PARAMETER_AVERAGE[Ltheta], CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT);
		// fprintf(COOP, "%d\t%f\n", Ltheta, CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta]);
		// fprintf(stdout, "theta = %d --> %.4f / %d = %.4f\n", Ltheta, CHAIN_ORIENTATION_ORDER_PARAMETER_SUM[Ltheta], CHAIN_ORIENTATION_ORDER_PARAMETER_COUNT, CHAIN_ORIENTATION_ORDER_PARAMETER_AVERAGE[Ltheta]);
		// fflush (stdout);
	// }

	// Printing intermolecular order parameter
	for (int r = 0; r < nbins_iop; ++r)
	{
		// fprintf(stdout, "r: %d --> iop: %f\n", iop, sum_iop[iop]);
		// printf("sum_iop[%d]: %f; denominator_iop[%d]: %d; avg_iop[%d]: %f\n", r, sum_iop[r], r, denominator_iop[r], r, sum_iop[r]/denominator_iop[r]);
		if (denominator_iop[r] != 0)
		{
			fprintf(IOP, "%d\t%f\t%d\n", r, sum_iop[r]/denominator_iop[r], denominator_iop[r]);
		}
		else
		{
			fprintf(IOP, "%d\t%f\t%d\n", r, 0.0, denominator_iop[r]);
		}
		// fflush(stdout);
	}


	return 0;
}

/*

54.7 to 125 means negative value in order parameter
-54.7 to 54.7 means positive value in order parameter

*/