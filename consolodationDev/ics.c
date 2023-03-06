#include <stdbool.h>
#include <stdlib.h>
#include <time.h>

#include "icsWOSMethods.h"

/*
icsWOSMethods.h includes funcs3d.h includes stdio.h and math.h
*/

#define PI 3.1415926535

int main(void){
	srand(time(NULL)); //yo okay, this needs to be at the beginning of main or else things won't work

	//length scales for defined problem
	double lzTM = 35 * pow(10, (-9));
	double lzDP = 50 * pow(10, (-9));
	double diffr = 15 * pow(10, (-9));
	double diff = diffr + 0.5 * (lzTM + lzDP);
	double lzBM = 100 * pow(10, (-9));
	double dSiN = 1500 * pow(10, (-9));
	double dr = (10 * pow(10, (-9)));
	double d = dr + (0.5 * lzTM);
	double lx = 5000 * pow(10, (-9));
	double ly = 200 * pow(10, (-9));
	double dtd = 3000 * pow(10, (-9));

	//initialize electron
	double electronPos[] = {0.0, 0.0, d};
	ElectronClass epectron;
	double q = 1.60217662 * pow(10, (-19)); //[C]
	electron(&epectron, -1.0 * q, electronPos);

	//initialize metal plates in system
	int numObs = 2;
	MetalSheetClass objects[numObs];

	//top metal
	MetalSheetClass topMetal;
	int topMetalNumHoles = 2;
	HoleClass topMetalHoles[topMetalNumHoles];
	double topMetalH0Pos[] = {-0.5 * dtd, 0.0, 0.0}; //note that this last arg wont matter
	hole(&(topMetalHoles[0]), topMetalH0Pos, 50 * pow(10, (-9)));
	double topMetalH1Pos[] = {0.5 * dtd, 0.0, 0.0}; //note that this last arg wont matter
	hole(&(topMetalHoles[1]), topMetalH1Pos, 50 * pow(10, (-9)));
	double topMetalPos[] = {0.0, 0.0, 0.0};
	double topMetalGeo[] = {lx, ly, lzTM};
	metalSheet(&topMetal, 0, topMetalPos, topMetalGeo, topMetalNumHoles, topMetalHoles);
	objects[0] = topMetal;

	//dot potential metal
	MetalSheetClass dotPot;
	int dotPotNumHoles = 0;
	HoleClass dotPotHoles[dotPotNumHoles];
	double dotPotPos[] = {0.0, 0.0, -diff};
	double dotPotGeo[] = {lx, ly, lzDP};
	metalSheet(&dotPot, 0, dotPotPos, dotPotGeo, dotPotNumHoles, dotPotHoles);
	objects[1] = dotPot;

	/*
	//bottom metal (adds marginal accuracy, can omit)
	MetalSheetClass bottomMetal;
	int bottomMetalNumHoles = 0;
	HoleClass bottomMetalHoles[bottomMetalNumHoles];

	double zBM = dSiN + 0.5 * (lzDP + lzBM);
	double lxBM = 15000 * pow(10, (-9));
	double bottomMetalPos[] = {0.0, 0.0, -zBM};
	double bottomMetalGeo[] = {lxBM, ly, lzBM};

	metalSheet(&bottomMetal, 0, bottomMetalPos, bottomMetalGeo, bottomMetalNumHoles, bottomMetalHoles);
	objects[2] = bottomMetal;
	//*/

	//main code

	//define bounding box
	double lBound = 1 * pow(10, (-11));
	double hBound = 1.0;

	//define the number of shots per data point
	int reps = pow(10, 4);
	double normer = 1 / ((double) reps);

	//give the parameters of the position sweep
	double x0 = (-0.5 * lx) * 1.25;
	double xf = (0.5 * lx) * 1.25;
	double perSide = 50; //defines the number of points per side

	double num = (2 * perSide) + 1; //makes sweep even
	double dx = (xf - x0) / (num - 1); //find the step size in the simulation

	int div = 1; //prints this many X's to show progress of the simulation

	//initialize the file for data collection
	FILE *printer;
	printer = fopen("ICSData.txt", "w");
	fprintf(printer, "Height (nm):\t%f\n", dr * 1000000000.0);
	fprintf(printer, "Position (nm)\tVoltage (mV)\n");

	//initialiize the electron for iteracting through the loop
	int j = 0;
	((&epectron)->pos)[0] = x0;

	for(; j < num; j++){ //for the number of points in the simulation
		double volt = 0; //initializes the res
		int i = 0;

		for(; i < reps; i++){ //for the number of monte carlo shots requested
			//calculate the potential shot and add it the the result for later
			volt += potentialShot(epectron, objects, numObs, lBound, hBound);
			if(((i + 1) % (reps / div)) == 0){ //if the current iteration is a multiple of the div variable
				printf("X"); //print a progress X
			}
		}
		printf("> #%3d/%d done\n", j + 1, (int) num);

		volt *= normer; //normalize to the number of shots calculated

		//output the calculation to a file
		fprintf(printer, "%f\t%f\n", ((&epectron)->pos)[0] * 1000000000.0, volt * 1000.0);

		((&epectron)->pos)[0] += dx; //move the electron to the next position for further calculation
	}

	//end of code, close and free all necessary objects and return 0
	fclose(printer);
	return 0;
}