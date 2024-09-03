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

	//START ORIENTATION DEFINITION

	//length scales for defined problem
	double lGate = 100 * pow(10, (-6));
	double wGate = 100 * pow(10, (-6));
	double hGate = 100 * pow(10, (-6));

	double dElec = (10) * pow(10, (-9)) + (hGate / 2);

	//initialize electron
	double electronPos[] = {0.0, 0.0, dElec};
	ElectronClass epectron;
	double q = 1.60217662 * pow(10, (-19)); //[C]
	electron(&epectron, -1.0 * q, electronPos);

	//initialize metal plates in system
	int numObs = 1;
	MetalSheetClass objects[numObs];
	HoleClass nullHoles[0];

	//defining gates

	double V0 = 0.0;

	//positive sloshing metal
	MetalSheetClass metaa;
	HoleClass metaaHoles[0];
	double metaaPos[] = {0.0, 0.0, 0.0};
	double metaaGeo[] = {lGate, wGate, hGate};
	metalSheet(&metaa, V0, metaaPos, metaaGeo, 0, metaaHoles);
	objects[0] = metaa;

	//END ORIENTATION DEFINITION

	//main code

	//define bounding box
	double lBound = 1.0 * pow(10, (-11));
	double hBound = 1.0;

	//define the number of shots per data point
	int reps = pow(10, 6);
	double normer = 1 / ((double) reps);

	//give the parameters of the position sweep
	double x0 = 0.95 * (hGate / 2);//100 * pow(10, (-6));
	double xf = 10 * (hGate / 2);//100 * pow(10, (-6));
	double perSide = 1500; //defines the number of points per side

	double num = (2 * perSide) + 1; //makes sweep even
	double dx = (xf - x0) / (num - 1); //find the step size in the simulation

	int div = 10; //prints this many X's to show progress of the simulation

	//initialize the file for data collection
	FILE *printer;
	printer = fopen("ICSData.txt", "w");
	fprintf(printer, "Height (nm):\t%f\n", (dElec - (hGate / 2)) * 1000000000.0);
	fprintf(printer, "Position (nm)\tVoltage (mV)\n");

	//initialiize the electron for iteracting through the loop
	int j;
	j = 0;
	((&epectron)->pos)[0] = x0;

	double dt_avg = 0;

	for(; j < num; j++){ //for the number of points in the simulation
		double volt = 0; //initializes the res
		int i = 0;

		double timer0 = (double) clock();
		for(; i < reps; i++){ //for the number of monte carlo shots requested
			//calculate the potential shot and add it the the result for later
			volt += potentialShot(epectron, objects, numObs, lBound, hBound);
			if(((i + 1) % (reps / div)) == 0){ //if the current iteration is a multiple of the div variable
				printf("X"); //print a progress X
			}
		}
		double dt = (((double) clock()) - timer0) / CLOCKS_PER_SEC;

		dt_avg = ((dt_avg * j) + dt) / (j + 1);
		int remainer = num - j - 1;
		double time_remainer = dt_avg * remainer / 60;

		printf("> #%3d/%d done - est. time reamaining: %f min\n", j + 1, (int) num, time_remainer);

		volt *= normer; //normalize to the number of shots calculated

		//output the calculation to a file
		fprintf(printer, "%f\t%f\n", ((&epectron)->pos)[0] * 1000000000.0, volt * 1000.0);

		((&epectron)->pos)[0] += dx; //move the electron to the next position for further calculation
	}

	//end of code, close and free all necessary objects and return 0
	fclose(printer);
	return 0;
}
