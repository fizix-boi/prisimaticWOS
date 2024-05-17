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
	double hGate = 25 * pow(10, (-9));
	double tGate = 100 * pow(10, (-9)); //spacing

	double dElec = (11 + 200) * pow(10, (-9)); //spacing

	//initialize electron
	double electronPos[] = {0.0, 0.0, dElec};
	ElectronClass epectron;
	double q = 1.60217662 * pow(10, (-19)); //[C]
	electron(&epectron, -1.0 * q, electronPos);

	//initialize metal plates in system
	int numObs = 3;
	MetalSheetClass objects[numObs];
	HoleClass nullHoles[0];

	//defining gates

	double V0 = 0.015;

	//positive sloshing metal
	MetalSheetClass posSlosh;
	int posSloshNumHoles = 0;
	HoleClass posSloshHoles[posSloshNumHoles];
	double posSloshPos[] = {-lGate - tGate, 0.0, -0.5 * hGate};
	double posSloshGeo[] = {lGate, wGate, hGate};
	metalSheet(&posSlosh, V0, posSloshPos, posSloshGeo, posSloshNumHoles, posSloshHoles);
	objects[0] = posSlosh;

	//ground sensing metal
	MetalSheetClass sense;
	int senseNumHoles = 0;
	HoleClass senseHoles[senseNumHoles];
	double sensePos[] = {0.0, 0.0, -0.5 * hGate};
	double senseGeo[] = {lGate, wGate, hGate};
	metalSheet(&sense, 0, sensePos, senseGeo, senseNumHoles, senseHoles);
	objects[1] = sense;

	//negative sloshing metal
	MetalSheetClass negSlosh;
	int negSloshNumHoles = 0;
	HoleClass negSloshHoles[negSloshNumHoles];
	double negSloshPos[] = {lGate + tGate, 0.0, -0.5 * hGate};
	double negSloshGeo[] = {lGate, wGate, hGate};
	metalSheet(&negSlosh, -1 * V0, negSloshPos, negSloshGeo, negSloshNumHoles, negSloshHoles);
	objects[2] = negSlosh;

	//END ORIENTATION DEFINITION

	//main code

	//define bounding box
	double lBound = 1.0 * pow(10, (-11));
	double hBound = 1.0;

	//define the number of shots per data point
	int reps = pow(10, 6);
	double normer = 1 / ((double) reps);

	//give the parameters of the position sweep
	double x0 = -51 * pow(10, (-6));
	double xf = -49 * pow(10, (-6));
	double perSide = 100; //defines the number of points per side

	double num = (2 * perSide) + 1; //makes sweep even
	double dx = (xf - x0) / (num - 1); //find the step size in the simulation

	int div = 10; //prints this many X's to show progress of the simulation

	//initialize the file for data collection
	FILE *printer;
	printer = fopen("ICSData.txt", "w");
	fprintf(printer, "Height (nm):\t%f\n", dElec * 1000000000.0);
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
