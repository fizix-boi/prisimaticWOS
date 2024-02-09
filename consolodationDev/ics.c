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
	double hBottom = 0.87 * pow(10, (-6));
	double hTop = 0.2 * pow(10, (-6));
	double dr =  100 * pow(10, (-9));

	//l=y, w=x

	double lGate = 2 * pow(10, (-6));
	double wGate = 0.5 * pow(10, (-6));

	double lOpen = (8.0 + 8.0 + 1.0) * wGate;
	double wOpen = wGate;

	double negXLength = 50 * pow(10, (-6));
	double negYLength = 100 * pow(10, (-6));

	double x0 = -3 * (lOpen * 0.5);

	//initialize electron
	double electronPos[] = {0.0, 0.0, dr};
	ElectronClass epectron;
	double q = 1.60217662 * pow(10, (-19)); //[C]
	electron(&epectron, -1.0 * q, electronPos);

	//initialize metal plates in system
	int numObs = 12;
	MetalSheetClass objects[numObs];
	HoleClass nullHoles[0];
	double poser[numObs][3];
	double geoer[numObs][3];

	//defining frame

	//negative x metal
	MetalSheetClass negX;
	int negXNumHoles = 0;
	HoleClass negXHoles[negXNumHoles];
	double negXPos[] = {-0.5 * (negXLength + lOpen), 0.0, -hTop * 0.5};
	double negXGeo[] = {negXLength, lGate, hTop};
	metalSheet(&negX, 0, negXPos, negXGeo, negXNumHoles, negXHoles);
	objects[8] = negX;

	//positive x metal
	MetalSheetClass posX;
	int posXNumHoles = 0;
	HoleClass posXHoles[posXNumHoles];
	double posXPos[] = {0.5 * (negXLength + lOpen), 0.0, -hTop * 0.5};
	double posXGeo[] = {negXLength, lGate, hTop};
	metalSheet(&posX, 0, posXPos, posXGeo, posXNumHoles, posXHoles);
	objects[9] = posX;

	//negative y metal
	MetalSheetClass negY;
	int negYNumHoles = 0;
	HoleClass negYHoles[negYNumHoles];
	double negYPos[] = {0.0, -0.5 * (negYLength + wOpen), -hTop * 0.5};
	double negYGeo[] = {2 * negXLength + lOpen, negYLength, hTop};
	metalSheet(&negY, 0, negYPos, negYGeo, negYNumHoles, negYHoles);
	objects[10] = negY;

	//positive y metal
	MetalSheetClass posY;
	int posYNumHoles = 0;
	HoleClass posYHoles[posYNumHoles];
	double posYPos[] = {0.0, 0.5 * (negYLength + wOpen), -hTop * 0.5};
	double posYGeo[] = {2 * negXLength + lOpen, negYLength, hTop};
	metalSheet(&posY, 0, posYPos, posYGeo, posYNumHoles, posYHoles);
	objects[11] = posY;

	//defining twiddles

	double numTwid = 8;

	int j = 0;
	for(; j < numTwid; j++){
		MetalSheetClass twiddle;
		poser[j][0] = (-7.0 + (2.0 * j)) * wGate;
		poser[j][1] = 0.0;
		poser[j][2] = -hBottom * 0.5 - hTop;
		double gateGeo[] = {wGate, lGate, hBottom};
		metalSheet(&(objects[j]), 0, poser[j], gateGeo, 0, nullHoles);
		//print(objects[j]);
	}

	//END ORIENTATION DEFINITION

	//main code

	//define bounding box
	double lBound = 1.0 * pow(10, (-11));
	double hBound = 1.0;

	//define the number of shots per data point
	int reps = pow(10, 4);
	double normer = 1 / ((double) reps);

	//give the parameters of the position sweep
	double xf = -1 * x0;
	double perSide = 500; //defines the number of points per side

	double num = (2 * perSide) + 1; //makes sweep even
	double dx = (xf - x0) / (num - 1); //find the step size in the simulation

	int div = 1; //prints this many X's to show progress of the simulation

	//initialize the file for data collection
	FILE *printer;
	printer = fopen("ICSData.txt", "w");
	fprintf(printer, "Height (nm):\t%f\n", dr * 1000000000.0);
	fprintf(printer, "Position (nm)\tVoltage (mV)\n");

	//initialiize the electron for iteracting through the loop
	j = 0;
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