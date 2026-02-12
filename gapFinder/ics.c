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
	double lz = 50 * pow(10, (-9));
	double lx = 300000 * pow(10, (-9));
	double ly = 500000 * pow(10, (-9));
	
	double d = (300 * pow(10, (-9))) + (0.5 * lz);
	double g = 500 * pow(10, (-9));

	//initialize electron
	double electronPos[] = {0.0, 0.0, 0.0};
	ElectronClass epectron;
	double q = 1.60217662 * pow(10, (-19)); //[C]
	electron(&epectron, -1.0 * q, electronPos);

	//initialize metal plates in system
	int numObs = 2;
	MetalSheetClass objects[numObs];
	HoleClass nullHoles[0];

	//defining gates

	double V0 = 0.0;

	//positive sloshing metal
	MetalSheetClass posSlosh;
	int posSloshNumHoles = 0;
	HoleClass posSloshHoles[posSloshNumHoles];
	double posSloshPos[] = {0, 0.0, 0.0};
	double posSloshGeo[] = {lx, ly, lz};
	metalSheet(&posSlosh, V0, posSloshPos, posSloshGeo, posSloshNumHoles, posSloshHoles);
	objects[0] = posSlosh;

	//negative sloshing metal
	MetalSheetClass negSlosh;
	int negSloshNumHoles = 0;
	HoleClass negSloshHoles[negSloshNumHoles];
	double negSloshPos[] = {0, 0.0, 0.0};;
	double negSloshGeo[] = {lx, ly, lz};
	metalSheet(&negSlosh, V0, negSloshPos, negSloshGeo, negSloshNumHoles, negSloshHoles);
	objects[1] = negSlosh;

	//END ORIENTATION DEFINITION

	//main code
	double lBound = 1 * pow(10, (-9));
	double hBound = 1.0;
	int reps = pow(10, 8);
	int div = 20;
	double normer = 1 / ((double) reps);
	
	double z0 = (100 * pow(10, (-9))) + (0.5 * lz);
	double zf = (600 * pow(10, (-9))) + (0.5 * lz);
	double num = 10 + 1;
	double dz = (zf - z0) / (num - 1);
	
	double g0 = (90 * pow(10, (-9)));
	double gf = (10 * pow(10, (-9)));
	double numg = 4 + 1;
	double dg = (gf - g0) / (numg - 1);
	
	FILE *printer;
	printer = fopen("ICSData.txt", "w");
	
	int i,j,l;
	
	double gapCurr = g0;
	
	l = 0;
	
	for(; l < numg; l++){
		printf("Gap (nm):\t%f\n", gapCurr * 1000000000.0);
		
		((&posSlosh)->pos)[0] = -0.5 * (lx + gapCurr);
		((&negSlosh)->pos)[0] = 0.5 * (lx + gapCurr);
		
		fprintf(printer, "Gap (nm):\t%f\n", gapCurr * 1000000000.0);
		fprintf(printer, "Position (nm)\tVoltage (mV)\n");
		
		j = 0;
		
		((&epectron)->pos)[2] = z0;
		
		for(; j < num; j++){
			((&epectron)->pos)[0] = ((&posSlosh)->pos)[0];
			double base = 0;
			i = 0;
			
			for(; i < reps; i++){
				base += potentialShot(epectron, objects, numObs, lBound, hBound);
				if(((i + 1) % (reps / div)) == 0){
					printf("X");
				}
			}
			printf("> #%3d/%d done plate\n", j + 1, (int) num);
			
			base *= normer;
			
			((&epectron)->pos)[0] = 0;
			double gapper = 0;
			i = 0;
			
			for(; i < reps; i++){
				gapper += potentialShot(epectron, objects, numObs, lBound, hBound);
				if(((i + 1) % (reps / div)) == 0){
					printf("X");
				}
			}
			printf("> #%3d/%d done gap\n", j + 1, (int) num);
			
			gapper *= normer;
			
			double barrier = base - gapper;
			
			fprintf(printer, "%f\t%f\n", (((&epectron)->pos)[2] - (0.5 * lz)) * 1000000000.0, barrier * 1000.0);
			
			((&epectron)->pos)[2] += dz;
		}
		
		gapCurr += dg;
		fprintf(printer, "\n");
	}
	
	fclose(printer);
	
	//ending
	return 0;
}
