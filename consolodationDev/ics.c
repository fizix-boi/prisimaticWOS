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

	//fundamental constants
	//double h = 6.62607015 * pow(10, (-34)); //[Js]
	//double k = 1.380649 * pow(10, (-23)); //[J/K]
	//double m = 9.10938356 * pow(10, (-31)); //[kg]
	double q = 1.60217662 * pow(10, (-19)); //[C]
	//double c = 299792458.0; //[m/s]
	//double u0 = 4.0 * PI * pow(10, (-7)); //[kg*m/C^2]
	//double e0 = 1.0 / ((c * c) * u0); //[F/m]
	//double fine = (q * q * c * u0) / (2 * h); //[unitless]
	//double gr = 9.8; //[m/s^2]

	//parameters of system
	//double rho = 124.8; //[kg/m^3]
	//double n = 1 * pow(10, 14);//[#/m^2]
	//double er = 1.057; //[unitless]

	//length scales for defined problem
	double lzTM = 35 * pow(10, (-9));
	double lzDP = 50 * pow(10, (-9));
	double diffr = 15 * pow(10, (-9));
	//double lzBM = 100 * pow(10, (-9));
	double diff = diffr + 0.5 * (lzTM + lzDP);
	//double dSiN = 1500 * pow(10, (-9));
	//double zBM = dSiN + 0.5 * (lzDP + lzBM);
	double dr = (10 * pow(10, (-9)));
	double d = dr + (0.5 * lzTM);
	double lx = 5000 * pow(10, (-9));
	//double lxBM = 15000 * pow(10, (-9));
	double ly = 200 * pow(10, (-9));
	double dtd = 3000 * pow(10, (-9));

	//initialize electron
	double electronPos[] = {0.0, 0.0, d};
	ElectronClass epectron;
	electron(&epectron, -1.0 * q, electronPos);

	//initialize metal plates in system
	int numObs = 2;
	MetalSheetClass objects[numObs];

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

	MetalSheetClass dotPot;
	int dotPotNumHoles = 0;
	HoleClass dotPotHoles[dotPotNumHoles];
	double dotPotPos[] = {0.0, 0.0, -diff};
	double dotPotGeo[] = {lx, ly, lzDP};
	metalSheet(&dotPot, 0, dotPotPos, dotPotGeo, dotPotNumHoles, dotPotHoles);
	objects[1] = dotPot;

	/*
	MetalSheetClass bottomMetal;
	int bottomMetalNumHoles = 0;
	HoleClass bottomMetalHoles[bottomMetalNumHoles];
	double bottomMetalPos[] = {0.0, 0.0, -zBM};
	double bottomMetalGeo[] = {lxBM, ly, lzBM};
	metalSheet(&bottomMetal, 0, bottomMetalPos, bottomMetalGeo, bottomMetalNumHoles, bottomMetalHoles);
	objects[2] = bottomMetal;
	*/

	//main code
	double lBound = 1 * pow(10, (-11));
	double hBound = 1.0;
	int reps = pow(10, 4);
	double normer = 1 / ((double) reps);

	double x0 = (-0.5 * lx) * 1.25;
	double xf = (0.5 * lx) * 1.25;
	double perSide = 50;

	double num = (2 * perSide) + 1;
	double dx = (xf - x0) / (num - 1);

	int div = 1;
	FILE *printer;
	printer = fopen("ICSData.txt", "w");
	fprintf(printer, "Height (nm):\t%f\n", dr * 1000000000.0);
	fprintf(printer, "Position (nm)\tVoltage (mV)\n");

	int j = 0;
	((&epectron)->pos)[0] = x0;

	for(; j < num; j++){
		double volt = 0;
		int i = 0;

		for(; i < reps; i++){
			volt += potentialShot(epectron, objects, numObs, lBound, hBound);
			if(((i + 1) % (reps / div)) == 0){
				printf("X");
			}
		}
		printf("> #%3d/%d done\n", j + 1, (int) num);

		volt *= normer;

		fprintf(printer, "%f\t%f\n", ((&epectron)->pos)[0] * 1000000000.0, volt * 1000.0);

		((&epectron)->pos)[0] += dx;
	}

	//ending
	fclose(printer);

	return 0;
}