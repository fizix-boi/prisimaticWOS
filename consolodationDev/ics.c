#include <stdio.h>
#define PI 3.1415926535
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>


//object declarations

typedef struct {
	double Q;
	double* pos;
} ElectronClass;

void electron(ElectronClass *elec, double Q, double* pos){
	elec->Q = Q;
	elec->pos = pos;
}

typedef struct {
	double* pos;
	double R;
} HoleClass;

typedef struct {
	double V;
	double* pos;
	double* geo;
	int numHoles;
	HoleClass* holes;
} MetalSheetClass;

void hole(HoleClass *hol, double* pos, double R){
	hol->pos = pos;
	hol->R = R;
}

void metalSheet(MetalSheetClass *meta, double V, double* pos, double* geo, int numHoles, HoleClass* holes){
	meta->V = V;
	meta->pos = pos;
	meta->geo = geo;
	meta->numHoles = numHoles;
	meta->holes = holes;
}

//helper functions
void print3Arr(double* threeArr){
	printf("(%e, %e, %e)\n", threeArr[0], threeArr[1], threeArr[2]);
}

void add3Arr(double* res, double* arr1, double* arr2){
	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = arr1[i] + arr2[i];
	}
}

void subtract3Arr(double* res, double* arr1, double* arr2){
	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = arr1[i] - arr2[i];
	}
}

void multiply3Arr(double* res, double* arr1, double* arr2){
	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = arr1[i] * arr2[i];
	}
}

void constMultiply3Arr(double* res, double num, double* arr){
	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = num * arr[i];
	}
}

void sign3Arr(double* res, double* arr){
	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = (arr[i] > 0) ? 1.0 : -1.0;
	}
}

void equal3Arr(double* res, double* arr){
	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = arr[i];
	}
}

double mag3Arr(double* arr){
	return(sqrt((arr[0] * arr[0]) + (arr[1] * arr[1]) + (arr[2] * arr[2])));
}

//function headers
void findClosestPoint(double*, double*, MetalSheetClass);
double potentialShot(ElectronClass, MetalSheetClass*, int, double, double);
int isAboveHole(double*, MetalSheetClass);

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
	double perSide = 10;

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

void findClosestPoint(double* closestPoint, double* posUse, MetalSheetClass objecti){
	double* center = (&objecti)->pos;
	double* geom = (&objecti)->geo;
	double halfGeom[3]; constMultiply3Arr(halfGeom, 0.5, geom);
	double norm[3]; subtract3Arr(norm, posUse, center);
	double unit[] = {1.0, 1.0, 1.0};

	int holy = isAboveHole(posUse, objecti);

	if(!holy){
		bool isInX = (fabs(norm[0]) < (halfGeom[0]));
		bool isInY = (fabs(norm[1]) < (halfGeom[1]));
		bool isInZ = (fabs(norm[2]) < (halfGeom[2]));
		double planeInners[] = {isInX?1.0:0.0, isInY?1.0:0.0, isInZ?1.0:0.0};
		double planeOuters[3]; subtract3Arr(planeOuters, unit, planeInners);
		double signs[3]; sign3Arr(signs, norm);

		double fromPos[3]; multiply3Arr(fromPos, planeInners, posUse);
		double fromGeo[3]; multiply3Arr(fromGeo, signs, halfGeom); add3Arr(fromGeo, fromGeo, center); multiply3Arr(fromGeo, planeOuters, fromGeo);
		add3Arr(closestPoint, fromPos, fromGeo);
	}else{
		HoleClass* holes = (&objecti)->holes;
		int ind = abs(holy) - 1;
		int isIn = (1 - (holy / abs(holy))) / 2;

		double* locCurr = (&(holes[ind]))->pos;
		double differ[3]; subtract3Arr(differ, norm, locCurr);
		differ[2] = 0;
		double mag = mag3Arr(differ);

		double holyLoc[3]; equal3Arr(holyLoc, differ);
		constMultiply3Arr(holyLoc, ((&(holes[ind]))->R) / mag, holyLoc);
		if(!isIn){
			holyLoc[2] = halfGeom[2] * (norm[2] / fabs(norm[2]));
		}else{
			holyLoc[2] = norm[2];
		}
		add3Arr(closestPoint, holyLoc, locCurr);
	}
}

int isAboveHole(double* posUse, MetalSheetClass objecti){
	double* center = (&objecti)->pos;
	double norm[3]; subtract3Arr(norm, posUse, center);

	int numHole = (&objecti)->numHoles;
	HoleClass* holes = (&objecti)->holes;
	int res = 0;

	int i;

	for(i = 0; i < numHole; i++){
		double* locCurr = (&(holes[i]))->pos;
		double differ[3]; subtract3Arr(differ, norm, locCurr);
		differ[2] = 0;
		double mag = mag3Arr(differ);
		int isAbove = (mag < (&(holes[i]))->R)?1:0;
		int isIn = (fabs(norm[2]) < (0.5 * ((&objecti)->geo)[2]))?1:0;
		if(isAbove){
			res = i + 1;
			if(isIn){
				res = -1 * res;
			}
			break;
		}
	}

	return(res);
}

double potentialShot(ElectronClass charge, MetalSheetClass* objects, int numObs, double lBound, double hBound){
	double c = 299792458.0; //[m/s]
	double u0 = 4.0 * PI * pow(10, (-7)); //[kg*m/C^2]
	double e0 = 1.0 / ((c * c) * u0); //[F/m]
	double posUse[3]; equal3Arr(posUse, (&charge)->pos);

	while(true){
		double rPoints[numObs];
		double rPointsEStats[numObs];
		int i;

		for(i = 0; i < numObs; i++){
			double point[3]; findClosestPoint(point, posUse, objects[i]);
			double diff[3]; subtract3Arr(diff, posUse, point);
			double diffEStat[3]; subtract3Arr(diffEStat, (&charge)->pos, point);
			rPoints[i] = mag3Arr(diff);
			rPointsEStats[i] = mag3Arr(diffEStat);
		}
		int mini = 0;
		double r = rPoints[mini];
		for (i = 0; i < numObs; ++i){
		    if (rPoints[i] < r){
		        r = rPoints[i];
		        mini = i;
		    }
		}
		if((r < lBound) || (r > hBound)){
			double eStatConst = 0.5 * ((&charge)->Q) / (4 * PI * e0);
			return(((&(objects[mini]))->V) - (eStatConst / rPointsEStats[mini]));
		}
		double phi, theta;
		double irm = 1 / ((double) RAND_MAX);

    	phi = (((double)rand()) * irm) * (2.0 * PI);
    	theta = acos((2.0 * (((double)rand()) * irm)) - 1.0);
    	double newDir[] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
		double newDisp[3]; constMultiply3Arr(newDisp, r, newDir);
		add3Arr(posUse, newDisp, posUse);
	}

	return(-1);
}
