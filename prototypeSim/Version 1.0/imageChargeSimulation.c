#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>


//object declatations

typedef struct {
	double Q;
	double* pos;
} ElectronClass;

void electron(ElectronClass *elec, double Q, double* pos){
	elec->Q = Q;
	elec->pos = pos;
}

typedef struct {
	double V;
	double* pos;
	double* geo;
} MetalSheetClass;

void metalSheet(MetalSheetClass *meta, double V, double* pos, double* geo){
	meta->V = V;
	meta->pos = pos;
	meta->geo = geo;
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

int main(void){
	srand(time(NULL)); //yo okay, this needs to be at the beginning of main or else things won't work
	
	//fundamental constants
	double h = 6.62607015 * pow(10, (-34)); //[Js]
	double k = 1.380649 * pow(10, (-23)); //[J/K]
	double m = 9.10938356 * pow(10, (-31)); //[kg]
	double q = 1.60217662 * pow(10, (-19)); //[C]
	double c = 299792458.0; //[m/s]
	double u0 = 4.0 * M_PI * pow(10, (-7)); //[kg*m/C^2]
	double e0 = 1.0 / ((c * c) * u0); //[F/m]
	double fine = (q * q * c * u0) / (2 * h); //[unitless]
	double gr = 9.8; //[m/s^2]
	
	//parameters of system	
	double rho = 124.8; //[kg/m^3]
	double n = 1 * pow(10, 14);//[#/m^2]
	double er = 1.057; //[unitless]
	
	//length scales for defined problem
	double lz = 300 * pow(10, (-9));
	double d = (200 * pow(10, (-9))) + (0.5 * lz);
	double lx = 300000 * pow(10, (-9));
	double ly = 500000 * pow(10, (-9));
	double g = 100 * pow(10, (-9));
	
	double hg = 100 * pow(10, (-9));
	double wc = 4000 * pow(10, (-9));
	double hc = d - (0.5 * hg);
	
	double lg = (2 * lx) + g;
	double wg = 0.5 * wc;
	
	double diph = ((wc * wc) / (8 * fine)) * (rho * gr * d);
	double dipe = ((wc * wc) / (8 * fine)) * ((n * n * q * q) / (2 * er * e0));
	double dip = diph + dipe;
	
	printf("%e\n%e=%e+%e\n\n", d, dip, dipe, diph);
	
	double disCent = 0.5 * (lx + g);
	
	//initialize electron
	double electronPos[] = {0.0, 0.0, d};
	ElectronClass epectron;
	electron(&epectron, -1.0 * q, electronPos);
	
	//initialize metal plates in system
	int numObs = 2;
	MetalSheetClass objects[numObs];
	
	double lead1Pos[] = {-disCent, 0.0, 0.0};
	double lead1Geo[] = {lx, ly, lz};
	MetalSheetClass lead1;
	metalSheet(&lead1, 0, lead1Pos, lead1Geo);
	objects[0] = lead1;
	
	///*
	double lead2Pos[] = {disCent, 0.0, 0.0};
	double lead2Geo[] = {lx, ly, lz};
	MetalSheetClass lead2;
	metalSheet(&lead2, 0, lead2Pos, lead2Geo);
	objects[1] = lead2;
	//*/
	
	/*
	double guard1Pos[] = {0.0, -0.5 * (wg + wc), hc};
	double guard1Geo[] = {lg, wg, hg};
	MetalSheetClass guard1;
	metalSheet(&guard1, 0, guard1Pos, guard1Geo);
	objects[2] = guard1;
	//*/
	
	/*
	double guard2Pos[] = {0.0, 0.5 * (wg + wc), hc};
	double guard2Geo[] = {lg, wg, hg};
	MetalSheetClass guard2;
	metalSheet(&guard2, 0, guard2Pos, guard2Geo);
	objects[3] = guard2;
	//*/
	
	//main code
	double lBound = 1 * pow(10, (-9));
	double hBound = 1.0;
	int reps = pow(10, 7);
	int div = 100;
	double normer = 1 / ((double) reps);
	
	double x0 = -g * 6;
	double xf = -x0;
	double perSide = 30;
	double num = (2 * perSide) + 1;
	double dx = (xf - x0) / (num - 1);
	
	FILE *printer;
	printer = fopen("ICSData.txt", "w");
	fprintf(printer, "Height (nm):\t%f\n", (d - (0.5 * lz)) * 1000000000.0);
	fprintf(printer, "Gap (nm):\t%f\n", g * 1000000000.0);
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
	double unit[] = {1.0, 1.0, 1.0};
	
	double norm[3]; subtract3Arr(norm, posUse, center);
	bool isInX = (fabs(norm[0]) < (halfGeom[0]));
	bool isInY = (fabs(norm[1]) < (halfGeom[1]));
	bool isInZ = (fabs(norm[2]) < (halfGeom[2]));
	double planeInners[] = {isInX?1.0:0.0, isInY?1.0:0.0, isInZ?1.0:0.0};
	double planeOuters[3]; subtract3Arr(planeOuters, unit, planeInners);
	double signs[3]; sign3Arr(signs, norm);
	
	double fromPos[3]; multiply3Arr(fromPos, planeInners, posUse);
	double fromGeo[3]; multiply3Arr(fromGeo, signs, halfGeom); add3Arr(fromGeo, fromGeo, center); multiply3Arr(fromGeo, planeOuters, fromGeo);
	add3Arr(closestPoint, fromPos, fromGeo);
}

double potentialShot(ElectronClass charge, MetalSheetClass* objects, int numObs, double lBound, double hBound){
	double c = 299792458.0; //[m/s]
	double u0 = 4.0 * M_PI * pow(10, (-7)); //[kg*m/C^2]
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
			double eStatConst = ((&charge)->Q) / (4 * M_PI * e0);
			return(((&(objects[mini]))->V) - (eStatConst / rPointsEStats[mini]));
		}
		double phi, theta;
		double irm = 1 / ((double) RAND_MAX);

    	phi = (((double)rand()) * irm) * (2.0 * M_PI);
    	theta = acos((2.0 * (((double)rand()) * irm)) - 1.0);
    	double newDir[] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
		double newDisp[3]; constMultiply3Arr(newDisp, r, newDir);
		add3Arr(posUse, newDisp, posUse);
	}
}
