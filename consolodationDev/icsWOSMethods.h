#include "funcs3d.h"
#include "icsPrimatives.h"

/*
funcs3d.h includes stdio.h and math.h
*/

#define PI 3.1415926535

//function headers
void findClosestPoint(double*, double*, MetalSheetClass);
int isAboveHole(double*, MetalSheetClass);
double potentialShot(ElectronClass, MetalSheetClass*, int, double, double);

void findClosestPoint(double* closestPoint, double* posUse, MetalSheetClass objecti){
    /* 
	 * FUNCTION: findClosestPoint
	 * ------------------------------
	 * 
	 * finds the point closest to a given position posUse on a given MetalSheetClass object objecti
	 * 
	 * closestPoint: the 3-array where the closest point will be cast at the end of the function
     * posUse: the position which will act as the origin of the search for the closest point
     * objecti: the object of which the closest point will be on
	 * 
	 * returns: void
	 */

    //I found that the position and geometry vector needed to be recast, probably
    //due to some indexing I wasn't being careful about
	double* center = (&objecti)->pos;
	double* geom = (&objecti)->geo;

	double halfGeom[3]; constMultiply3Arr(halfGeom, 0.5, geom); //this will be useful to know for later
	double norm[3]; subtract3Arr(norm, posUse, center); //this sets the metal sheet's center as the (0, 0, 0) of the problem
	double unit[] = {1.0, 1.0, 1.0}; //unitary vector for later use

	int holy = isAboveHole(posUse, objecti); //tells the index of the hole (i + 1), whether or not you are in the hole (positive or neagtive) (0 means not in or above any hole)

	if(!holy){ //if the electric charge is NOT above or in a hole
        /*
        These functions see if the x, y, and z coordinates are within the bounds of the geometry.
        This helps to know if the position for the charge can simply be projected onto the sheet or not
        */
		bool isInX = (fabs(norm[0]) < (halfGeom[0]));
		bool isInY = (fabs(norm[1]) < (halfGeom[1]));
		bool isInZ = (fabs(norm[2]) < (halfGeom[2]));

		double planeInners[] = {isInX?1.0:0.0, isInY?1.0:0.0, isInZ?1.0:0.0}; //converts coordinate bools to a 1-0 logic 3-array
		double planeOuters[3]; subtract3Arr(planeOuters, unit, planeInners); //logical counter of planeInners
		double signs[3]; sign3Arr(signs, norm); //finds which side of the sheet the charge is on and sends those signs to an array

        /*
        Constructs the closest point via the position components and the geometry components.

        If the charge's position in one of the directions can directly be projected onto the metal sheet,
        than that position can simply be used in the closest point description.

        If the charge's position in one of the directions cannot be directly projected, then it must use the geometry
        to fill in the closest point for the direction. It does this by taking the half-geometry vector, multiplying it
        by the signs vector (to get the approriate signs on the geometry based on what octet of the prisms coordinates
        the charge is in), adding the center of the prism to revert to the original coordates of the problem, and then
        multiplies by planeOuters to ensure that only the components that still need the geometric components are
        filled.

        The sum of these two, then, yields the final resultant closest point
        */
		double fromPos[3]; multiply3Arr(fromPos, planeInners, posUse);
		double fromGeo[3]; multiply3Arr(fromGeo, signs, halfGeom); add3Arr(fromGeo, fromGeo, center); multiply3Arr(fromGeo, planeOuters, fromGeo);
		add3Arr(closestPoint, fromPos, fromGeo);
	}else{ //if the electric charge IS above or in a hole
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
    /* 
	 * FUNCTION: electron
	 * ------------------------------
	 * 
	 * creates an elcetric charge of charge Q at position pos using the ElectronClass struct
	 * 
	 * elec: the charge object that the parameters will be set to
     * Q: the charge of the electric charge
     * pos: the position as a 3-array of the electric charge
	 * 
	 * returns: void
	 */

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
    /* 
	 * FUNCTION: electron
	 * ------------------------------
	 * 
	 * creates an elcetric charge of charge Q at position pos using the ElectronClass struct
	 * 
	 * elec: the charge object that the parameters will be set to
     * Q: the charge of the electric charge
     * pos: the position as a 3-array of the electric charge
	 * 
	 * returns: void
	 */

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