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
		HoleClass* holes = (&objecti)->holes; //set the given hole array to a local array
		int ind = abs(holy) - 1; //find the array index of the hole you are above
		int isIn = (1 - (holy / abs(holy))) / 2; //find if you are in the hole (negative holy)

		double* locCurr = (&(holes[ind]))->pos; //sets the location of the center of the hole 
		double differ[3]; subtract3Arr(differ, norm, locCurr); //finds the difference between the location of the hole and charge
        /*
        Note that, the hole locations are centered with the metal sheet, so both vectors here do, indeed, have the same origin
        */

		differ[2] = 0; //this sets the z-component to 0 and essentially makes this a 2-dimensional problem from here
		double mag = mag3Arr(differ); //finds the magnitude of the 

		double holyLoc[3]; equal3Arr(holyLoc, differ); //idk why, but this needed to be cast to a different variable to work
		constMultiply3Arr(holyLoc, ((&(holes[ind]))->R) / mag, holyLoc); //this ensures that the closest position is normalized to the radius of the hole
		if(!isIn){ //if the position is in the hole
			holyLoc[2] = halfGeom[2] * (norm[2] / fabs(norm[2])); //this uses the edge os the hole as the position in z
		}else{ //if the position is not in the hole
			holyLoc[2] = norm[2]; // uses the current position of the charge
		}
		add3Arr(closestPoint, holyLoc, locCurr); //renormalizes the origin of the closest position
	}
}

int isAboveHole(double* posUse, MetalSheetClass objecti){
    /* 
	 * FUNCTION: isAboveHole
	 * ------------------------------
	 * 
	 * determines whether or not a position is: above a hole, in a hole, or neither
	 * 
	 * posUse: the position which the will be tested for its 'holiness'
     * objecti: the object (presumable with holes) which the position will be tested against
	 * 
	 * returns: gives an integer whose magnitude is the index (starting at 1) and whose sign tells whether or not
     *      the position is inside the hole (negative) or simply above (positive). 0 means that it is not above
     *      any hole.
	 */

	double* center = (&objecti)->pos;
	double norm[3]; subtract3Arr(norm, posUse, center); //set origin to the center of the object

	int numHole = (&objecti)->numHoles;
	HoleClass* holes = (&objecti)->holes;
	int res = 0; //initialize result (will stay 0 if no hole is encountered)

	int i;

	for(i = 0; i < numHole; i++){ //for each hole in the object
		double* locCurr = (&(holes[i]))->pos; //sets position to a location instantiation
		double differ[3]; subtract3Arr(differ, norm, locCurr); //nomalizes location to be centered on the hole
		differ[2] = 0; //moves problem to 2 dimensions
		double mag = mag3Arr(differ);
		int isAbove = (mag < (&(holes[i]))->R)?1:0; //looks to see if project distance of the position from the center of the hole is within the hole
		int isIn = (fabs(norm[2]) < (0.5 * ((&objecti)->geo)[2]))?1:0; //checks if the position is in the hole
		if(isAbove){ //if the position is, indeed, within the bounds of the radius of the hole 
			res = i + 1; //set the result equal to the index of the hole plus one (since 0 is the no-hole case)
			if(isIn){ //if the position is in the hole
				res = -1 * res; //set the result negative if the position is inside of the hole itself
			}
			break; //break if a hole is found (given the 2-dimension nature of the sheets, a position cannot be above more than one hole)
		}
	}

	return(res);
}

double potentialShot(ElectronClass charge, MetalSheetClass* objects, int numObs, double lBound, double hBound){
    /* 
	 * FUNCTION: potentialShot
	 * ------------------------------
	 * 
     * gives a Monte Carlo walk-on-spheres shot of a given charge in with a set of boundary conditions provided by the objects,
     * using the lower bound lBound as the epsilon value and the upper bound hBound as a large bounding box to protect against
     * a runaway charge. (This is extremely unlikely to happen because all surfaces in an unbound set of boundary conditions are
     * attractive, and thus make it incredibly hard for the charge to climb all the way up to the upper bounding box.)
	 * 
	 * charge: the electric charge being tested in question
     * objects: the objects of the system providing the boundary conditions
     * numObs: the number of objects in the system (needed)
     * lBound: the lower bound under which the simulation will be stopped (often considered the epsilon value)
     * hBound: the upper bound that creates a bounding box to the problem
     * 
	 * returns: the resultant Monte Carlo shot value from the walk-on-spheres method
	 */

	double c = 299792458.0; //speed of light [m/s]
	double u0 = 4.0 * PI * pow(10, (-7)); // vacuum permeability [kg*m/C^2]
	double e0 = 1.0 / ((c * c) * u0); //vacuum permitivity[F/m]
	double posUse[3]; equal3Arr(posUse, (&charge)->pos); //declares posUse in new variable

	while(true){ //this loop will run until it breaks from the inside
		double rPoints[numObs]; //array of the distances from the current position to the closest point on an object
		double rPointsEStats[numObs]; //array of distance from the charge to their closest point on an object
		int i;

		for(i = 0; i < numObs; i++){ //for each object...
			double point[3]; findClosestPoint(point, posUse, objects[i]);
			double diff[3]; subtract3Arr(diff, posUse, point);
			double diffEStat[3]; subtract3Arr(diffEStat, (&charge)->pos, point);
			/*
			diff and diffEStat are extrememly important to keep separate. diff is simply
			how far away the closest point is from the position you're looking at (THIS
			IS ONLY THE CHARGE POSITION ON THE FIRST LOOP.) The difference between the
			closest position and the charge is stored in diffEStat (EStat meaning
			electrostatic here). This helps to keep the distances correct so that
			the correct value is used in the correct application
			*/
			rPoints[i] = mag3Arr(diff);
			rPointsEStats[i] = mag3Arr(diffEStat);
		}
		int mini = 0;
		double r = rPoints[mini];
		for (i = 0; i < numObs; ++i){ //for each object...
		    if (rPoints[i] < r){ //find which object is closest to the position being used
		        r = rPoints[i];
		        mini = i;
		    }
		}
		if((r < lBound) || (r > hBound)){ //if the position terminates or is out of bounds
			double eStatConst = 0.5 * ((&charge)->Q) / (4 * PI * e0);
			return(((&(objects[mini]))->V) - (eStatConst / rPointsEStats[mini])); //return the voltage
		}
		/*
		If the position being used is still larger than the lower bound and within
		the upper bound, the systems finds a new position on the bounding sphere
		and continues the Monte Carlo simulation, a la the walk on spheres method.
		*/
		double phi, theta;
		double irm = 1 / ((double) RAND_MAX);

    	phi = (((double)rand()) * irm) * (2.0 * PI); //generates a random phi
    	theta = acos((2.0 * (((double)rand()) * irm)) - 1.0); //generates random theta (the acos is needed to keep the distribution uniform)
    	double newDir[] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)}; //generates a new unit direction vector
		double newDisp[3]; constMultiply3Arr(newDisp, r, newDir); //multiples unit direction vector by distance to boundary
		add3Arr(posUse, newDisp, posUse); //adds new position component to vector
	}

	return(-1);
}
