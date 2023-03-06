//primative object declarations

typedef struct {
    /* 
     * STRUCT: ElectronClass
     * ------------------------------
     * 
     * Models an electric charge (colloquially labeled as 'electron') at some position
     * 
     * Q: the charge of the electric charge
     * pos: the position as a 3-array of the electric charge
     */

	double Q;
	double* pos;
} ElectronClass;

void electron(ElectronClass *elec, double Q, double* pos){
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

	elec->Q = Q;
	elec->pos = pos;
}

typedef struct {
    /* 
     * STRUCT: HoleClass
     * ------------------------------
     * 
     * Models a hole that is in a metal sheet or object
     * 
     * pos: the position as a 3-array of the hole
     * R: the radius of the hole
     */

	double* pos;
	double R;
} HoleClass;

typedef struct {
    /* 
     * STRUCT: MetalSheetClass
     * ------------------------------
     * 
     * Models a rectangular prisimatic piece of a perfectly conducting metal
     * 
     * V: the voltage of the sheet (acts as the equipotential boundary conditions for the Laplace enviornment)
     * pos: the position as a 3-array of the center of the sheet
     * geo: a 3-array with the length of the prism in x, y, and z
     * numHoles: the number of holes in the metal objects
     * holes: the array of hole objects that give the holes in the metal sheet
     */

	double V;
	double* pos;
	double* geo;
	int numHoles;
	HoleClass* holes;
} MetalSheetClass;

void hole(HoleClass *hol, double* pos, double R){
    /* 
	 * FUNCTION: hole
	 * ------------------------------
	 * 
	 * creates a hole of radius R at position pos using the HoleClass struct
	 * 
	 * hol: the hole object that the parameters will be set to
     * R: the radius of the hole
     * pos: the position as a 3-array of the hole
	 * 
	 * returns: void
	 */

	hol->pos = pos;
	hol->R = R;
}

void metalSheet(MetalSheetClass *meta, double V, double* pos, double* geo, int numHoles, HoleClass* holes){
    /* 
	 * FUNCTION: metalSheet
	 * ------------------------------
	 * 
	 * creates a metal sheet with voltage V, position p, geometry geo, and hole array holes (with numHoles number of holes) using the MetalSheetClass struct
	 * 
	 * meta: the metal sheet object that the parameters will be set to
     * V: the voltage the metal sheet will have
     * pos: the position of the center of the metal sheet
     * geo: the geometry of the metal sheet
     * numHoles: the number of holes of the metal sheet
     * holes: the holes in question for the metal sheet
	 * 
	 * returns: void
	 */

	meta->V = V;
	meta->pos = pos;
	meta->geo = geo;
	meta->numHoles = numHoles;
	meta->holes = holes;
}