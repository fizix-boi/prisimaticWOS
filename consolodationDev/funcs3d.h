#include <stdio.h>
#include <math.h>

//function declarations

void print3Arr(double*);
void add3Arr(double*, double*, double*);
void subtract3Arr(double*, double*, double*);
void multiply3Arr(double*, double*, double*);
void constMultiply3Arr(double*, double, double*);
void sign3Arr(double*, double*);
void equal3Arr(double*, double*);
double mag3Arr(double*);

//functions

void print3Arr(double* threeArr){
	/* 
	 * FUNCTION: print3Arr
	 * ------------------------------
	 * 
	 * prints a 3-array in a format of a row vector
	 * 
	 * threeArr: the 3-array to be printed by the function
	 * 
	 * returns: void
	 */

	printf("(%e, %e, %e)\n", threeArr[0], threeArr[1], threeArr[2]);
}

void add3Arr(double* res, double* arr1, double* arr2){
	/* 
	 * FUNCTION: add3Arr
	 * ------------------------------
	 * 
	 * adds two arrays and sets the result equal to another given result array
	 * note that this function is commutative in arr1 and arr2
	 * 
	 * res: the resultant array which the sum will be cast to
	 * arr1: the first array that will be in the sum
	 * arr2: the second array in the sum 
	 * 
	 * returns: void
	 */

	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = arr1[i] + arr2[i];
	}
}

void subtract3Arr(double* res, double* arr1, double* arr2){
	/* 
	 * FUNCTION: subtract3Arr
	 * ------------------------------
	 * 
 	 * subtracts two arrays and sets the result equal to another given result array
	 * note that this function is NOT commutative in arr1 and arr2
	 * 
	 * res: the resultant array which the difference will be cast to
	 * arr1: the first array that will be in the difference
	 * arr2: the second array in the difference 
	 * 
	 * returns: void
	 */
	
	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = arr1[i] - arr2[i];
	}
}

void multiply3Arr(double* res, double* arr1, double* arr2){	
	/* 
	 * FUNCTION: multiply3Arr
	 * ------------------------------
	 * 
	 * multiplies two arrays and sets the result equal to another given result array
	 * note that this function is commutative in arr1 and arr2
	 * 
	 * res: the resultant array which the product will be cast to
	 * arr1: the first array that will be in the product
	 * arr2: the second array in the product 
	 * 
	 * returns: void
	 */

	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = arr1[i] * arr2[i];
	}
}

void constMultiply3Arr(double* res, double num, double* arr){
	/* 
	 * FUNCTION: constMultiply3Arr
	 * ------------------------------
	 * 
	 * multiplies a 3-array by a constant value and sets the result equal to another given result array
	 * 
	 * res: the resultant array which the product will be cast to
	 * num: the constant which the array will be multiplued by
	 * arr: the array that will have a constant multipied to it
	 * 
	 * returns: void
	 */

	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = num * arr[i];
	}
}

void sign3Arr(double* res, double* arr){
	/* 
	 * FUNCTION: sign3Arr
	 * ------------------------------
	 * 
	 * gives the sign of each of the elements in a 3-array and casts the result to a different array
	 * 
	 * res: the resulting sign 3-array
	 * arr: the 3-array of which the sign function will be taken
	 * 
	 * returns: void
	 */

	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = (arr[i] > 0) ? 1.0 : -1.0;
	}
}

void equal3Arr(double* res, double* arr){
	/* 
	 * FUNCTION: equal3Arr
	 * ------------------------------
	 * 
	 * casts a 3-array to another 3-array
	 * 
	 * res: 3-array which will take in the values of the given 3-arr
	 * arr: given 3-arr to cast to the resultant 3-array
	 * 
	 * returns: void
	 */

	int i;

	for(i = 0; i <= 2; ++i){
		res[i] = arr[i];
	}
}

double mag3Arr(double* arr){
	/* 
	 * FUNCTION: mag3Arr
	 * ------------------------------
	 * 
	 * gives the magnitude of a 3-array
	 * 
	 * arr: the 3-array of which the magnitude will be given
	 * 
	 * returns: magnitude of the given vector as a double
	 */

	return(sqrt((arr[0] * arr[0]) + (arr[1] * arr[1]) + (arr[2] * arr[2])));
}