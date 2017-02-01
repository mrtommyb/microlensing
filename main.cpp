// Sample main project file.
// This code has been developed by Valerio Bozza, University of Salerno.
// Any use of this code for scientific publications should be acknowledged by a citation to
// V. Bozza, MNRAS 408 (2010) 2188

#include "stdafx.h" //(if you are in Visual Studio, you may require this line

#include<stdio.h>
#include <math.h>
#include"..\\..\\VBBinaryLensingLibrary\\VBBinaryLensingLibrary\\VBBinaryLensingLibrary.h"

int main()
{
//////////////////////////////////////////
	// First declare an instance to the VBBinaryLensing class. 
	// This can be done once in your code and then you can refer to this instance 
	// whenever you need to use the binary magnification calculation
	VBBinaryLensing *VBBL=new VBBinaryLensing;

//////////////////////////////////////////
// Let us calculate the magnification for a specific binary lens and a particular source position
	double Mag,s,q,y1,y2,Rs,accuracy;

	s=0.8; //separation between the two lenses [should be < 6]
	q=0.1; // mass ratio: mass of the lens on the right divided by mass of the lens on the left [should be > 1.e-5]
	
	// Position of the center of the source with respect to the center of mass. [ |y| should be < 10 ]
	y1=0.01; // y1 is the source coordinate along the axis parallel to the line joining the two lenses 
	y2=0.01; // y2 is the source coordinate orthogonal to the first one
	Rs=0.01; // Source radius in Einstein radii of the total mass. [should be < 0.1]

	accuracy=1.e-2; // Absolute accuracy of the result. The result will be: Mag +- accuracy 
	Mag=VBBL->BinaryMag(s,q,y1,y2,Rs,accuracy); // Call to the BinaryMag function with these parameters
	printf("Magnification = %lf\n",Mag); // Output should be 18.28....

//////////////////////////////////////////
// Let us calculate the magnification for a linear limb-darkened source
	double a1;
	a1=0.51; // Linear limb darkening coefficient
	Mag=VBBL->BinaryMagDark(s,q,y1,y2,Rs,a1,accuracy); // Call to the BinaryMagDark function with these parameters
	printf("Magnification with limb darkened source = %lf\n",Mag);  // Output should be 18.27.....

//////////////////////////////////////////
// Several functions are also available in the form VBBL->function(double *pr,double t);
// where pr is an array of parameters and t is the time where we want to calculate the function.
// Note that the accuracy of these functions is specified separately acting on VBBL->Tol
// By default these function are without limb darkening.

	double pr[15];
	double u0,t0,tE,alpha,t;

	u0=0.09; // Impact parameter
	alpha=0.53; // Angle between a vector pointing to the left and the source velocity
	t0=7550.4; // Time of closest approach to the center of mass
	tE=100.3; // Einstein time

	// Putting all parameters in the array
	pr[0]=log(s); // Note that when fitting it is sometimes convenient to fit the log, as in this case
	pr[1]=log(q);
	pr[2]=u0;
	pr[3]=alpha;
	pr[4]=log(Rs);
	pr[5]=log(tE);
	pr[6]=t0;
	VBBL->Tol=1.e-2; // Setting accuracy 
	t=7561.2; 
	Mag=VBBL->BinaryLightCurve(pr,t);
	printf("Another magnification calculation using BinaryLightCurve = %lf\n",Mag);  // Output should be 24.33.....

//////////////////////////////////////////
// Parallax calculation
// First you need to initialize the event coordinates

	char coordinatefile[256]="OB151212coords.txt"; // Text file containing the event coordinates in J2000.0
	// The format should be HH:MM:SS.SSS +DD:MM:SS.SSS (see the sample file provided)
	char sattabledir[256]="."; // Directory where satellite positions tables lie
	// Only important if you are using spacecraft observations (see below).
	double year=2015; // Year of the event (needed for equinox precession calculation)

	VBBL->SetObjectCoordinates(coordinatefile,sattabledir); 
	// This function sets the event coordinates in the VBBL library and 
	// loads the satellite position tables (if present).
	// You should call it only at the begginning of your analysis or when you move to another event.

	double pai1,pai2;
	pai1=0.3; // Parallax component parallel to the Earth acceleration (let us call it \alpha).
	pai2=0.13; // Parallax component orthogonal to the Earth acceleration (directed toward \alpha \wedge Object direction)

	pr[7]=pai1; // Include these two additional parameters in the parameter array
	pr[8]=pai2;

	// let us calculate the value at a time far from t0, just to se the effect:

	t=7651.2;

	Mag=VBBL->BinaryLightCurve(pr,t);
	printf("Using BinaryLightCurve = %lf\n",Mag); // Output should be 1.39....

	Mag=VBBL->BinaryLightCurveParallax(pr,t);
	printf("Using BinaryLightCurveParallax = %lf\n",Mag); // Output should be 1.27....

	// In alternative, you can use the parallax North and East components.
	// If you prefer this coordinate system, set
	VBBL->parallaxsystem=1;
	// Now pr[7] and pr[8] are the North and East components respectively.
	Mag=VBBL->BinaryLightCurveParallax(pr,t);
	printf("Using BinaryLightCurveParallax = %lf\n",Mag); // Output should be 1.30....


//////////////////////////////////////////
// Calculation of magnification as seen by a satellite
// The satellite position table should be in the format generated by http://ssd.jpl.nasa.gov/horizons.cgi
// In particular, we assume five columns:

// JD
// RA (degrees)
// Dec (degrees)
// Distance from Earth (AU)
// Distance rate change (not really needed but included by default in Horizons).

// See the satellite tables attached as examples.
// The table file names should be "satellite*.txt" (with * replaced by a single character). 
// These tables are sorted alphabetically and assigned a satellite number.

// If you want the magnification as seen from satellite 1, then just set VBBL->satellite to 1 before the calculation.

	VBBL->satellite=1;
	t=7561.2;
	Mag=VBBL->BinaryLightCurveParallax(pr,t);
	printf("Using BinaryLightCurveParallax for satellite 1: %lf\n",Mag); // Output should be 2.35....

// If you want to return to the ground do not forget to set VBBL->satellite back to 0
	VBBL->satellite=0;


//////////////////////////////////////////
// Other functions also implemented in VBBL:

// PSPLCurve(double *pr,t)						(Parameters are log_u0, log_tE, t0)
// PSPLParallaxCurve(double *pr,t)				(u0, log_tE, t0, pai1, pai2)
// ESPLCurve(double *pr,t)						(log_u0, log_tE, t0, logRs)
// ESPLParallaxCurve(double *pr,t)				(u0, log_tE, t0, logRs, pai1, pai2)
// BinaryLightCurveOrbital(double *pr,t)		(log_s, log_q, u0, alpha_0, log_Rs, log_tE, t0, pai1, pai2, w1, w2, w3)
// Orbital parameters are in the hypothesis of circular motion (no eccentricity)
// w1=(ds/dt)/s
// w2=dalpha/dt
// w3=(dsz/dt)/s

// BinSourceMag(double *pr,t)					(log_tE, log_FluxRatio, u0_1, u0_2, t0_1, t0_2)
// BinSourceParallaxMag(double *pr,t)			(log_tE, log_FluxRatio, u0_1, u0_2, t0_1, t0_2, pai1, pai2)
// BinSourceXallarapMag(double *pr,t)			(log_tE, log_FluxRatio, u0_1, u0_2, t0_1, t0_2, pai1, pai2, w1, w2, w3)

//////////////////////////////////////////
	delete VBBL; // At the end of the program you should delete the instance to the VBBinaryLensing class you used.

	getchar();

	return 0;
}

