#ifndef __MESH_H
#define __MESH_H

#define dfloat double
#define dint int

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"
#include <sys/stat.h>

#define mymax(a,b) ( ((a)>(b)) ? (a):(b) )
#define mymin(a,b) ( ((a)<(b)) ? (a):(b) )

typedef struct{

	clock_t CPUstarttime;
	clock_t CPUelapsedtime;
	dfloat CPUelapsedtimeInSeconds;

	//------------------------------------------------------------------!
	// Declarations for meshReadStatistics
	//------------------------------------------------------------------!

	dint ncellsx; /* number of cells along x-direction*/
	dint ncellsy; /* number of cells along y-direction*/
	dint ncellsz; /* number of cells along z-direction*/

	dint ncellVertices; /* number of vertices per cell */
	dint ncellFaces;    /* number of faces per cell */
	dint nfaceVertices; /* number of vertices per Face */
	dint ndimension;    /*  3 =>  x,y,z */

	dint nicells; /* number of inner cells */
	dint nbcells; /* number of boundary cells */
	dint ncells;  /* number of cells in mesh */
	dint nivertices; /* number of vertices of inner cells in mesh */
	dint nvertices;  /* number of vertices in mesh */

	dint dim_x;
	dint dim_y;
	dint dim_z;

	//------------------------------------------------------------------!
	//  Declarations for meshSettings
	//------------------------------------------------------------------!
	dfloat SMALL;
	dfloat VSMALL;
	dfloat VVSMALL;

	dfloat millisecond;
	dfloat microsecond;
	dfloat nanosecond;
	dfloat picosecond;
	dfloat femtosecond;

	dfloat minute;
	dfloat hour;
	dfloat day;

	dfloat PI;
	dfloat permeabilityofFreeSpace;

	dint vertA;
	dint vertB;
	dint vertC;
	dint vertD;
	dint vertE;
	dint vertF;
	dint vertG;
	dint vertH;

	dint faceAEFB; // leftFace
	dint faceEGHF; // frontFace
	dint faceGCDH; // rightFace
	dint faceCABD; // backFace
	dint faceFHDB; // topFace
	dint faceACGE; // bottomFace

	dint leftFace;
	dint frontFace;
	dint rightFace;
	dint backFace;
	dint topFace;
	dint bottomFace;


	//------------------------------------------------------------------!
	// Declarations for meshConstructVerticesFields
	//------------------------------------------------------------------!

	dfloat *vertices; /* coordinates of vertices */

	char *variableName ;
	char *directoryName ;
	char *subdirectoryName ;
	char *fileName;

	//------------------------------------------------------------------!
	//  Declarations for meshConstructConnectivityFields
	//------------------------------------------------------------------!

	dint *cellToVertex; /* Global cell to Global Vertex connectivity */
	dint *cellToCells; /* Global cell to Global cell connectivity -1 for boundaryface*/
	dint *cellToFaces; /* Global cell to adjacent Global cell Face connectivity */
	dint *faceToVertices; /* local face to Local vertices connectivity */


	//------------------------------------------------------------------!
	//  Declarations for meshCalculateNoOfBoundaryandInnerFaces
	//------------------------------------------------------------------!

	dint nInnerFaces; /* No of inner faces */
	dint nBoundaryFaces; /* No of boundary faces */

	//------------------------------------------------------------------!
	//  Declarations for meshConstructFaceToCellConnectivityFields
	//------------------------------------------------------------------!

	dint *innerFaces; /* Cell Id of inner faces */
	dint *boundaryFaces; /* Cell Id of boundary faces */

	dint *faceToOwnerCell;    /* Face to Owner Cell connectivity */
	dint *faceToNeighborCell; /* Face to Neighbor Cell connectivity */


	//------------------------------------------------------------------!
	//  Declarations for meshConstructMetricsFields
	//------------------------------------------------------------------!

	dfloat *C; /* coordinates of cell centers */
	dfloat *Cf; /* coordinates of cell face centers */
	dfloat *Sf; /* cell normal face Areas - outward normal*/
	dfloat *magSf; /* Magnitude of Face areas */
	dfloat *Sn; /* Surface normal Unit vector */

	dfloat *V; /* cell volumes (in 2D Cell Area) */
	dfloat *rV; /* reciprocal of cell volumes */

	dfloat *d; /* cell center to cell center distance , Vector from from cell P to cell N  */
	dfloat *magd; /* distance from P to N */
	dfloat *rmagd; /* reciprocal of magd */

	dfloat *df;     /* Cell to Face distance Vector  */
	dfloat *magdf;  /* Cell to Face distance */
	dfloat *rmagdf; /* reciprocal of magdf */

	dfloat *delta;     /* Combining magd and magdf for inner and boundary cells */
	dfloat *rdelta;     /* Combining rmagd and rmagdf for inner and boundary cells */

	dfloat *fx; /* Cell Face interpolation Factor */

	dint *boundaryFacePatchID; /* Name each boundary face as inlet, outlet, top wall, bottom wall */

	//------------------------------------------------------------------!
	//  Declarations for meshCountNoOfBoundaryPatches
	//------------------------------------------------------------------!

	dint  npatchs; /* number of patches */

	//------------------------------------------------------------------!
	//  Declarations for meshConstructPatchConnectivityFields
	//------------------------------------------------------------------!

	dint *patchID;
	dint *npatchFaces; /* array of size npatches with entries of no of boundary Faces in each patch */
	dint *patchStartFaceIndex; /* start Face Index of Each Patch */
	dint *patchToBoundaryFace; /* Patch to Boundary Face Connectivity */

	dint *patchType;

	dint *periodicNeighborPatchID; /* Periodic Neighbor Patch */
	dint *periodicNeighborPatch; /* Periodic Neighbor Patch */

	dint *faceToBoundaryCell; /* Face to Boundary Cell connectivity */

	//------------------------------------------------------------------!
	//  Declarations for meshGeneratePeriodicConnectivity
	//------------------------------------------------------------------!

	dint nInnerAndPeriodicFaces; /* No of inner faces */

	//------------------------------------------------------------------!
	//  Declarations for meshConstructInnerAndPeriodicFacesFields
	//------------------------------------------------------------------!

	dint *innerAndPeriodicFaces; /* Face Id of inner and Periodic faces */

	//------------------------------------------------------------------!
	//  Declarations for meshConstructDomainFields
	//------------------------------------------------------------------!

	dint *domain;


	//------------------------------------------------------------------!
	//  Declarations for meshConstructReconstructFields
	//------------------------------------------------------------------!
	dfloat *SnDyadSf;
	dfloat *summedSnDyadSf;
	dfloat *inverseOfSummedSnDyadSf;


	//------------------------------------------------------------------!
	//  Declarations for meshConstructConstantsFields
	//------------------------------------------------------------------!

	dfloat *g;     		// acceleration due to gravity
	dfloat *gf;    		// acceleration due to gravity, face values

	//------------------------------------------------------------------!
	//  Declarations for setConstants
	//------------------------------------------------------------------!

	dfloat rhoG;      	// Density of Gas
	dfloat rhoL;      	// Density of Liquid

	dfloat muG;      	// Dynamic viscosity of Gas
	dfloat muL;      	// Dynamic viscosity of Liquid

	dfloat Re;      	// Reynolds number
	dfloat Fr;     	    // Froude number
	dfloat Ha;      	// Reynolds number

	dfloat alphaG;      // Volume fraction of Gas
	dfloat alphaL;      // Volume fraction of Liquid

	dfloat sigmaG;      // conductivity of Gas
	dfloat sigmaL;      // conductivity of Liquid
	//------------------------------------------------------------------!
	//  Declarations for setScaleFactors
	//------------------------------------------------------------------!

	dfloat L_0; 	// Length Scale
	dfloat rho_0;   // Density scale
	dfloat g_0; 	// Acceleration due to gravity Scale
	dfloat U_0;     // Velocity Scale
	dfloat t_0;     // Time Scale
	dfloat p_0; 	// Pressure Scale
	dfloat mu_0;    // viscosity scale

	dfloat sigma_0; // Electrical conductivity Scale
	dfloat Phi_0;   // Electric Potential Scale
	dfloat B_0; 	// Magnetic Field scale
	dfloat J_0; 	// Current Density scale

	//------------------------------------------------------------------!
	//  Declarations for meshConstructGravityMetricsFields
	//------------------------------------------------------------------!
	dfloat *Xref;   // reference height
	dfloat *X;      // Cell center minus reference height
	dfloat *Xf;     // Cell face center Height minus reference height
	dfloat *gDotX;  // acceleration due to gravity dot product with position vector
	dfloat *gDotXf; // acceleration due to gravity dot product with position vector

	dfloat *dh;       // differential reference height
	dfloat *gDotdh;

	//------------------------------------------------------------------!
	//  Declarations for constructSolverFields
	//------------------------------------------------------------------!

	// For Scalar Equation
	dfloat *aP; // diagonal coefficients for scalar matrix eqn
	dfloat *aN; // neighbor cell coefficients for scalar matrix eqn
	dfloat *bP; // source term for scalar matrix eqn

	// For Vector Equation
	dfloat *AP; // diagonal coefficients for vector matrix eqn
	dfloat *AN; // neighbor cell coefficients for vector matrix eqn
	dfloat *BP; // source term for vector matrix eqn

	// For Computing HbyA, store coefficients without grad P term for Uequation
	dfloat *AAP; // diagonal coefficients for vector matrix eqn
	dfloat *AAN; // neighbor cell coefficients for vector matrix eqn
	dfloat *BBP; // source term for vector matrix eqn

	//------------------------------------------------------------------!
	//  Declarations for solverSetControls
	//------------------------------------------------------------------!

	dint starttimeStep;
	dint timeStepMax;
	dfloat finalTime;
	dfloat deltatime;
	dint WriteFrequency;

	bool timeStepSizeControl;
	dfloat Cot; // target courant number: user input
	dfloat Kl;
	dfloat lambdat;

	dint 	divScheme;
	dfloat 	divSchemeblendFactor;

	dint   targetNCorrector;

	dfloat 	UEqnlambda;  	 // Under relaxation parameter
	dfloat 	PEqnlambda;  	 // Under relaxation parameter
	dfloat 	alphaEqnlambda;  // Under relaxation parameter



	bool includeGravity;
	bool includeLorentzForce;

	//------------------------------------------------------------------!
	//  Declarations for mainCompressibleNavierStokesSolve
	//------------------------------------------------------------------!
	dint timeStep;
	dfloat time;
	bool stopCriterion;
	dint currentNCorrector; // counter for no of outer corrector

	//------------------------------------------------------------------!
	//  Declarations for constructSolutionFields
	//------------------------------------------------------------------!

	dfloat *U;	   		// Velocity
	dfloat *p;     		// Pressure
	dfloat *alpha; 	    // volume fraction

	dfloat *rho;  	  	// Density
	dfloat *mu;  	  	// Dynamic viscosity

	dfloat *rhof;  	  	// Density face values
	dfloat *muf;  	  	// Dynamic viscosity face values

	dfloat *phiU;    	// Velocity Flux stored at cell faces
	dfloat *rhophiU;    // Mass Flux stored at cell faces

	dfloat *Cof; 		// Courant number at cell faces
	dfloat *Co;  		// Courant number at cell center: Maximum of cell face values

	dfloat *Ustar;      //  Reference Velocity
	dfloat *pstar;      //  Reference Pressure

	dfloat *rAP;     	// reciprocal of diagonal coefficients of momentum eqn matrix
	dfloat *rAPf;    	// face values

	dfloat *HbyA;    	// derived parameter from discretized momentum equation
	dfloat *phiHbyA; 	// Flux of HbyA

	dfloat *ph;     		// hydrostatic Pressure

	//Pressure gradient variables
	dfloat *gradPbyA;
	dfloat *phigradPbyA;


	// Body Forces
	dfloat *BF1;		// Body Force rho*g term
	dfloat *BF2;		// Body Force hydrostatic pressure gradient
	dfloat *BF3;		// Body Force due to Lorentz force




	dfloat *phi;     	// Electric Potential

	dfloat *B;		 	// Applied External Magnetic field Intensity
	dfloat *Bf;	 		// Applied External Magnetic field face values
	dfloat *phiB;    	// Applied Magnetic field flux stored at cell faces

	dfloat *sigma;     // conductivity
	dfloat *sigmaf;	   // conductivity face values

	dfloat *sigmaUcrossBFlux;	 // Flux of Sigma times U cross B
	dfloat *sigmagradPhiFlux;	 // Flux of Sigma times Gradient of Potential
	dfloat *phiJ;    			 // Total current Density flux stored at cell faces

	dfloat *Jphi;    			// Current density due to electric potential
	dfloat *Ju;		 			// Current density due to motion of conducting liquid
	dfloat *J;		 			// Total Current Density

	dfloat *A;		 			// Magnetic Vector Potential due to applied current
	dfloat *B_J;		 	    // Magnetic Field due to applied current


	//------------------------------------------------------------------!
	//  Declarations for computePostProcessingParameters
	//------------------------------------------------------------------!

	dfloat iMass;

	dfloat iKE;

	dfloat iKEx;
	dfloat iKEy;
	dfloat iKEz;

	//------------------------------------------------------------------!
	//  Declarations for computeDeltaTime
	//------------------------------------------------------------------!

	dfloat CoU;
	dfloat Co0; // Courant Number from Previous time step
	dfloat deltatt;

	//------------------------------------------------------------------!
	//  Declarations for computeResidualsScalarEqn
	//------------------------------------------------------------------!

	dfloat Residual_max;
	dfloat Residual_rms;
	dfloat Residual_max_scaled;
	dfloat Residual_rms_scaled;

	//------------------------------------------------------------------!
	//  Declarations for computeResidualsVectorEqn
	//------------------------------------------------------------------!

	dfloat Residual_max_Xeqn;
	dfloat Residual_max_Yeqn;
	dfloat Residual_max_Zeqn;

	dfloat Residual_rms_Xeqn;
	dfloat Residual_rms_Yeqn;
	dfloat Residual_rms_Zeqn;

	dfloat Residual_max_scaled_Xeqn;
	dfloat Residual_max_scaled_Yeqn;
	dfloat Residual_max_scaled_Zeqn;

	dfloat Residual_rms_scaled_Xeqn;
	dfloat Residual_rms_scaled_Yeqn;
	dfloat Residual_rms_scaled_Zeqn;

	//------------------------------------------------------------------!
	//  Declarations for computeResidualUEquation
	//------------------------------------------------------------------!
	dfloat UEqnResidual_max_Xeqn;
	dfloat UEqnResidual_max_Yeqn;
	dfloat UEqnResidual_max_Zeqn;

	dfloat UEqnResidual_rms_Xeqn;
	dfloat UEqnResidual_rms_Yeqn;
	dfloat UEqnResidual_rms_Zeqn;

	dfloat UEqnResidual_max_scaled_Xeqn;
	dfloat UEqnResidual_max_scaled_Yeqn;
	dfloat UEqnResidual_max_scaled_Zeqn;

	dfloat UEqnResidual_rms_scaled_Xeqn;
	dfloat UEqnResidual_rms_scaled_Yeqn;
	dfloat UEqnResidual_rms_scaled_Zeqn;


	//------------------------------------------------------------------!
	//  Declarations for computeResidualPEquation
	//------------------------------------------------------------------!
	dfloat PEqnResidual_max;
	dfloat PEqnResidual_rms;
	dfloat PEqnResidual_max_scaled;
	dfloat PEqnResidual_rms_scaled;

	//------------------------------------------------------------------!
	//  Declarations for computeResidualPhiEquation
	//------------------------------------------------------------------!
	dfloat PhiEqnResidual_max;
	dfloat PhiEqnResidual_rms;
	dfloat PhiEqnResidual_max_scaled;
	dfloat PhiEqnResidual_rms_scaled;

}mesh_t;

typedef enum {
	vertexScalar		 =	0,
	vertexVector		 =	1,
	volumeScalar		 =	2,
	volumeVector		 =	3,
	surfaceScalar		 =	4,
	surfaceVector		 =	5,
	scalar			 	 =	6,
	vector			 	 =	7,
	cellFaceVertexScalar =	8,
	cellFaceVertexVector =	9,
	charField 			 =	10,
	volumeScalar2		 =	11,
	volumeScalar3		 =	12,
	volumeVector2		 =	13,
	volumeVector3		 =	14,
	timeProbeScalar      =  15,
	timeProbeVector      =  16,
	volumeTensor		 =  17,
	surfaceTensor		 =  18,
} field_type;

typedef enum {
	emptyPatch   		 =	0,
	boundaryPatch		 =	1,
	periodicPatch		 =	2,
} patch_type;

typedef enum {
	initialTimeStep 	=	0,
	targetCourantNumber =	1,
} timeStep_control;

typedef enum {
	linear				 =	0,
	upwind  			 =	1,
	interfaceCompression =  2,
	Gamma  				 =	3,
	VanLeer  		     =	4,
	MUSCL				 =  5,
	Minmod				 =  6,
	SuperBee			 =  7
} interpolation_schemes;

typedef enum {
	euler		    =	0,
	eulerBackward  	=	1,
} time_schemes;

typedef enum {
	fixedValue			=	0,
	zeroNormalGradient  =	1,
	periodic			=	2,
	fixedNormalGradient	=	3,
	totalPressureInletPBC    =	4,
	totalPressureOutletPBC    =	5,
	emptyBoundary       =	6,
} bc_type;



// meshmetrics

void meshReadStatistics(mesh_t *mesh);
void meshSettings(mesh_t *mesh);

void createField(dfloat** phi, dint variableType, mesh_t *mesh);
void createCharField(char** phi, mesh_t *mesh);
void createIntField(dint** phi, dint variableType, mesh_t *mesh);


void meshConstructVerticesFields(mesh_t *mesh);
void meshDestructVerticesFields(mesh_t *mesh);

void meshReadVertices(mesh_t *mesh);

void printVertexVector(dfloat* Phi, mesh_t *mesh);
dint computeLinearIndex(dint i, dint j, dint k, dint nx, dint ny, dint nz);

void meshConstructConnectivityFields(mesh_t *mesh);
void meshDestructConnectivityFields(mesh_t *mesh);

void meshGenerateCellConnectivity(mesh_t *mesh);

void meshCalculateNoOfBoundaryandInnerFaces(mesh_t *mesh);

void meshConstructFaceToCellConnectivityFields(mesh_t *mesh);
void meshDestructFaceToCellConnectivityFields(mesh_t *mesh);

void meshGenerateFaceToCellConnectivity(mesh_t *mesh);

void meshConstructMetricsFields(mesh_t *mesh);
void meshDestructMetricsFields(mesh_t *mesh);

void meshGenerateCellMetrics(mesh_t *mesh);

void printVolumeVector(dfloat* Phi, mesh_t *mesh);
void printVolumeScalar(dfloat* Phi, mesh_t *mesh);
void printSurfaceVector(dfloat* Phi, mesh_t *mesh);
void printSurfaceScalar(dfloat* Phi, mesh_t *mesh);


void meshCountNoOfBoundaryPatches(mesh_t *mesh);

void meshConstructPatchConnectivityFields(mesh_t *mesh);

void meshDestructPatchConnectivityFields(mesh_t *mesh);

void meshGenerateBoundaryPatchConnectivity(mesh_t *mesh);

void meshSetUpPatchType(mesh_t *mesh);

void meshGeneratePeriodicConnectivity(mesh_t *mesh);

void meshConstructInnerAndPeriodicFacesFields(mesh_t *mesh);

void meshDestructInnerAndPeriodicFacesFields(mesh_t *mesh);

void meshCalculateInnerFaceIDIncludingPeriodicFaces(mesh_t *mesh);

void meshConstructDomainFields(mesh_t *mesh);

void meshDestructDomainFields(mesh_t *mesh);

void printVolumeScalarInt(dint* Phi, mesh_t *mesh);

void writeMesh(mesh_t *mesh);

dint getNdimension(dint fieldType);

dint getNdatapercell(dint fieldType);

void writeMeshData(dfloat* Phi, int fieldType,mesh_t *mesh);
void printMeshStencil(mesh_t *mesh);

dint computeLocalFace(dint fid, mesh_t *mesh);

// SetUpCase
const char* getPatchID(dint variableName);
const char* getPatchType(dint variableName);
void meshGenerateBoundaryPatchID(mesh_t *mesh);
void meshSetUpPeriodicNeighborPatches(mesh_t *mesh);
void setUpDomains(mesh_t *mesh);
void setConstants(mesh_t *mesh);
void setScaleFactors(mesh_t *mesh);
void setReferenceHeight(mesh_t *mesh);
void solverSetControls(mesh_t *mesh);
void setInitialConditions(mesh_t *mesh);
void computePostProcessingParameters(mesh_t *mesh);
void computeAndWriteAnalyticalSolution(mesh_t *mesh);
dint getBoundaryConditionType(dint patch, dint phiName, mesh_t *mesh);
void getBoundaryConditionValue(dint patch, dint phiName,  dfloat* bcValue, mesh_t *mesh);
dint getFieldType(dint phiName);
void computeInitialVelocityFlux(mesh_t *mesh);



// mathFunctions
dfloat magnitude(dfloat a);
dfloat magVector(dfloat *X);
dfloat stabilise(dfloat a, mesh_t *mesh);
void copyDataToFrom(dfloat* copyTo,dfloat* copyFrom,int fieldType,mesh_t *mesh);
dfloat dotProduct(dfloat *A,dfloat *B);
dfloat sgn(dfloat F);
dfloat sqr(dfloat phi);
void reinitializeField(dfloat* phi,int fieldType,mesh_t *mesh);



//fvs
void meshConstructReconstructFields(mesh_t *mesh);
void meshDestructReconstructFields(mesh_t *mesh);
void generateMetricsForFluxReconstruction(mesh_t *mesh);
void meshConstructConstantsFields(mesh_t *mesh);
void meshDestructConstantsFields(mesh_t *mesh);
void meshConstructGravityMetricsFields(mesh_t *mesh);
void meshDestructGravityMetricsFields(mesh_t *mesh);
void meshGenerateCellMetricsGravity(mesh_t *mesh);
void constructSolverFields(mesh_t *mesh);
void destructSolverFields(mesh_t *mesh);
void constructSolutionFields(mesh_t *mesh);
void destructSolutionFields(mesh_t *mesh);
void computeDensityFaceValues(mesh_t *mesh);
void computeDynamicViscosityFaceValues(mesh_t *mesh);
void storeOldTimeValues(mesh_t *mesh);
void writeResults(mesh_t *mesh);
void writeResults(mesh_t *mesh);
void writeResultsToFile(dfloat* phi, int fieldType, dfloat lengthScaleFactor, dfloat fieldScaleFactor,mesh_t *mesh);
void writeTimeProbe(dfloat *phi,int fieldType,mesh_t *mesh);
void writeXLineDataToFile(dfloat* Phi, int fieldType, dint K, dfloat lengthScaleFactor, dfloat fieldScaleFactor,mesh_t *mesh);
void writeZLineDataToFile(dfloat* Phi, int fieldType, dint I, dfloat lengthScaleFactor, dfloat fieldScaleFactor,mesh_t *mesh);
void computeDeltaTime(mesh_t *mesh);
void computeCourantNumberUsingConvectiveVelocity(mesh_t *mesh);
void computeDerivedFieldsAlphaEquation(dfloat* alphaDivU,mesh_t *mesh);
int PetscMatrixAssemblyAndSolveScalarEqn(dfloat* alpha,mesh_t *mesh);
int PetscMatrixAssemblyAndSolveVectorEqn(dfloat* alpha,mesh_t *mesh);
void updateBoundaryValues(dfloat* phi, int phiName, mesh_t *mesh);
void constructFluxofSurfaceVector(dfloat* PhiFace, dfloat* Flux, mesh_t *mesh);
void solverUpdateStoppingCriterion(mesh_t *mesh);
void storeReferenceVelocity(mesh_t *mesh);
void storeReferencePressure(mesh_t *mesh);
void constructUEquation(mesh_t *mesh);
void storeUEqnMatrixCoefficients(mesh_t *mesh);
void relaxUEqn(mesh_t *mesh);
void relaxPEqn(mesh_t *mesh);
void computeResidualsScalarEqn(dfloat* phi,dfloat* Residual,mesh_t *mesh);
void computeResidualsVectorEqn(dfloat* phi,dfloat* Residual,mesh_t *mesh);
void computeResidualUEquation(mesh_t *mesh);
void constructrAP(mesh_t *mesh);
void constructHbyA(mesh_t *mesh);
void extrapolateToBoundaryCellScalar(dfloat* phiCell,mesh_t *mesh);
void extrapolateToBoundaryCellVector(dfloat* phiCell,mesh_t *mesh);
void correctExtrapolationToAxisymmetricBoundaryCellScalar(dfloat* phiCell,mesh_t *mesh);
void correctExtrapolationToAxisymmetricBoundaryCellVector(dfloat* phiCell,mesh_t *mesh);
void constructPEquation(mesh_t *mesh);
void computeModifiedFaceArea(dfloat* alpha_f, dfloat* Sfm,dfloat* magSfm,mesh_t *mesh);
void computeResidualPEquation(mesh_t *mesh);
void constructGradPbyA(mesh_t *mesh);
void updateMixtureDensity(mesh_t *mesh);
void updateMixtureViscosity(mesh_t *mesh);
void computeBodyForces(mesh_t *mesh);
void reconstructFluxFromCellFaceToCellCenter(dfloat* Flux, dfloat* cellCenterVector,mesh_t *mesh);
void constructBFbyAFlux(mesh_t *mesh);
void setReferencePressure(mesh_t *mesh);
void updateHydroStaticPressure(mesh_t *mesh);
void setUpExternalMagneticField(mesh_t *mesh);
void computeExternalMagneticFieldFlux(mesh_t *mesh);
void constructPhiEquation(mesh_t *mesh);
void updateMixtureConductivity(mesh_t *mesh);
void computeConductivityFaceValues(mesh_t *mesh);
void computeResidualPhiEquation(mesh_t *mesh);
void computeCurrentDensityFlux(mesh_t *mesh);
void computeCurrentDensityAtCellCenter(mesh_t *mesh);


//fvc

dfloat fvcVolumetricSourceScalar
(
	dfloat GV,
	dfloat *phi,
	mesh_t *mesh
);

dfloat fvcVolumetricSourceVector
(
	dfloat GV,
	dfloat *phi,
	mesh_t *mesh
);

dfloat fvcDivergenceOfScalar
(
	dfloat GD,
	dfloat* Phif,
	dfloat* Flux,
	mesh_t *mesh
);

dfloat fvcDivergenceOfFlux
(
	dfloat GD,
	dfloat* Flux,
	mesh_t *mesh
);


//fvd
void fvdGradientOfScalar
(
	dfloat GG,
	dfloat *phif,
	mesh_t *mesh,
	dfloat* GradientOfScalar
);

//fvm
dfloat fvmDdtOfScalar
(
	dfloat GT,
	dfloat* Phi,
	dint PhiName,
	dint timeScheme,
	mesh_t *mesh
);
dfloat fvmDivergenceOfScalar
(
	dfloat GD,
	dfloat* Phi,
	dint PhiName,
	dfloat* Flux,
	dint interpolationScheme,
	dfloat blendingCoeff,
	mesh_t *mesh
);
dfloat fvmddtVectorWithScalarCoefficient
(
	dfloat GT,
	dfloat* alpha,
	dfloat* Phi,
	dint PhiName,
	dint timeScheme,
	mesh_t *mesh
);
dfloat fvmLaplacianOfVectorWithScalarCoefficient
(
	dfloat GL,
	dfloat* alpha_face,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
);
dfloat fvmDivergenceOfVector
(
	dfloat GD,
	dfloat* Phi,
	dint PhiName,
	dfloat* Flux,
	dint interpolationScheme,
	dfloat blendingCoeff,
	mesh_t *mesh
);
dfloat fvmLaplacianOfScalarWithVectorCoefficient
(
	dfloat GL,
	dfloat* alpha_f,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
);

dfloat fvmLaplacianOfScalar1DZ
(
	dfloat GL,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
);

dfloat fvmLaplacianOfVectorWithScalarCoefficient
(
	dfloat GL,
	dfloat* alpha_face,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
);


dfloat fvmLaplacianOfScalarWithScalarCoefficient
(
	dfloat GL,
	dfloat* alpha_face,
	dfloat* phi,
	dint phiName,
	mesh_t *mesh
);

//solvers
void mainNavierStokesSolve(mesh_t *mesh);

//interpolation
void interpolateCellToFaceLinear(dfloat* phiCell,dfloat* phiFace,dint fieldType,mesh_t*	mesh);
void computeCellFaceValues(
	dfloat* phiCell,
	dfloat* phiFace,
	dint 	fieldType,
	dfloat* Flux,
	dint 	interpolationScheme,
	dfloat 	blendingCoeff,
	mesh_t*	mesh
);


//limiter

dfloat computeLimiter
(
	dfloat Flux,
	dfloat phiP,
	dfloat phiN,
	dfloat* gradCP,
	dfloat* gradCN,
	dfloat* d,
	dint interpolationScheme,
	dfloat blendingCoeff,
	mesh_t*	mesh
);

dfloat computePhiC
(
	dfloat Flux,
	dfloat phiP,
	dfloat phiN,
	dfloat* gradCP,
	dfloat* gradCN,
	dfloat* d,
	mesh_t*	mesh
);

dfloat compute_r // Compute r ratio
(
	dfloat Flux,
	dfloat phiP,
	dfloat phiN,
	dfloat* gradCP,
	dfloat* gradCN,
	dfloat* d,
	mesh_t*	mesh
);

#endif
