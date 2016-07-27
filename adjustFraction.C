/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Author                                                               
    Mauro Bracconi <mauro.braccon@polimi.it>                             
    Department of Energy                                                  
    Politecnico di Milano                                                
    via La Masa, 34 - 20156 - Milano, Italy   
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    foamAdjustFractions

Description
    Utility to move a patch along the normal direction in order to get
    the user provided solid or void fraction
    

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dynamicFvMesh.H"
#include "pimpleControl.H"
#include "pointPatchField.H"
#include "surfaceFields.H"
#include "pointMesh.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createDynamicFvMesh.H"
	#include "readOptions.H"
	
	// Ident
   	label patchID = mesh.boundaryMesh().findPatchID(surfacePatch_);
   	 
   	pointVectorField& displacement = const_cast<pointVectorField&>
	(
		mesh.objectRegistry::lookupObject<pointVectorField>
		(
				"pointDisplacement"
		)
	);
			
	while (runTime.loop())
   	{
		Info<< "Time = " << runTime.timeName() << endl;
		
		vectorField& pointPosition = refCast<vectorField>(displacement.boundaryField()[patchID]);
		vectorField pointNormal   = mesh.boundaryMesh()[patchID].pointNormals();

        //Find the relevant size of the vector and declare a vectorField.
        label size = pointPosition.size();
        vectorField dispPoints(size);

		// Evaluate total volume of the bouding box
		const boundBox& box = mesh.bounds();
		Vector<double> min = box.min();
		Vector<double> max = box.max();
		scalar bbVolume = (max[0]-min[0])*(max[1]-min[1])*(max[2]-min[2]);

		// Evaluate mesh volume
		const scalarField& vol = mesh.V();
        scalar meshVolume = gSum(vol); 
		
		// Evaluate void fraction
		scalar voidFrac = 0.0;
		if( fractionType_ )
		{
			voidFrac = meshVolume/bbVolume;
		} 
		else
		{
			voidFrac = 1. - meshVolume/bbVolume;
		}
		
		// Evaluate error on the void fraction fraction
		scalar epsiError = (voidFrac - epsiTarget_) / epsiTarget_;
		
		// Define the max displacement of a node, equal to 1/4 of minimum cell dimension		
		scalar maxBuf = 1./(4.*gMin(mesh.surfaceInterpolation::deltaCoeffs()));
		
		// Define displacement of nodes
		scalar buf = Foam::min(maxBuf,(fractionType_ ? -1.0 : 1.0) * fractionMove_ * epsiError);

		// Loop over nodes
		forAll(dispPoints, index)
		{
			dispPoints[index].x() = pointPosition[index].x() + buf * pointNormal[index].x();
			dispPoints[index].y() = pointPosition[index].y() + buf * pointNormal[index].y();
			dispPoints[index].z() = pointPosition[index].z() + buf * pointNormal[index].z();
		}

		//Once the values have been assigned to dispPoints, assign them to PointDisplacement boundaryField
		displacement.boundaryField()[patchID] == dispPoints;
					
		// Update mesh points
        mesh.update();
 
		if( checkMesh_ )
		{ 
			mesh.checkMesh(true);	
		}

		const scalarField& vol1 = mesh.V();
		meshVolume = gSum(vol1);
		
		voidFrac = 0.0;
		if( fractionType_ )
		{
			voidFrac = meshVolume/bbVolume;
		} 
		else
		{
			voidFrac = 1. - meshVolume/bbVolume;
		}
		
		epsiError = (voidFrac - epsiTarget_) / epsiTarget_;

		const scalarField& Ap = mesh.magSf().boundaryField()[patchID];
        scalar patchArea = gSum(Ap);
				
		Info << endl;
		Info << "Foam geometrical properties: " << endl;
		Info << " * Sv            	= " <<  patchArea/bbVolume << "	[L2/L3]" << endl;
		Info << " * Void fraction 	= " <<  voidFrac << "	[-]" << endl;
		Info << "\nError on the void fraction: " <<  epsiError << nl << endl;

		runTime.write();

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;     
    }

	mesh.write();
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
