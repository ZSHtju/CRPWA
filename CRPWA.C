/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    chemical reaction path way analyze

Description
	Shenghui Zhong
	* 	July 03, 2019.
	
	This postProcessing tool had been successfully applied in following papers.
	Zhong, S; et al. Fuel 2018, 234, 1044-1054.
	Zhong, S; et al. Applied Energy 2020, 275, 115320.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiReactionThermo.H"
#include "BasicChemistryModel.H"
#include "reactingMixture.H"
#include "chemistrySolver.H"
#include "thermoPhysicsTypes.H"
#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	timeSelector::addOptions();
	argList::addOption
	(
		"species",
		"list",
		"specify a list of species to be processed, e.g. '( O2 CH4 H2)' -"
		"otherwise all species will be processed"
	);
	argList::addBoolOption
	(
		"noElementaryReaction",
		"No output for elementary reaction information, only species reaction rate"
	);
	
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

	HashSet<word> selectedSpecies;
	if(args.optionFound("species"))
	{
		args.optionLookup("species" )() >> selectedSpecies;
	}
	const bool noElementaryReaction = args.optionFound("noElementaryReaction");
    instantList timeDirs = timeSelector::select0(runTime, args);
	#include "postProcess.H"
	forAll(timeDirs, timeI)
	{
		#include "createMesh.H"
		#include "createFields.H"
		#include "createNewFields.H"

		volScalarField chemistryhsSource
		(
			IOobject
			(
				"chemistryhsSource",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("chemistryhsSource", dimEnergy/dimTime/dimVolume, 0.0)
		);
		
		
		thermo.correct();
		chemistry.calculate();
		chemistryhsSource = chemistry.Qdot();
		chemistryhsSource.write();
		forAll(Y,i)
		{
			rr_[i] = chemistry.RR(i);
			if(selectedSpecies.found(Y[i].name()) || !args.optionFound("species"))
			{
				rr_[i].write();
			}
		}
       

//**************************************************************//


		label II = 0;
		forAll(reactions,m)
		{   
			const Reaction<gasHThermoPhysics>& R = reactions[m];
			label lhsn = -1;
			forAll(Y,n)
			{ 
				if(selectedSpecies.found(Y[n].name()) || !args.optionFound("species"))
				{
					scalar countl = 0.;
					scalar countr = 0.;
					forAll(R.lhs(),s)
					{

					if (n == R.lhs()[s].index  )
					{
						II++;
						const label a = m;
						const label b = n;
						RR_[II] =  chemistry.calculateRR(a,b);   
						RR_[II].write();
						const scalar sl = R.lhs()[s].stoichCoeff;
						countl -= sl; 
						lhsn = n;     
						}
						}

						forAll(R.rhs(),s)
						{
							if (n == R.rhs()[s].index )
							{
								II++;
								const label a = m;
								const label b = n;
								RR_[II] =  chemistry.calculateRR(a,b);  
								RR_[II].write();
								const scalar sr = R.rhs()[s].stoichCoeff;
								countr += sr;             
							}
						}
						if (countl != 0 || countr != 0)
						{
							const label c = m;
							const label d = n;
							const DimensionedField<scalar, volMesh> TTT =  chemistry.calculateRR(c,d);
							SH_[m] -= TTT*hcSp[n]*(countl+countr);
						}
					}
			}
			//const label e = m;
			//const label f = lhsn;
			//const DimensionedField<scalar, volMesh> TTTT =  chemistry.calculateRR(e,f)/specieData[f].W();       
			if(!noElementaryReaction)
			{
				SH_[m].write();
			}
		}





		//***********************************************//
		forAll(Y,n)
		{   
			if(selectedSpecies.found(Y[n].name()) || !args.optionFound("species"))
			{
				autoPtr<OFstream> speciesNamePtr;
				if (speciesNamePtr.empty())
				{
					// File update
					fileName output;
					word name_ = "speciesNamePtr";

					// Open new file at start up
					speciesNamePtr.reset(new OFstream(Y[n].name()));	 
				}

				forAll(reactions,m)
				{ 
					const Reaction<gasHThermoPhysics>& R = reactions[m];
					scalar countl = 0.;
					scalar countr = 0.;
					forAll(R.lhs(),s)
					{
						if (n == R.lhs()[s].index  )
						{
							const scalar sl = R.lhs()[s].stoichCoeff;
							countl -= sl;       
						}
					}

					forAll(R.rhs(),s)
					{
						if (n == R.rhs()[s].index )
						{
							const scalar sr = R.rhs()[s].stoichCoeff;
							countr += sr;            
						}
					}

					if (countl != 0 || countr != 0)
					{
						//const label e = m;
						//const label f = n;
						speciesNamePtr() << m << tab << reactions[m] << endl;
						//const DimensionedField<scalar, volMesh> SS = chemistry.calculateRR(e,f);
						//DimensionedField<scalar, volMesh>  SSS = SS*(countl+countr);

	/* 					forAll(mesh.C(),celli)
						{
							if(celli == 0)
							{
								speciesNamePtr()<< reactions[m];
								speciesNamePtr()<< SSS[celli]  << endl;
							}
							else
							{ 
								speciesNamePtr()<<  SSS[celli] << endl; 
							}
						} */
					}
				}
			}
		}
}
//*****************************************************//

	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
		<< nl << endl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //;C()
