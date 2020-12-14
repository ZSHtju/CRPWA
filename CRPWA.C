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
    Instaneous chemical reaction path way analyze

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
		"reactionRateSpecies",
		"list",
		"specify a list of species to be processed with reaction rate, e.g. '( O2 CH4 H2)' -"
		"otherwise all species will be processed"
	);
	argList::addOption
	(
		"pathFluxSpecies",
		"list",
		"specify a list of species to be processed with path flux analysis, e.g. '( O2 CH4 H2)' -"
		"otherwise all species will be processed"
	);
	argList::addBoolOption
	(
		"noElementaryReaction",
		"No output for elementary reaction information, only species reaction rate"
	);
	argList::addOption
	(
		"positions",
		"vectorList" ,
		"path flux position by specified <vectorList> - eg, '((1 0 0)(2 0 0))' "
	);
	
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
	List<vector> positions;
	HashSet<word> reactionRateSpecies;
	HashSet<word> pathFluxSpecies;
	if(args.optionFound("positions"))
	{
		args.optionLookup("positions" )() >> positions;
	}
	if(args.optionFound("reactionRateSpecies"))
	{
		args.optionLookup("reactionRateSpecies" )() >> reactionRateSpecies;
	}
	if(args.optionFound("pathFluxSpecies"))
	{
		args.optionLookup("pathFluxSpecies" )() >> pathFluxSpecies;
	}
	const bool noElementaryReaction = args.optionFound("noElementaryReaction");
    instantList timeDirs = timeSelector::select0(runTime, args);
	#include "postProcess.H"
	// output the average heat release rate of each elementary reaction
	OFstream  HRROutput( "log.Reactions");
	OFstream  HRROutputNormalized( "log.ReactionsNormalized");
	//
	label timeIndex = 0;
	forAll(timeDirs, timeI)
	{
		#include "createMesh.H"
		#include "createFields.H"
		#include "createNewFields.H"
		timeIndex++;
		
		labelList cellID;
		if(args.optionFound("positions"))
		{
			cellID.setSize(positions.size(),0);
			forAll(positions,i)
			{
				cellID[i] = mesh.findCell(positions[i]);
			}
		}
	
		if(timeIndex == 1)
		{
			HRROutput << "Time [s or CA]" << tab << "HRTotal [J/s]" << tab; 
			forAll(reactions,m)
			{
				HRROutput << "R" << m << "[J/s]" << tab;
			}
			HRROutput << endl;
			
			HRROutputNormalized << "Time [s or CA]" << tab << "HRTotal" << tab; 
			forAll(reactions,m)
			{
				HRROutputNormalized << "R" << m << tab;
			}
			HRROutputNormalized << endl;
		}
		
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
			rr_[i].ref() = chemistry.RR(i);
			if(reactionRateSpecies.found(Y[i].name()) || !args.optionFound("reactionRateSpecies"))
			{
				rr_[i].write();
			}
		}
       

//**************************************************************//


		//label II = 0;
		forAll(reactions,m)
		{   
			const Reaction<gasHThermoPhysics>& R = reactions[m];
			forAll(mesh.C(),celli)
			{
				scalar pf,cf,pr,cr;
				label lRef,rRef;
				scalar Ti = thermo.T()()[celli];
				scalar pi = thermo.p()()[celli];
				scalar rhoi= thermo.rho()()[celli];
				scalarField c(Y.size(),0);
				forAll(Y,i)
				{
					c[i] = rhoi*Y[i][celli]/specieData[i].W();
				}
				scalar omegai = R.omega(pi,Ti,c,pf,cf,lRef,pr,cr,rRef);// get the instantaneous reaction rates
	
				forAll(Y,n)
				{ 
					if(reactionRateSpecies.found(Y[n].name()) || !args.optionFound("reactionRateSpecies"))
					{
						forAll(R.lhs(),s)
						{
							if (n == R.lhs()[s].index  )
							{
								label II = m*(n-1)+n;
								RR_[II][celli] -=  omegai*R.lhs()[s].stoichCoeff*specieData[n].W();   
								SH_[m][celli] -= RR_[II][celli]*hcSp[n].value();
							}
						}

						forAll(R.rhs(),s)
						{
							if (n == R.rhs()[s].index )
							{
								label II = m*(n-1)+n;
								RR_[II][celli] +=  omegai*R.rhs()[s].stoichCoeff*specieData[n].W();   
								SH_[m][celli] -= RR_[II][celli]*hcSp[n].value();
							}
						}
					}
				}
			}
		}
		// species path flux analysis	
				//***********************************************// 
		// global
		forAll(Y,n)
		{   
			if(pathFluxSpecies.found(Y[n].name()) || !args.optionFound("pathFluxSpecies"))
			{
				autoPtr<OFstream> speciesRRGlobal;
				autoPtr<OFstream> speciesRRLocal;
				if (speciesRRGlobal.empty())
				{
					// File update
					//fileName output;
					//word name_ = "speciesNamePtr";

					// Open new file at start up
					speciesRRGlobal.reset(new OFstream("log.pathFlux_"+Y[n].name()));	 
				}
				if(timeIndex == 1)
				{
					speciesRRGlobal() << Y[n].name() << tab; 
					forAll(reactions,m)
					{ 
						const Reaction<gasHThermoPhysics>& R = reactions[m];
						forAll(R.lhs(),s)
						{
							if (n == R.lhs()[s].index  )
							{
								speciesRRGlobal() << "R" << m << "[kg/s]" << tab;
							}
						}

						forAll(R.rhs(),s)
						{
							if (n == R.rhs()[s].index )
							{
								speciesRRGlobal() << "R" << m << "[kg/s]" << tab;        
							}
						}
					}
					speciesRRGlobal() << endl;
				}
				scalar volIntegrateRRi = 0.0;
				forAll(mesh.C(),celli)
				{
					volIntegrateRRi += rr_[n][celli]*mesh.V()[celli];
				}
				speciesRRGlobal() << scientific << volIntegrateRRi << tab;
				forAll(reactions,m)
				{ 
					const Reaction<gasHThermoPhysics>& R = reactions[m];
					forAll(R.lhs(),s)
					{
						if (n == R.lhs()[s].index  )
						{
							label II = m*(n-1)+n;
							scalar volIntegrateRR = 0.0;
							forAll(mesh.C(),celli)
							{
								volIntegrateRR += RR_[II][celli]*mesh.V()[celli];
							}
							speciesRRGlobal() << scientific << volIntegrateRR << tab;
						}
					}

					forAll(R.rhs(),s)
					{
						if (n == R.rhs()[s].index )
						{
							label II = m*(n-1)+n;
							scalar volIntegrateRR = 0.0;
							forAll(mesh.C(),celli)
							{
								volIntegrateRR += RR_[II][celli]*mesh.V()[celli];
							}
							speciesRRGlobal() << scientific << volIntegrateRR << tab;   
						}
					}
				}
				speciesRRGlobal() << endl;
			}
		}
		//local 
		forAll(Y,n)
		{   
			if(pathFluxSpecies.found(Y[n].name()) || !args.optionFound("pathFluxSpecies"))
			{
				forAll(cellID,i)
				{
					autoPtr<OFstream> speciesRRLocal;
					if (speciesRRLocal.empty())
					{
						// File update
						//fileName output;
						//word name_ = "speciesNamePtr";

						// Open new file at start up
						speciesRRLocal.reset(new OFstream("log.pathFlux_"+Y[n].name()+"_"+Foam::name(positions[i])));	 
					}
					if(timeIndex == 1)
					{
						speciesRRLocal() << Y[n].name() << tab; 
						forAll(reactions,m)
						{ 
							const Reaction<gasHThermoPhysics>& R = reactions[m];
							forAll(R.lhs(),s)
							{
								if (n == R.lhs()[s].index  )
								{
									speciesRRLocal() << "R" << m << "[kg/s]" << tab;
								}
							}

							forAll(R.rhs(),s)
							{
								if (n == R.rhs()[s].index )
								{
									speciesRRLocal() << "R" << m << "[kg/s]" << tab;        
								}
							}
						}
						speciesRRLocal() << endl;
					}
					forAll(mesh.C(),celli)
					{
						if(celli == cellID[i])
						{
							speciesRRLocal() << scientific << rr_[n][celli]*mesh.V()[celli] << tab;
						}
					}

					forAll(reactions,m)
					{ 
						const Reaction<gasHThermoPhysics>& R = reactions[m];
						forAll(R.lhs(),s)
						{
							if (n == R.lhs()[s].index  )
							{
								label II = m*(n-1)+n;
								forAll(mesh.C(),celli)
								{
									if(celli == cellID[i])
									{
										speciesRRLocal() << scientific << RR_[II][celli]*mesh.V()[celli] << tab;
									}
								}
							}
						}

						forAll(R.rhs(),s)
						{
							if (n == R.rhs()[s].index )
							{
								label II = m*(n-1)+n;
								forAll(mesh.C(),celli)
								{
									if(celli == cellID[i])
									{
										speciesRRLocal() << scientific << RR_[II][celli]*mesh.V()[celli] << tab;
									}
								}
							}
						}
					}
					speciesRRLocal() << endl;
				}
			}
		}
		//global
		scalar volIntegrateChemistryhsSource = 0.0;
		forAll(mesh.C(), celli)
		{
			volIntegrateChemistryhsSource += chemistryhsSource[celli]*mesh.V()[celli];
		}
		HRROutput << runTime.timeName() << tab << scientific << volIntegrateChemistryhsSource << tab;
		HRROutputNormalized << runTime.timeName() << tab << scientific << volIntegrateChemistryhsSource/(volIntegrateChemistryhsSource+SMALL) << tab;
		forAll(reactions,m)
		{   
			const Reaction<gasHThermoPhysics>& R = reactions[m];	
			forAll(Y,n)
			{ 
				if(reactionRateSpecies.found(Y[n].name()) || !args.optionFound("reactionRateSpecies"))
				{
					forAll(R.lhs(),s)
					{
						if (n == R.lhs()[s].index  )
						{
							label II = m*(n-1)+n;
							RR_[II].write();
						}
					}

					forAll(R.rhs(),s)
					{
						if (n == R.rhs()[s].index )
						{
							label II = m*(n-1)+n;
							RR_[II].write();  
						}
					}
				}
			}
			if(!noElementaryReaction)
			{
				SH_[m].write();
			}
			scalar volIntegrateSH = 0.0;
			forAll(mesh.C(), celli)
			{
				volIntegrateSH += SH_[m][celli]*mesh.V()[celli];
			}			
			HRROutput << scientific << volIntegrateSH << tab; 
			HRROutputNormalized << scientific << volIntegrateSH/(volIntegrateChemistryhsSource+SMALL) << tab; 
		}
		HRROutput << endl;
		HRROutputNormalized << endl;

		//local
		forAll(Y,n)
		{   
			forAll(cellID,i)
			{
				autoPtr<OFstream> HRLocal;
				if (HRLocal.empty())
				{
					// File update
					//fileName output;
					//word name_ = "speciesNamePtr";

					// Open new file at start up
					HRLocal.reset(new OFstream("log.Reactions_"+Foam::name(positions[i])));	 
				}
				if(timeIndex == 1)
				{
					forAll(reactions,m)
					{ 
						const Reaction<gasHThermoPhysics>& R = reactions[m];
						forAll(R.lhs(),s)
						{
							if (n == R.lhs()[s].index  )
							{
								HRLocal() << "R" << m << "[kg/s]" << tab;
							}
						}

						forAll(R.rhs(),s)
						{
							if (n == R.rhs()[s].index )
							{
								HRLocal() << "R" << m << "[kg/s]" << tab;        
							}
						}
					}
					HRLocal() << endl;
				}
				forAll(mesh.C(),celli)
				{
					if(celli == cellID[i])
					{
						HRLocal() << scientific << chemistryhsSource[celli]*mesh.V()[celli] << tab;
					}
				}

				forAll(reactions,m)
				{ 
					forAll(mesh.C(),celli)
					{
						if(celli == cellID[i])
						{
							HRLocal() << scientific << SH_[m][celli]*mesh.V()[celli] << tab;
						}
					}
				}
				HRLocal() << endl;
			}
		}


		//***********************************************//
		forAll(Y,n)
		{   
			if(reactionRateSpecies.found(Y[n].name()) || !args.optionFound("reactionRateSpecies") || pathFluxSpecies.found(Y[n].name()) || !args.optionFound("pathFluxSpecies"))
			{
				autoPtr<OFstream> speciesInReaction;
				if (speciesInReaction.empty())
				{
					// File update
					//fileName output;
					//word name_ = "speciesNamePtr";

					// Open new file at start up
					speciesInReaction.reset(new OFstream(Y[n].name()));	 
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
						speciesInReaction() << m << tab << reactions[m] << endl;
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
