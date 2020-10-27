
label N = Y.size();
label M = reactions.size();
label NM = N*M; 

PtrList<DimensionedField<scalar, volMesh>> RR_(NM);//reaction rate for species n in reaction m
PtrList<DimensionedField<scalar, volMesh>> SH_(M);//enthorpy source in each reaction
PtrList<DimensionedField<scalar, volMesh>> rr_(N);//species consumption for each species, calculated from chemistry model rr_[i] = chemistry.RR(i);


PtrList<dimensionedScalar> hcSp(N);//The formation enthapy of each species
forAll(Y, i)
{
	const scalar hi = specieData[i].Hc();
    hcSp.set
    (
		i,
        new dimensionedScalar("hi",dimEnergy/dimMass,hi)
    );
}


forAll(Y,i)
{
	word name("RR_"  + Y[i].name());
	rr_.set
	(
		i,
		new DimensionedField<scalar , volMesh>
		(
			IOobject
			(
				name,
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
		)
	);
}



forAll(reactions,m)
{
	word name("SH_"  + reactions[m].name());
	SH_.set
	(
		m,
		new DimensionedField<scalar , volMesh>
		(
			IOobject
			(
				name,
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
		)
	);
}




label I = 0;
forAll(reactions,m)
{
	forAll(Y,n)
	{ 
		forAll(reactions[m].lhs(),s)
		{
			if (n == reactions[m].lhs()[s].index  )
			{
				I++;
				word name( Y[n].name() + "_" + reactions[m].name());
				RR_.set
				(
					I,
					new DimensionedField<scalar, volMesh>
					(
						IOobject
						(
							name,
							mesh.time().timeName(),
							mesh,
							IOobject::NO_READ,
							IOobject::AUTO_WRITE
						),
						mesh,
						dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
					)
				);
			}
		}

		forAll(reactions[m].rhs(),s)
		{
			if (n == reactions[m].rhs()[s].index )
			{
				I++;
				word name( Y[n].name() + "_" + reactions[m].name());
				RR_.set
				(
					I,
					new DimensionedField<scalar, volMesh>
					(
						IOobject
						(
						name,
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
						),
						mesh,
						dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
					)
				);
			}
		}
	}
}
