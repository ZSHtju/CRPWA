
label N = Y.size();
label M = reactions.size();
label NM = N*M; 


PtrList<volScalarField> RR_(NM);//reaction rate for species n in reaction m
PtrList<volScalarField> SH_(M);//enthorpy source in each reaction
PtrList<volScalarField> rr_(N);//species consumption for each species, calculated from chemistry model rr_[i] = chemistry.RR(i);

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
		new volScalarField
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
	word name("SH_r" +Foam::name(m));
	SH_.set
	(
		m,
		new volScalarField
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



forAll(reactions,m)
{
	forAll(Y,n)
	{ 
		forAll(reactions[m].lhs(),s)
		{
			if (n == reactions[m].lhs()[s].index  )
			{
				label II = N*m + n;
				word name( Y[n].name() + "_R" +Foam::name(m));
				RR_.set
				(
					II,
					new volScalarField
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
				label II = N*m + n;
				word name( Y[n].name() + "_R" +Foam::name(m));
				RR_.set
				(
					II,
					new volScalarField
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
