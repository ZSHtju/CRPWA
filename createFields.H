Info<< nl << "Reading thermophysicalProperties" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
autoPtr< BasicChemistryModel<psiReactionThermo>> pChemistry
(
  BasicChemistryModel<psiReactionThermo>::New(thermo)
);
BasicChemistryModel<psiReactionThermo>& chemistry = pChemistry();
thermo.validate(args.executable(), "h", "e");

//*********************delete*********************************//
//basicMultiComponentMixture& composition = thermo.composition();
//************************************************************//
basicSpecieMixture& composition = thermo.composition();
//****************************************************//
PtrList<volScalarField>& Y = composition.Y();
 
PtrList<gasHThermoPhysics> specieData(Y.size());
forAll(specieData, i)
{
    specieData.set
    (
        i,
        new gasHThermoPhysics
        (
            dynamic_cast<const reactingMixture<gasHThermoPhysics>&>
                (thermo).speciesData()[i]
        )
    );
}

const PtrList<Reaction<gasHThermoPhysics>> reactions = dynamic_cast<const reactingMixture<gasHThermoPhysics>&>(thermo);
