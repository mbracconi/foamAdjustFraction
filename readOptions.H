Info<< "\nReading foamAdjustDict options\n" << endl;
IOdictionary foamOptions
(
	IOobject
	(
		"foamAdjustDict",
		runTime.constant(),
		runTime.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

// Surface patch
word  surfacePatch_(foamOptions.lookup("foamSurfacePatch"));

// Switch -> off solid fraction | on void fraction
Switch  fractionType_(foamOptions.lookup("voidFraction"));

// Switch for checkMesh
Switch  checkMesh_(foamOptions.lookupOrDefault("checkMesh",false));

// Target void fraction
scalar epsiTarget_ = readScalar(foamOptions.lookup("epsiTarget"));

scalar fractionMove_ = readScalar(foamOptions.lookup("fractionMove"));

