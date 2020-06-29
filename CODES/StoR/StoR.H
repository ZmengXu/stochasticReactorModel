/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

Class
    Foam::combustionModels::StoR

Description
    Eddy Dissipation Concept (StoR) turbulent combustion model.

    This model considers that the reaction occurs in the regions of the flow
    where the dissipation of turbulence kinetic energy takes place (fine
    structures). The mass fraction of the fine structures and the mean residence
    time are provided by an energy cascade model.

    There are many versions and developments of the StoR model, 4 of which are
    currently supported in this implementation: v1981, v1996, v2005 and
    v2016.  The model variant is selected using the optional \c version entry in
    the \c StoRCoeffs dictionary, \eg

    \verbatim
        StoRCoeffs
        {
            version v2016;
        }
    \endverbatim

    The default version is \c v2015 if the \c version entry is not specified.

    Model versions and references:
    \verbatim
        Version v2005:

            Cgamma = 2.1377
            Ctau = 0.4083
            kappa = gammaL^exp1 / (1 - gammaL^exp2),

            where exp1 = 2, and exp2 = 2.

            Magnussen, B. F. (2005, June).
            The Eddy Dissipation Concept -
            A Bridge Between Science and Technology.
            In ECCOMAS thematic conference on computational combustion
            (pp. 21-24).

        Version v1981:

            Changes coefficients exp1 = 3 and exp2 = 3

            Magnussen, B. (1981, January).
            On the structure of turbulence and a generalized
            eddy dissipation concept for chemical reaction in turbulent flow.
            In 19th Aerospace Sciences Meeting (p. 42).

        Version v1996:

            Changes coefficients exp1 = 2 and exp2 = 3

            Gran, I. R., & Magnussen, B. F. (1996).
            A numerical study of a bluff-body stabilized diffusion flame.
            Part 2. Influence of combustion modeling and finite-rate chemistry.
            Combustion Science and Technology, 119(1-6), 191-217.

        Version v2016:

            Use local constants computed from the turbulent Da and Re numbers.

            Parente, A., Malik, M. R., Contino, F., Cuoci, A., & Dally, B. B.
            (2016).
            Extension of the Eddy Dissipation Concept for
            turbulence/chemistry interactions to MILD combustion.
            Fuel, 163, 98-111.
    \endverbatim

SourceFiles
    StoR.C

\*---------------------------------------------------------------------------*/

#ifndef StoR_H
#define StoR_H

#include "../laminar/laminar.H"
#include "Enum.H"
#include "laplaceFilter.H"	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
/*
//- StoR model versions
enum class StoRversions
{
    v1981,
    v1996,
    v2005,
    v2016
};

extern const Enum<StoRversions> StoRversionNames;
extern const StoRversions StoRdefaultVersion;

const scalar StoRexp1[] = {3, 2, 2, 2};
const scalar StoRexp2[] = {3, 3, 2, 2};
*/
/*---------------------------------------------------------------------------*\
                            Class StoR Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class StoR
:
    public laminar<ReactionThermo>
{
    // Private data

        //- The selected model version
//        StoRversions version_;

	//- The followings are for StoR combustion model
	//- The parameter for presumed PDF 
	word PDF_Name_;
	scalar spanZoneForPDF_;
	scalar truncationForPDF_;
	scalarField P_i_;//The posibility for Zone i in spanZoneForPDF;
	scalarField T_i_;//The normalized temperature value for Zone i;

    //- Filter parameter
	laplaceFilter spaceFilter_;
	scalar deviationSimilarCoeff_;//similarity model for deviation;

	//- volScalarField
	volScalarField Tsgs_;
	volScalarField Qdot_;
	PtrList<fvScalarMatrix> R_;


    // Private Member Functions

	//- The followings are for StoR combustion model
        //- The functions for presumed PDF
	tmp<scalar> PDF_Norm(scalar phi, scalar phi_average, scalar sigma);
	tmp<scalar> PDF_Density(scalar normPhi, word PDFname);
	tmp<scalar> PDF_Accumulate(scalar normPhi, word PDFname);
	tmp<scalar> PDF_Expection(scalar normPhi, word PDFname);
	tmp<scalar> PDF_Expection(scalar normPhi_left, scalar normPhi_right, word PDFname);
	tmp<scalar> PDF_Density(scalar normPhi, scalar leftTruncation, scalar rightTruncation, word PDFname);
	tmp<scalar> PDF_Accumulate(scalar normPhi, scalar leftTruncation, scalar rightTruncation, word PDFname);
	tmp<scalar> PDF_Expection(scalar normPhi, scalar leftTruncation, scalar rightTruncation, word PDFname);
	tmp<scalar> PDF_Expection(scalar normPhi_left, scalar normPhi_right, scalar leftTruncation, scalar rightTruncation, word PDFname);
	tmp<scalarField> Alpha_Calculation( scalar SpanZone_Phi, scalar Left_Truncation, scalar Right_Truncation, word PDF_Name);
	tmp<scalarField> NormT_Calculation( scalar SpanZone_Phi, scalar Left_Truncation, scalar Right_Truncation, word PDF_Name);



        //- No copy construct
        StoR(const StoR&) = delete;

        //- No copy assignment
        void operator=(const StoR&) = delete;


public:

    //- Runtime type information
    TypeName("StoR");


    // Constructors

        //- Construct from components
        StoR
        (
            const word& modelType,
            ReactionThermo& type,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~StoR();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix.
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "StoR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
