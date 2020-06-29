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

\*---------------------------------------------------------------------------*/

#include "StoR.H"
//#include "zeroGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::StoR<ReactionThermo>::StoR
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    laminar<ReactionThermo>(modelType, thermo, turb, combustionProperties),
	PDF_Name_(this->coeffs().template lookupOrDefault<word>("temperaturePDF", "NormalDistribution")),
	spanZoneForPDF_(3),
	truncationForPDF_(this->coeffs().lookupOrDefault("truncationForPDF",  1)),
	P_i_( spanZoneForPDF_ ),//The posibility for Zone i in spanZoneForPDF;
	T_i_( spanZoneForPDF_ ),//The normalized temperature value for Zone i;
	spaceFilter_( this->mesh(), 12 ),//coeff_=12, means that use near 6 points to calculate the filterd data
	deviationSimilarCoeff_(this->coeffs().lookupOrDefault("deviationSimilarCoeff", 1)),
	Tsgs_
	(
		IOobject
		(
			"Tsgs",
			this->mesh().time().timeName(),
			this->mesh(),
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		this->mesh(),
		dimensionedScalar("Tsgs", dimTemperature, 0.0 ),
		fixedValueFvPatchField<scalar>::typeName
	),
	Qdot_
    (
		IOobject
		(
			"Qdot",
			this->mesh().time().timeName(),
			this->mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE,
			false
		),
		this->mesh(),
		dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
    ),
	R_(this->thermo().composition().Y().size())
{
	forAll(R_, i)
    {
        R_.set
        (
            i,
            new fvScalarMatrix
            (
				this->thermo().composition().Y()[i],
				dimMass/dimTime
            )
        );
    }

	forAll(R_, speciesI) R_[speciesI] = 0.0*R_[speciesI];   


    truncationForPDF_ = 1; //Default Number

    P_i_[0] = 0.24; P_i_[1] = 0.52; P_i_[2] = 0.24;//Default Number
    T_i_[0] = -1;   T_i_[1] = 0;    T_i_[2] = 1;   //Default Number
//	P_i_ = this->Alpha_Calculation( spanZoneForPDF_, (-1)*truncationForPDF_, truncationForPDF_, PDF_Name_);
//	T_i_ = this->NormT_Calculation( spanZoneForPDF_, (-1)*truncationForPDF_, truncationForPDF_, PDF_Name_);	

    Info<< "StoR Combustion Model constructed:" << endl
		<< " presumed PDF name is " << PDF_Name_ << endl
        << " and Truncations are " << truncationForPDF_ << endl
		<< " normalized variances are  " << T_i_ << endl
		<< " posibility coeff array are " << P_i_ << endl;	

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::StoR<ReactionThermo>::~StoR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::StoR<ReactionThermo>::correct()
{

	forAll(R_, speciesI) R_[speciesI] = 0.0*R_[speciesI];
	Qdot_ *= 0.0;    
	Tsgs_ *= 0.0;

    if (this->active())
    {
		volScalarField& T = this->thermo().T();
		const volScalarField T_org = this->thermo().T();
	
		volScalarField T_Filter = 0*T_org;
		T_Filter = spaceFilter_(T_org);    
	//	volScalarField T_rms2T  = 0*Tsgs_*Tsgs_;
	//	T_rms2T  = (T_Filter - MZ_1)*(T_Filter - MZ_1);
		Tsgs_ = sqrt( deviationSimilarCoeff_ * mag(spaceFilter_((T_Filter - T_org)*(T_Filter - T_org))));
		Tsgs_.correctBoundaryConditions();
	

		
		scalar deltaT = 0.;
		deltaT = this->mesh().time().deltaTValue();	   
	
		for( label i=0; i< spanZoneForPDF_; i++)
		{
//			Info<< "4444spanZoneForPDF_ " << i << endl;	
			T = T_i_[i]*Tsgs_ + T_org;									//Update the stochastic temperature
//			Info<< "5555spanZoneForPDF_ " << i << endl;	
			this->chemistryPtr_->solve(deltaT);							//Run the reaction calculation, get reaction rate
//			Info<< "6666spanZoneForPDF_ " << i << endl;	
			Qdot_.ref() += P_i_[i] * this->chemistryPtr_->Qdot();		//Get the weighted reaction rate for a grid
			forAll(R_, speciesI)										//Get the weighted reaction source term for each species
				R_[speciesI] += P_i_[i] * this->chemistryPtr_->RR(speciesI);				
		}
		
		T = T_org;
    }

}

template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::StoR<ReactionThermo>::R(volScalarField& Y) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    
    fvScalarMatrix& Su = tSu.ref();

    if (this->active())
    {
        const label specieI =
            this->thermo().composition().species()[Y.member()];
        Su += R_[specieI];
    }

    return tSu;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::StoR<ReactionThermo>::Qdot() const
{
    return Qdot_;
}


template<class ReactionThermo>
bool Foam::combustionModels::StoR<ReactionThermo>::read()
{
    if (laminar<ReactionThermo>::read())
    {
		PDF_Name_ = this->coeffs().template lookupOrDefault<word>("temperaturePDF", "NormalDistribution");
		truncationForPDF_ = this->coeffs().lookupOrDefault("truncationForPDF",  1);
		deviationSimilarCoeff_ = this->coeffs().lookupOrDefault("deviationSimilarCoeff", 1);
		spanZoneForPDF_ = 3;

//		P_i_ = Alpha_Calculation( spanZoneForPDF_, (-1)*truncationForPDF_, truncationForPDF_, PDF_Name_);
//		T_i_ = NormT_Calculation( spanZoneForPDF_, (-1)*truncationForPDF_, truncationForPDF_, PDF_Name_);

		Info<< "StoR Combustion Model constructed:" << endl
			<< " presumed PDF name is " << PDF_Name_ << endl
			<< " and Truncations are " << truncationForPDF_ << endl
			<< " normalized variances are  " << T_i_ << endl
			<< " posibility coeff array are " << P_i_ << endl;	

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //



template<class Type>
Foam::tmp<Foam::scalar>
Foam::combustionModels::StoR<Type>::PDF_Norm
(
scalar phi,
scalar phi_average,
scalar sigma
)
{
	return	((phi-phi_average)/sigma);
}

template<class Type>
Foam::tmp<Foam::scalar>
Foam::combustionModels::StoR<Type>::PDF_Density(scalar normPhi, word PDFname)
{
	if(PDFname == "NormalDistribution")
		return	(1.0/Foam::sqrt(2.0*Foam::constant::mathematical::pi))*Foam::exp(-1*(Foam::sqr(normPhi)/2.0));
	else
		return 0;
}

template<class Type>
Foam::tmp<Foam::scalar>
Foam::combustionModels::StoR<Type>::PDF_Accumulate(scalar normPhi, word PDFname)
{
	if(PDFname == "NormalDistribution")
		return	(0.5*(1+Foam::erf(normPhi/Foam::sqrt(2.0))));	//accumulate PDF for normal distribution
	else
		return 0;
}

template<class Type>
Foam::tmp<Foam::scalar>
Foam::combustionModels::StoR<Type>::PDF_Expection(scalar normPhi, word PDFname)
{
	if(PDFname == "NormalDistribution")
		return ((-1)*this->PDF_Density(normPhi, PDFname));
	else
		return 0;
}

template<class Type>
Foam::tmp<Foam::scalar>
Foam::combustionModels::StoR<Type>::PDF_Expection(scalar normPhi_left, scalar normPhi_right, word PDFname)
{
	return ( this->PDF_Expection(normPhi_right, PDFname) - this->PDF_Expection(normPhi_left, PDFname) );
}

template<class Type>
Foam::tmp<Foam::scalar>
Foam::combustionModels::StoR<Type>::PDF_Density(scalar normPhi, scalar leftTruncation, scalar rightTruncation, word PDFname)
{
	if( (this->PDF_Accumulate(rightTruncation,PDFname)-this->PDF_Accumulate(leftTruncation,PDFname)) != 0 )
		return	(this->PDF_Density(normPhi, PDFname))/(this->PDF_Accumulate(rightTruncation,PDFname)-this->PDF_Accumulate(leftTruncation,PDFname));
	else
		return	0;
}

template<class Type>
Foam::tmp<Foam::scalar>
Foam::combustionModels::StoR<Type>::PDF_Accumulate(scalar normPhi, scalar leftTruncation, scalar rightTruncation, word PDFname)
{	
	if( (this->PDF_Accumulate(rightTruncation,PDFname)-this->PDF_Accumulate(leftTruncation,PDFname)) != 0 )
		return	(this->PDF_Accumulate(normPhi,PDFname)-this->PDF_Accumulate(leftTruncation,PDFname))/(this->PDF_Accumulate(rightTruncation,PDFname)-this->PDF_Accumulate(leftTruncation,PDFname));
	else
		return	0;
}

template<class Type>
Foam::tmp<Foam::scalar>
Foam::combustionModels::StoR<Type>::PDF_Expection(scalar normPhi, scalar leftTruncation, scalar rightTruncation, word PDFname)
{
	if( (PDF_Accumulate(rightTruncation,PDFname)-PDF_Accumulate(leftTruncation,PDFname)) != 0 )
		return (PDF_Expection(normPhi, PDFname))/(PDF_Accumulate(rightTruncation,PDFname)-PDF_Accumulate(leftTruncation,PDFname));
	else
		return	0;
}

template<class Type>
Foam::tmp<Foam::scalar>
Foam::combustionModels::StoR<Type>::PDF_Expection(scalar normPhi_left, scalar normPhi_right, scalar leftTruncation, scalar rightTruncation, word PDFname)
{
	if( (this->PDF_Accumulate(rightTruncation,PDFname)-this->PDF_Accumulate(leftTruncation,PDFname)) != 0 )
		return (this->PDF_Expection(normPhi_left, normPhi_right, PDFname))/(this->PDF_Accumulate(rightTruncation,PDFname)-this->PDF_Accumulate(leftTruncation,PDFname));
	else
		return	0;
}
/*
template<class Type>
Foam::tmp<Foam::scalarField>
Foam::combustionModels::StoR<Type>::Alpha_Calculation
(
scalar SpanZone_Phi,
scalar Left_Truncation,
scalar Right_Truncation,
word PDF_Name
)
{

	scalarField	Norm_Phi(SpanZone_Phi, 0.0);
	scalarField	Alpha_Phi(SpanZone_Phi, 0.0);
	Norm_Phi[0] = this->PDF_Norm( 0, 0, -1);
	
	if(SpanZone_Phi == 1)
	{
		Norm_Phi[0]  = 0;
		Alpha_Phi[0] = 1;
	}
	else
	{
		for(label i = 0; i < SpanZone_Phi; i++)
		{
			label ispan = i;
			Norm_Phi[i] = Left_Truncation + ispan*(Right_Truncation-Left_Truncation)/(SpanZone_Phi-1);
		}
//Info << " 1111111111111111111111111111111111	" <<"	SpanZone_Phi	"<<SpanZone_Phi<< "	Norm_Phi	"<<Norm_Phi<< endl;
		for(label i = 0; i < SpanZone_Phi; i++)
		{

			if( i==0 )
			{
//Info << " 2222222222222222222222222222222222	" <<"Norm_Phi[i]	"<<Norm_Phi[i]<<"	Norm_Phi[i+1]"<<Norm_Phi[i+1]<<"	Sigma_Phi	"<<Sigma_Phi<< endl;
				Alpha_Phi[i] = ( PDF_Expection(Norm_Phi[i], Norm_Phi[i+1], Left_Truncation, Right_Truncation, PDF_Name) \
						- Norm_Phi[i+1]*(PDF_Accumulate(Norm_Phi[i+1], Left_Truncation, Right_Truncation, PDF_Name)-PDF_Accumulate(Norm_Phi[i], Left_Truncation, Right_Truncation, PDF_Name))) \
							/(Norm_Phi[i]-Norm_Phi[i+1]);
			}		
			if( i==(SpanZone_Phi-1) )
			{
//Info << " 3333333333333333333333333333333333	" <<"Norm_Phi[i-1]	"<<Norm_Phi[i-1]<<"	Norm_Phi[i]	"<<Norm_Phi[i]<<"	Sigma_Phi	"<<Sigma_Phi<< endl;
				Alpha_Phi[i] = (-1)*((PDF_Expection(Norm_Phi[i-1], Norm_Phi[i], Left_Truncation, Right_Truncation, PDF_Name) \
						- Norm_Phi[i-1]*(PDF_Accumulate(Norm_Phi[i], Left_Truncation, Right_Truncation, PDF_Name)-PDF_Accumulate(Norm_Phi[i-1], Left_Truncation, Right_Truncation, PDF_Name))) \
							/(Norm_Phi[i-1]-Norm_Phi[i]));
			}
			if( (i!=0) && (i!=(SpanZone_Phi-1)) )
			{	
//Info << " 4444444444444444444444444444444444    " <<"Norm_Phi[i-1]"<<Norm_Phi[i-1]<<"	Norm_Phi[i]	"<<Norm_Phi[i]<<"	Norm_Phi[i+1]	"<<Norm_Phi[i+1]<<"	Norm_Phi[i+1]	"<<Sigma_Phi<< endl;
				Alpha_Phi[i] = ( PDF_Expection(Norm_Phi[i], Norm_Phi[i+1], Left_Truncation, Right_Truncation, PDF_Name) \
						- Norm_Phi[i+1]*(PDF_Accumulate(Norm_Phi[i+1], Left_Truncation, Right_Truncation, PDF_Name)-PDF_Accumulate(Norm_Phi[i], Left_Truncation, Right_Truncation, PDF_Name))) \
							/(Norm_Phi[i]-Norm_Phi[i+1])
						+ \
								(-1)*((PDF_Expection(Norm_Phi[i-1], Norm_Phi[i], Left_Truncation, Right_Truncation, PDF_Name) \
						- Norm_Phi[i-1]*(PDF_Accumulate(Norm_Phi[i], Left_Truncation, Right_Truncation, PDF_Name)-PDF_Accumulate(Norm_Phi[i-1], Left_Truncation, Right_Truncation, PDF_Name))) \
							/(Norm_Phi[i-1]-Norm_Phi[i]));
			}
		
		}				
	}

	return Alpha_Phi;

}
*/

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::combustionModels::StoR<Type>::NormT_Calculation
(
scalar SpanZone_Phi,
scalar Left_Truncation,
scalar Right_Truncation,
word PDF_Name
)
{

	scalarField	Norm_Phi(SpanZone_Phi, 0.0);
//	scalarField	Alpha_Phi(SpanZone_Phi, 0.0);

	if(SpanZone_Phi == 1)
	{
		Norm_Phi[0]  = 0;
//		Alpha_Phi[0] = 1;
	}
	else
	{
		for(label i = 0; i < SpanZone_Phi; i++)
		{
			label ispan = i;
			Norm_Phi[i] = Left_Truncation + ispan*(Right_Truncation-Left_Truncation)/(SpanZone_Phi-1);
		}
		
	}

	return Norm_Phi;
}



