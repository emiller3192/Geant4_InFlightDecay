//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4NeutronDecay.cc                                                 //
//  Author: L.G. Sarmiento (Lund)                                             //
//  Date:   10 October 2015                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4NeutronDecay.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <iomanip>

//try to use the RandBreitWigner module
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandBreitWigner.h"

#include <fstream>

G4NeutronDecay::G4NeutronDecay(const G4ParticleDefinition* theParentNucleus,
                            const G4double& branch, const G4double& Qvalue,
                            const G4double& excitationE)
 : G4NuclearDecay("neutron decay", Neutron, excitationE), transitionQ(Qvalue)
{
  SetParent(theParentNucleus);  // Store name of parent nucleus, delete G4MT_parent
  SetBR(branch);

  SetNumberOfDaughters(2);
  G4IonTable* theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  G4int daughterZ = theParentNucleus->GetAtomicNumber();
  G4int daughterA = theParentNucleus->GetAtomicMass() - 1;
  SetDaughter(0, theIonTable->GetIon(daughterZ, daughterA, excitationE) );
  SetDaughter(1, "neutron");
}


G4NeutronDecay::~G4NeutronDecay()
{}


G4DecayProducts* G4NeutronDecay::DecayIt(G4double)
{



  // Get random number from Breit Wigner
  //TRandom3* randen       = new TRandom3(0);

  // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)  
  CheckAndFillParent();

  // Fill G4MT_daughters with neutron and residual nucleus (stored by SetDaughter)  
  CheckAndFillDaughters();

  G4double neutronMass = G4MT_daughters[1]->GetPDGMass();
  // Excitation energy included in PDG mass
  G4double nucleusMass = G4MT_daughters[0]->GetPDGMass();

  // Q value was calculated from atomic masses.
  // Use it to get correct neutron energy.

  double switchQ = 1000. * transitionQ;
  // double Q = transitionQ;

  // initialize a random engine
  // CLHEP::HepRandomEngine* anEngine;
  // CLHEP::RandBreitWigner BreitWigner = CLHEP::RandBreitWigner(anEngine);

  /* replaced with a function call
  * The BreitWignerQ function is called recursively to prevent over populating
  * the bins at upper and lower cuts

// below is a switch statement to determine the fuzziness of the energy
  switch((int)switchQ) {

    case 816:
    case 6718:
    case 10086:

        //KE = randen->BreitWigner(KE, 210);
        // Q = G4RandGauss::shoot(transitionQ,0.21);
      Q = BreitWigner.shoot(transitionQ,0.21);

      break;

    case 938:
    case 2046:

        //KE = randen->BreitWigner(KE, 200);
        // Q = G4RandGauss::shoot(transitionQ,0.2);
      Q = BreitWigner.shoot(transitionQ,0.2);
      break;

    case 1347:
    case 1532:
    case 4158:
    case 7516:

        //KE = randen->BreitWigner(KE, 230);
        // Q = G4RandGauss::shoot(transitionQ,0.23);
        Q = BreitWigner.shoot(transitionQ,0.23);

      break;

    case 3158:

        //KE = randen->BreitWigner(KE, 300);
        // Q = G4RandGauss::shoot(transitionQ,0.3);
        Q = BreitWigner.shoot(transitionQ,0.3);

      break;

    case 1368:

        //KE = randen->BreitWigner(KE, 045);
        // Q = G4RandGauss::shoot(transitionQ,0.045);
        Q = BreitWigner.shoot(transitionQ,0.045);

      break;

    case 88:
    case 3451:
    case 558:

        //KE = randen->BreitWigner(KE, 015);
        // Q = G4RandGauss::shoot(transitionQ,0.015);
        Q = BreitWigner.shoot(transitionQ,0.015);

      break;

    case 18:
    case 3385:

        //KE = randen->BreitWigner(KE, 010);
        // Q = G4RandGauss::shoot(transitionQ,0.010);
        Q = BreitWigner.shoot(transitionQ,0.010);

      break;

    case 2896:

        //KE = randen->BreitWigner(KE, 125);
        // Q = G4RandGauss::shoot(transitionQ,0.125);
        Q = BreitWigner.shoot(transitionQ,0.125);

      break;

    case 2186:

        //KE = randen->BreitWigner(KE, 200);
        // Q = G4RandGauss::shoot(transitionQ,0.2);
        Q = BreitWigner.shoot(transitionQ,0.2);

      break;

    default: 

        //KE = randen->BreitWigner(KE, 200);
        // Q = G4RandGauss::shoot(transitionQ,0.2);
        Q = BreitWigner.shoot(transitionQ,0.2);

     break;

  }

  if (Q <= 0) {

    Q = 0;
  } else if (Q > 20.6)
  {
    Q = 20.6;
  }
*/
  //double Q = G4RandGauss::shoot(transitionQ,0.2);
  double Q = G4NeutronDecay::BreitWignerQ(transitionQ);

  G4double cmMomentum = std::sqrt(Q*(Q + 2.*neutronMass)*
                                (Q + 2.*nucleusMass)*
                                (Q + 2.*neutronMass + 2.*nucleusMass) )/
                                (Q + neutronMass + nucleusMass)/2.; 

  // Set up final state
  // parentParticle is set at rest here because boost with correct momentum 
  // is done later
  G4DynamicParticle parentParticle(G4MT_parent, G4ThreeVector(0,0,0), 0.0);
  G4DecayProducts* products = new G4DecayProducts(parentParticle);

  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double sintheta = std::sqrt(1.0 - costheta*costheta);
  G4double phi  = twopi*G4UniformRand()*rad;
  G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),
                          costheta);

  G4double KE = std::sqrt(cmMomentum*cmMomentum + neutronMass*neutronMass)
              - neutronMass;

  // add 100 MeV of KE to test if the particle still comes to rest before decaying
  // G4double KE = std::sqrt(cmMomentum*cmMomentum + neutronMass*neutronMass)
              // - neutronMass + 100;


 // G4double KE = std::sqrt(cmMomentum*cmMomentum + neutronMass*neutronMass)
           //   - neutronMass;


  G4DynamicParticle* daughterparticle =
    new G4DynamicParticle(G4MT_daughters[1], direction, KE, neutronMass);
  products->PushProducts(daughterparticle);

  KE = std::sqrt(cmMomentum*cmMomentum + nucleusMass*nucleusMass) - nucleusMass;

  // add 100 MeV of KE to test if the particle still comes to rest before decaying
  // KE = std::sqrt(cmMomentum*cmMomentum + nucleusMass*nucleusMass) - nucleusMass + 100;

  
  daughterparticle =
    new G4DynamicParticle(G4MT_daughters[0], -1.0*direction, KE, nucleusMass);
  products->PushProducts(daughterparticle);

  // energies << direction.x() << "    " << direction.y() << "    " << direction.z() 
  // << "    " << KE << "    " << nucleusMass << std::endl;

  // Energy conservation check
  // For neutron decays, do final energy check against reaction Q value
  // which is well-measured using atomic mass differences.  Nuclear masses
  // should not be used since they are not usually directly measured and we
  // always decay atoms and not fully stripped nuclei.
  /*
  G4int nProd = products->entries();
  G4DynamicParticle* temp = 0;
  G4double Esum = 0.0;
  for (G4int i = 0; i < nProd; i++) {
    temp = products->operator[](i);
    Esum += temp->GetKineticEnergy();
  }
  G4double eCons = (transitionQ - Esum)/keV;
  if (eCons > 1.e-07) G4cout << " Neutron decay check: Ediff (keV) = " << eCons << G4endl;
  */

  //delete randen;
  return products;
}

G4double G4NeutronDecay::BreitWignerQ(G4double transitionQ)
{
  CLHEP::HepRandomEngine* anEngine;
  CLHEP::RandBreitWigner BreitWigner = CLHEP::RandBreitWigner(anEngine);
  G4double Q = transitionQ;
  double switchQ = 1000. * transitionQ;
  switch((int)switchQ) {

    case 816:
    case 6718:
    case 10086:

        //KE = randen->BreitWigner(KE, 210);
        // Q = G4RandGauss::shoot(transitionQ,0.21);
      Q = BreitWigner.shoot(transitionQ,0.21);

      break;

    case 938:
    case 2046:

        //KE = randen->BreitWigner(KE, 200);
        // Q = G4RandGauss::shoot(transitionQ,0.2);
      Q = BreitWigner.shoot(transitionQ,0.2);
      break;

    case 1347:
    case 1532:
    case 4158:
    case 7516:

        //KE = randen->BreitWigner(KE, 230);
        // Q = G4RandGauss::shoot(transitionQ,0.23);
        Q = BreitWigner.shoot(transitionQ,0.23);
        //Q = BreitWigner.shoot(transitionQ,0.09);

      break;

    case 3158:

        //KE = randen->BreitWigner(KE, 300);
        // Q = G4RandGauss::shoot(transitionQ,0.3);
        Q = BreitWigner.shoot(transitionQ,0.3);

      break;

    case 1368:

        //KE = randen->BreitWigner(KE, 045);
        // Q = G4RandGauss::shoot(transitionQ,0.045);
        Q = BreitWigner.shoot(transitionQ,0.045);

      break;

    case 88:
    case 3451:
    case 558:

        //KE = randen->BreitWigner(KE, 015);
        // Q = G4RandGauss::shoot(transitionQ,0.015);
        Q = BreitWigner.shoot(transitionQ,0.015);

      break;

    case 18:
    case 3385:

        //KE = randen->BreitWigner(KE, 010);
        // Q = G4RandGauss::shoot(transitionQ,0.010);
        Q = BreitWigner.shoot(transitionQ,0.010);

      break;

    case 2896:

        //KE = randen->BreitWigner(KE, 125);
        // Q = G4RandGauss::shoot(transitionQ,0.125);
        Q = BreitWigner.shoot(transitionQ,0.125);

      break;

    case 2186:

        //KE = randen->BreitWigner(KE, 200);
        // Q = G4RandGauss::shoot(transitionQ,0.2);
        Q = BreitWigner.shoot(transitionQ,0.2);

      break;

    default: 

        //KE = randen->BreitWigner(KE, 200);
        // Q = G4RandGauss::shoot(transitionQ,0.2);
        Q = BreitWigner.shoot(transitionQ,0.2);

     break;

  }
  //*/
  if (Q <= 0 || Q > 20.6) {

    Q = BreitWignerQ(transitionQ);
  }
  return Q;
}


void G4NeutronDecay::DumpNuclearInfo()
{
  G4cout << " G4NeutronDecay for parent nucleus " << GetParentName() << G4endl;
  G4cout << " decays to " << GetDaughterName(0) << " + " << GetDaughterName(1)
         << " with branching ratio " << GetBR() << "% and Q value "
         << transitionQ << G4endl;
}

