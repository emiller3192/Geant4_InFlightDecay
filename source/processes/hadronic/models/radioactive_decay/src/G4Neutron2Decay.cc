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

#include "G4Neutron2Decay.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <iomanip>

#include <fstream>
#include <cmath>


G4Neutron2Decay::G4Neutron2Decay(const G4ParticleDefinition* theParentNucleus,
                            const G4double& branch, const G4double& Qvalue,
                            const G4double& excitationE)
 : G4NuclearDecay("neutron2 decay", Neutron2, excitationE), transitionQ(Qvalue)
{
  SetParent(theParentNucleus);  // Store name of parent nucleus, delete G4MT_parent
  SetBR(branch);

  SetNumberOfDaughters(3);
  G4IonTable* theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  G4int daughterZ = theParentNucleus->GetAtomicNumber();
  G4int daughterA = theParentNucleus->GetAtomicMass() - 2;
  SetDaughter(0, theIonTable->GetIon(daughterZ, daughterA, excitationE) );
  SetDaughter(1, "neutron");
  SetDaughter(2, "neutron");
}


G4Neutron2Decay::~G4Neutron2Decay()
{}


G4DecayProducts* G4Neutron2Decay::DecayIt(G4double)
{
  // Find probability of correlation and randomly assign boolean with proper
  // weight. 
  bool Correlated = true;
  G4int mode = 0;

  std::fstream angles;
  angles.open("/home/eric/Desktop/angles.txt", std::ios::out | std::ios::app);

  // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)  
  CheckAndFillParent();

  // Fill G4MT_daughters with neutron and residual nucleus (stored by SetDaughter)  
  CheckAndFillDaughters();

  G4double neutronMass = G4MT_daughters[1]->GetPDGMass();
  // Excitation energy included in PDG mass
  G4double nucleusMass = G4MT_daughters[0]->GetPDGMass();

  G4double parentMass = G4MT_parent->GetPDGMass();

  // Q value was calculated from atomic masses.
  // Use it to get correct neutron energy.

  // Redo the cmMomentum appropriately

  //G4double switchQ = 1000. * transitionQ;
  G4double Q = 0;

  G4double corr_prob = 0.5;  // used to tune the 2n decay mode ratio

  if (G4UniformRand() >= corr_prob) Correlated = false;

  //G4double cmMomentum = std::sqrt(Q*(Q + 2.*neutronMass)*
  //                              (Q + 2.*nucleusMass)*
  //                              (Q + 2.*neutronMass + 2.*nucleusMass) )/
  //                              (Q + neutronMass + nucleusMass)/2.; 

  // Set up final state
  // parentParticle is set at rest here because boost with correct momentum 
  // is done later
  G4DynamicParticle parentParticle(G4MT_parent, G4ThreeVector(0,0,0), 0.0);
  G4DecayProducts* products = new G4DecayProducts(parentParticle);

  if (!Correlated) {

      // Electron, neutrino and daughter nucleus energies
    Q = G4RandGauss::shoot(transitionQ,0.2);
    G4double n1KE = Q*G4UniformRand();
    G4double n1Momentum = std::sqrt(n1KE*(n1KE + 2.*neutronMass) );
    G4double n2KE = Q - n1KE;

    G4double cosThetaENu = 2.*G4UniformRand() - 1.;
    G4double n1TE = neutronMass + n1KE;
    G4double n2Energy = ((Q - n1KE)*(parentMass + nucleusMass - n1TE)
            - n1Momentum*n1Momentum)/(parentMass - n1TE + n1Momentum*cosThetaENu)/2.;

    // Neutron 1 4-vector, isotropic angular distribution
    G4double cosTheta = 2.*G4UniformRand() - 1.0;
    G4double sinTheta = std::sqrt(1.0 - cosTheta*cosTheta);

    G4double phi = twopi*G4UniformRand()*rad;
    G4double sinPhi = std::sin(phi);
    G4double cosPhi = std::cos(phi);

    G4ParticleMomentum n1Direction(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta);

    if (n1Direction.z() < 0.) n1Momentum = 0.;

    G4DynamicParticle* dynamicNeutron1
      = new G4DynamicParticle(G4MT_daughters[1], n1Direction*n1Momentum);
    products->PushProducts(dynamicNeutron1);

    // Neutron 2 4-vector
    G4double sinThetaENu = std::sqrt(1.0 - cosThetaENu*cosThetaENu);
    phi = twopi*G4UniformRand()*rad;
    G4double sinPhiNu = std::sin(phi);
    G4double cosPhiNu = std::cos(phi);

    G4ParticleMomentum n2Direction;
    n2Direction.setX(sinThetaENu*cosPhiNu*cosTheta*cosPhi -
                     sinThetaENu*sinPhiNu*sinPhi + cosThetaENu*sinTheta*cosPhi);
    n2Direction.setY(sinThetaENu*cosPhiNu*cosTheta*sinPhi +
                     sinThetaENu*sinPhiNu*cosPhi + cosThetaENu*sinTheta*sinPhi);
    n2Direction.setZ(-sinThetaENu*cosPhiNu*sinTheta + cosThetaENu*cosTheta);

    if (n2Direction.z() < 0.) n2Energy = 0.;


    G4DynamicParticle* dynamicNeutron2
      = new G4DynamicParticle(G4MT_daughters[2], n2Direction*n2Energy);
    products->PushProducts(dynamicNeutron2);

    // Daughter nucleus 4-vector
    // p_D = - p_e - p_nu
    G4DynamicParticle* dynamicDaughter =
      new G4DynamicParticle(G4MT_daughters[0],
                            -n1Direction*n1Momentum - n2Direction*n2Energy);
    products->PushProducts(dynamicDaughter);

    mode = 1;

    angles << mode << "    " << n1Direction.x() << "  " << n1Direction.y() << "  " << n1Direction.z() << "  "
  << n2Direction.x() << "  " << n2Direction.y() << "  " << n2Direction.z() << "  " 
  << acos(n1Direction.x()*n2Direction.x() + n1Direction.y()*n2Direction.y() + n1Direction.z()*n2Direction.z()) << "  "
  << n1KE << "  " << n2KE << "  " << std::endl;
      // UNCORRELATED:  Direction of second neutron is random but =/= to the 
      //                direction of the first neutron.  The nucleus then needs
      //                to be in the opposite direction in the COM frame so that
      //                the COM stays at 0.

  // This method should be fine now: 8 July 2019

  }
  else {

    // Put correlated mode here

    Q = G4RandGauss::shoot(transitionQ,0.2);

    G4double cmMomentum = std::sqrt(Q*(Q + 4.*neutronMass)*
                                (Q + 2.*nucleusMass)*
                                (Q + 4.*neutronMass + 2.*nucleusMass) )/
                                (Q + 2.*neutronMass + nucleusMass)/2.;

    G4double costheta = 2.*G4UniformRand()-1.0;
    G4double sintheta = std::sqrt(1.0 - costheta*costheta);
    G4double phi  = twopi*G4UniformRand()*rad;
    G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),
                          costheta);

    G4double KE = std::sqrt(cmMomentum*cmMomentum + nucleusMass*nucleusMass) - nucleusMass;

    G4DynamicParticle* daughterparticle =
    new G4DynamicParticle(G4MT_daughters[0], -1.0*direction, KE, nucleusMass);
    products->PushProducts(daughterparticle);

    Q = std::sqrt(cmMomentum*cmMomentum + 2.*neutronMass*2.*neutronMass)
              - 2.*neutronMass;

    G4double Q_2n = 0.100; // 1, 10, 100 keV 

    G4double cmMomentum2 = std::sqrt(Q_2n*(Q_2n + 2.*neutronMass)*
                          (Q_2n + 2.*neutronMass)*
                          (Q_2n + 2.*neutronMass + 2.*neutronMass) )/
                          (Q_2n + neutronMass + neutronMass)/2.;


    G4double KE2 = std::sqrt(cmMomentum2*cmMomentum2 + neutronMass*neutronMass) - neutronMass;

    // Use Q value of 2n pair decay (1keV, 10keV, 100keV) rather than what's used
    // to conserve Energy in the lab frame. 

    G4double costheta2 = 2.*G4UniformRand()-1.0;
    G4double sintheta2 = std::sqrt(1.0 - costheta*costheta);
    G4double phi2  = twopi*G4UniformRand()*rad;
    G4ThreeVector direction2(sintheta2*std::cos(phi2),sintheta2*std::sin(phi2),
                          costheta2);

    G4double norm = std::sqrt(std::pow(cmMomentum, 2) + std::pow(cmMomentum2, 2));

    G4double K = 0.5 * (cmMomentum/cmMomentum2);

    G4double Tn1 = KE2 * (1. + K*K + 2.*K*costheta2);  // transform kinectic energy to lab frame
    G4double Tn2 = KE2 * (1. + K*K - 2.*K*costheta2);

    G4ThreeVector final_direction1 = (cmMomentum/norm) * direction + (cmMomentum2/norm) * direction2;
    G4ThreeVector final_direction2 = (cmMomentum/norm) * direction - (cmMomentum2/norm) * direction2;

    if (final_direction1.z() < 0.) Tn1 = 0.;
    if (final_direction2.z() < 0.) Tn2 = 0.;

    G4DynamicParticle* dynamicNeutron1
      = new G4DynamicParticle(G4MT_daughters[1], final_direction1, Tn1, neutronMass);
    products->PushProducts(dynamicNeutron1);

    // if final_directionx.z() < 0, don't PushProducts

    G4DynamicParticle* dynamicNeutron2
      = new G4DynamicParticle(G4MT_daughters[2], final_direction2, Tn2, neutronMass);
    products->PushProducts(dynamicNeutron2);

    mode = 2;

    angles << mode << "    " << final_direction1.x() << "  " << final_direction1.y() << "  " << final_direction1.z() << "  "
     << final_direction2.x() << "  " << final_direction2.y() << "  " << final_direction2.z() << "  " 
     << acos(final_direction1.x()*final_direction2.x() + final_direction1.y()*final_direction2.y() + final_direction1.z()*final_direction2.z()) << "  "
     << Tn1 << "  " << Tn2 << "  " << KE << "  " << std::endl;

    // CORRELATED:  Two 2-body decay.  First the 2n pair (d) and daughter, then d -> n + n
    // Most the Q value is probably used up in the first decay, and the lab angles for the neutron
    // pair after splitting is small compared to initial direction.

  }
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


  //angles << n1Direction.X() << "  " << n1Direction.Y() << "  " << n1Direction.Z() << "  "
  //<< n2Direction.X() << "  " << n2Direction.Y() << "  " << n2Direction.Z() << "  " << acos(n1Direction.X()*n2Direction.X() + n1Direction.Y()*n2Direction.Y() + n1Direction.Z()*n2Direction.Z());
  angles.close();

  return products;
}


void G4Neutron2Decay::DumpNuclearInfo()
{
  G4cout << " G4Neutron2Decay for parent nucleus " << GetParentName() << G4endl;
  G4cout << " decays to " << GetDaughterName(0) << " + " << GetDaughterName(1) << " + " 
         << GetDaughterName(2) 
         << " with branching ratio " << GetBR() << "% and Q value "
         << transitionQ << G4endl;
}

