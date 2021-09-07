#include "G4InFlightProcess.hh"
// #include "G4RadioactiveDecaymessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4ParticleChangeForRadDecay.hh"
#include "G4ITDecay.hh"
#include "G4BetaDecayType.hh"
#include "G4BetaMinusDecay.hh"
#include "G4BetaPlusDecay.hh"
#include "G4ECDecay.hh"
#include "G4AlphaDecay.hh"
#include "G4ProtonDecay.hh"
#include "G4NeutronDecay.hh"
#include "G4Neutron2Decay.hh"
#include "G4VDecayChannel.hh"
#include "G4NuclearDecay.hh"
#include "G4RadioactiveDecayMode.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4BetaDecayType.hh"
#include "Randomize.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NuclearLevelManager.hh"
#include "G4NuclearLevelStore.hh"
#include "G4ThreeVector.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Alpha.hh"
#include "G4Proton.hh"

#include "G4HadronicProcessType.hh"
#include "G4HadronicException.hh"
#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"

#include <vector>
#include <sstream>
#include <algorithm>
#include <fstream>

using namespace CLHEP;

G4InFlightProcess::G4InFlightProcess(const G4String& processName)
 : G4VRestContinuousDiscreteProcess(processName, fDecay)
{
	G4cout << "G4InFlightProcess Constructor: processName = " << processName << G4endl;
}

G4InFlightProcess::~G4InFlightProcess()
{};

G4double PostStepGetPhysicalInteractionLength(const G4Track &track, G4double previousStepSize, G4ForceCondition *condition)
{
	G4double interactionLength = G4VRestContinuousDiscreteProcess.PostStepGetPhysicalInteractionLength(track, previousStepSize, condition);
	G4cout << "Post Step Get Physical Interaction Length" << G4endl;
	return interactionLength;
}

G4double AlongStepGetPhysicalInteractionLength (const G4Track &track, G4double previousStepSize, G4double currentMinimumStep, G4double &currentSafety, G4GPILSelection *selection)
{
	G4double interactionLength = G4VRestContinuousDiscreteProcess.AlongStepGetPhysicalInteractionLength(track, previousStepSize, currentMinimumStep, currentSafety, selection);
	G4cout << "Along Step Get Physical Interaction Length" << G4endl;
	return interactionLength;
}

G4double AtRestGetPhysicalInteractionLength (const G4Track &, G4ForceCondition *condition)
{
	G4double interactionLength = G4VRestContinuousDiscreteProcess.AtRestGetPhysicalInteractionLength(track, condition);
	G4cout << "At Rest Get Physical Interaction Length" << G4endl;
	return interactionLength;
}

G4double GetMeanLifeTime(const G4Track &aTrack, G4ForceCondition *condition)
{
	G4cout << "GetMeanLifeTime entered." << g4endl;
	return 0.1*ns;
}

G4double GetContinuousStepLimit(const G4Track &aTrack, G4double previousStepSize, G4double currentMinimumStep, G4double &currentSafety)
{
	G4cout << "GetContinuousStepLimit entered." << G4endl;
	return 1.0;
}

G4double GetMeanFreePath(const G4Track &aTrack, G4double previousStepSize G4ForceCondition *condition)
{
	G4cout << "GetMeanFreePath entered." << G4endl;
	return 2.72*ns;
}