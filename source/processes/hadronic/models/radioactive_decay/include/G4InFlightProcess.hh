#ifndef G4InFlightProcess_h
#define G4InFlightProcess_h 1

// creating a dummy process to check how to implement a new process

#include <vector>
#include <map>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestContinuousDiscreteProcess.hh"
#include "G4ParticleChangeForRadDecay.hh"
// #include "G4RadioactiveDecaymessenger.hh"  

#include "G4NucleusLimits.hh"
#include "G4RadioactiveDecayRate.hh"
#include "G4RadioactiveDecayRateVector.hh"
#include "G4RIsotopeTable.hh"
#include "G4RadioactivityTable.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include "G4SystemOfUnits.hh"

class G4InFlightProcess : public G4VRestContinuousDiscreteProcess
{
	// public member functions
	public:
		G4InFlightProcess(const G4String& processName="InFlightProcess");
		~G4InFlightProcess();

		G4double PostStepGetPhysicalInteractionLength(const G4Track &track, G4double previousStepSize, G4ForceCondition *condition);
		// G4VParticleChange *PostStepDoIt(const G4Track &track, const G4Step &step);

		G4double AlongStepGetPhysicalInteractionLength (const G4Track &track, G4double previousStepSize, G4double currentMinimumStep, G4double &currentSafety, G4GPILSelection *selection);
		// G4VParticleChange *AlongStepDoIt(const G4Track &track, const G4Step &step);

		G4double AtRestGetPhysicalInteractionLength (const G4Track &, G4ForceCondition *condition);
		// G4VParticleChange *AtRestDoIt(const G4Track &track, const G4Step &step);

	protected:
		G4double GetMeanLifeTime(const G4Track &aTrack, G4ForceCondition *condition);
		G4double GetContinuousStepLimit(const G4Track &aTrack, G4double previousStepSize, G4double currentMinimumStep, G4double &currentSafety);
		// void SetGPILSelection(G4GPILSelection selection);
		// G4GPILSelection GetGPILSelection() const;
		G4double GetMeanFreePath(const G4Track &aTrack, G4double previousStepSize G4ForceCondition *condition);
}