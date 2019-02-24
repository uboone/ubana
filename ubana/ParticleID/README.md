# Particle ID module

This package is intended to be used for analyses on MicroBooNE. It's currently under active development. The intention of this framework is to centralise the PID efforts, and make it easy to create and implement additional algorithms in an experiment agnostic way. In addition to this, we have implemented a number of new particle ID algorithms in addition to the older, more tested algorithms.

## Algorithm details

__Bragg_Likelihood Estimator__
This algorithm uses B. Baller's theory, along with reconstructed hits in the Bragg peak of each track to estimate the likelihood of each data point belonging to each particle species, and then combines those values to get a single likelihood value for the track.

This can be thought of as an extension to the chi-square method. The underlying dE/dx distribution is assumed to be a Landau-Gaussian with widths measured in simulation and data rather than the Gaussian assumption of the chi-square method. There are a number of other additions which are expected to help give better separation.

__PIDA__
This makes use of B. Baller's home-brewed PIDA method which was used to great success on ArgoNEUT. There are thee implementations available here:
1. Mean (B. Baller)
2. Median (T. Yang + V. Meddage)
3. Kernel density estimator (A. Lister)

__Truncated Mean dE/dx Versus Track Length__
Makes use of D. Caratelli's "Truncated Mean" algorithm to plot the mean dE/dx of a track versus its length.

__Deposited Energy vs Energy By Range__
Calculates deposited energy, and compares against the energy by range under different assumptions for separation.

## Package Information and Dependencies

This package has been validated against uboonecode v06_26_01_13, but it should be "simple" to take this to later production versions. It should also be reasonably east to port this to develop releases, with one caveat: New data files as of ubooneocde v06_78_00 are using the 2d-deconvolution, and the plan is to move to using this for simulated data in the near future. This would require significantly more work to re-calibrate the data products and re-measure the Landau and Gaussian widths. **If you're using develop proceed with caution**. 

For v06_26_01_13, this package has two dependencies: 
- lardataobj feature branch `feature/feature/kduffy_pidrefactor_v1_11_00_04` 
  - This contains an extension of the structure of the anab::ParticleID class, as shown below.
- larana feature branch `origin/feature/alister1_TruncatedMeanPort` 
  - This contains a port of D. Caratelli's TruncatedMean class.

It is likely that if using a newer version of the production code will mean you do not need the larana feature branch.

### Structure of Updated anab::ParticleID Class

```cpp
struct sParticleIDAlgScores {

  std::string fAlgName;
  kVariableType fVariableType;
  int fAssumedPdg;
  float fValue;
  geo::PlaneID fPlaneID;

}

```

which holds the output of any generic PID algorithm. Here, `kVariableType` is an enum which can take the following values:

```cpp
enum kVariableType{
  kGOF,
  kLikelihood,
  kLikelihood_fwd,
  kLikelihood_bwd,
  kLogL,
  kLogL_fwd,
  kLogL_bwd,
  kScore,
  kPIDA,
  kdEdxtruncmean,
  kdQdxtruncmean,
  kTrackLength,
  kEdeposited,
  kEbyRange,
  kNotSet
}
```

Below is a minimal example of how to access the result of a given PID algorithm (in this case, the result of Bragg_Likelihood_Estimator under a muon assumption):

```cpp
art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, e, fPIDLabel);
if (!trackPIDAssn.isValid()){
  std::cout << "[ParticleIDValidation] trackPIDAssn.isValid() == false. Skipping track." << std::endl;
  continue;
}
    
std::vector<art::Ptr<anab::ParticleID>> trackPID = trackPIDAssn.at(track->ID());
if (trackPID.size() == 0){
  std::cout << "[ParticleIDValidation] No track-PID association found for trackID " << track->ID() << ". Skipping track." << std::endl;
  continue;
}

std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
    
double Bragg_fwd_mu = -999;

// Loop through AlgScoresVec and find the variables we want
for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){

  anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
  int planeid = AlgScore.fPlaneID.Plane;

  if (planeid < 0 || planeid > 2){
    std::cout << "[ParticleIDValidation] No information for planeid " << planeid << std::endl;
    continue;
  }
  
  if (AlgScore.fAlgName == "BraggPeakLLH"){
    if (anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood_fwd){
       if (AlgScore.fAssumedPdg == 13)   Bragg_fwd_mu = AlgScore.fValue;
    }
  }
}
```
