//Class created by Afroditi Papadopoulou (apapadop@mit.edu)

// _________________________________________________________________________________________________________________________________________________________________________________________________

#ifndef BACKTRACKERTRUTHMATCH_CXX
#define BACKTRACKERTRUTHMATCH_CXX

#include "BackTrackerTruthMatch.h"
// _________________________________________________________________________________________________________________________________________________________________________________________________

void BackTrackerTruthMatch::MatchToMCParticle(const art::Handle<std::vector<recob::Hit> >& hit_handle, const art::Event& e, std::vector<art::Ptr<recob::Hit> >& trk_hits_ptrs)
{

	art::InputTag MCParticleModuleLabel { "largeant" };
	art::InputTag BackTrackerLabel { "gaushitTruthMatch" };
	art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hit_handle,e,BackTrackerLabel);

	std::unordered_map<int,double> trkide;
	std::vector<art::Ptr<simb::MCParticle>> particle_vec;
	std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

	double max_dQinTruthMatchedHits = -1., dQinAllHits = 0.;
	std::unordered_map<int,double> trackid_dQinTruthMatchedHits;

	// Loop only over the recob::Hits
	for (int i_h = 0; i_h < int(trk_hits_ptrs.size()); ++i_h) {

		float dQinHit = trk_hits_ptrs[i_h]->Integral(); // Charge deposition of a recob::Hit
		dQinAllHits += dQinHit; // Total charge deposition of all recob::Hit

		particle_vec.clear();  // vector of MCParticles
		match_vec.clear();
		particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);

		// To avoid matching the same particle more than once, we introduce a vector of matched TrackId-s
		// and require that the matched particle has not been matched already for this recob::Hit

		std::vector <int> ParticlesMatchedInThisHit;

		// Loop over MCParticles that match this hit and ask which one deposited the most energy
		for(int i_p = 0; i_p < int(particle_vec.size()); ++i_p) {
			art::Handle<std::vector<simb::MCParticle>> mcparticle_handle; 
			e.getByLabel(MCParticleModuleLabel,mcparticle_handle);
			//art::FindOneP<simb::MCTruth> MCParticleToMCTruth(mcparticle_handle,e,MCParticleModuleLabel);
			//const art::Ptr<simb::MCTruth> mctruth = MCParticleToMCTruth.at( particle_vec[i_p].key());
                        //std::cout<<"MC truth pdg: "<<mctruth->GetParticle()->PdgCode()<<std::endl;
			//if (mctruth->Origin() != 1) { continue; }

			float Edep_particle = match_vec[i_p]->energy; // Energy deposited by ionization by this track ID [MeV]

			trkide[particle_vec[i_p]->TrackId()] += Edep_particle; // Store energy [MeV] deposited by track id

			ftote += Edep_particle; // Calculate total energy deposited by all recob::Hits

			if (trkide[particle_vec[i_p]->TrackId()] > fmaxe) { // Keep track of maximum

				fmaxe = trkide[ particle_vec[i_p]->TrackId() ];
				fmaxp_me = particle_vec[i_p];
				fMCParticleID = particle_vec[i_p].key();

				if (!ParticleAlreadyMatchedInThisHit( ParticlesMatchedInThisHit , (int)particle_vec[i_p]->TrackId()) ) {

					ParticlesMatchedInThisHit.push_back( (int)particle_vec[i_p]->TrackId() );
					trackid_dQinTruthMatchedHits[ particle_vec[i_p]->TrackId() ] += dQinHit; // store the integral on the hit by the track id
					max_dQinTruthMatchedHits = trackid_dQinTruthMatchedHits[ particle_vec[i_p]->TrackId() ];

				}
		
			} // End of the if-statement

		} // End of the loop over particles per hit

	} // End of the loop over the hits


	fpurity = max_dQinTruthMatchedHits / dQinAllHits;

	fcompleteness = fmaxe / ftote;

	double kMin_dQ_inTruthMatchedHits = 0.;

	if ( fpurity < kMin_dQ_inTruthMatchedHits) { 

		art::Ptr< simb::MCParticle > fmaxp_me_empty;
		fmaxp_me = fmaxp_me_empty;	

	}
        
	return;
}
// _________________________________________________________________________________________________________________________________________________________________________________________________

art::Ptr< simb::MCParticle > BackTrackerTruthMatch::ReturnMCParticle() {

	return fmaxp_me;
}
// _________________________________________________________________________________________________________________________________________________________________________________________________

bool BackTrackerTruthMatch::ParticleAlreadyMatchedInThisHit(std::vector<int> AlreadyMatched_TrackIDs ,int cTrackID )
{
	// to avoid from matching the same particle more than once, we introduce a vector of matched TrackId-s for each hit 
	// and require that the matched particle has not been mathced already for this hit  

	for (auto trk_id:AlreadyMatched_TrackIDs) { 

		if (trk_id == cTrackID) return true; 
	}

	return false; 
}
// _________________________________________________________________________________________________________________________________________________________________________________________________

double BackTrackerTruthMatch::ReturnPurity() {

	return fpurity;
}
// _________________________________________________________________________________________________________________________________________________________________________________________________

double BackTrackerTruthMatch::ReturnCompleteness() {

	return fcompleteness;
}
// _________________________________________________________________________________________________________________________________________________________________________________________________

int BackTrackerTruthMatch::ReturnMCParticleID() {

	return fMCParticleID;
}

// _________________________________________________________________________________________________________________________________________________________________________________________________

double BackTrackerTruthMatch::ReturnTrueAssDepositedEnergy() {

	return fmaxe;
}

// _________________________________________________________________________________________________________________________________________________________________________________________________

double BackTrackerTruthMatch::ReturnTotalTrueDepositedEnergy() {

	return ftote;
}

// ____________________________________________________________________________________________________________________________________________________________________________________________________

#endif
