What do all of these fhiccls do and what do they mean?

config_single_proton_filter.fcl -> configures all the common things like space charge
config_single_proton_filter_wgt.fcl -> just like above but also configures the GENIE reweighting
my_particle_recovery.fcl -> No clue
my_wgt.fcl  -> No Clue

run_single_proton_*_filter.fcl-> configs specific to whatever data produt
run_single_proton_*_crt_filter.fcl -> configs specific to whatever data product, but we do run the crt

run_single_proton_overlay_filter_wgt.fcl -> any fhicl with wgt at the end (should only be for the overlay products) runs the genie reweighting
