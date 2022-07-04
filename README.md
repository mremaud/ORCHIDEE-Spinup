# ORCHIDEE-Spinup
The repository contain Python scripts aiming to initialize the forested variables (e.g. Diameter, height) in the ORCHIDEE Land Surface Model.

More specifically, the scripts nudge the ORCHIDEE initial state of the forest toward observations from Pucher Inventory (https://www.mdpi.com/2072-4292/14/2/395). All scripts require at least 100 years of an ORCHIDEE simulation with a balanced carbon content (NBP close to 0). These years are obtained after a forest clearcut.

Nudged-spinup.py: Nudge the simulated diameter toward the diameter from the Pucher Inventory. It creates new restart files. 

Nudged-spinup_4ac.py: The script is similar to Nudge-spinup.py except it was adapted to an ORCHIDEE version that simulates 4 diameter classes over forested areas. Nudge the fraction of diameter classes toward the one from the Pucher Inventory. The fraction of vegetation is taken from the ORCHIDEE simulation.
