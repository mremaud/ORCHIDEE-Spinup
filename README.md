# Forest-Initialization
The repository contain Python scripts aiming to initialize the forested variables (e.g. Diameter, height) in the ORCHIDEE Land Surface Model.

More specifically, the scripts nudge the ORCHIDEE initial state of the forest toward observations from Pucher Inventory (https://www.mdpi.com/2072-4292/14/2/395). 

PYTHON SCRIPTS:

(1) Read_inventory.py: Read the observed variables (Diameter, fraction of age classes, fraction of vegetation) located in the repertory ForestInventory. This repertory contains data of the Pucher forest inventory and can be downloaded through a link written in the publication https://www.mdpi.com/2072-4292/14/2/395. The script generates the file DIAMETER.nc that will be used as an input to the script Nudged-spinup_*.py 

(2) This script nudges the simulated diameter toward the observations from the Pucher Inventory. It creates new restart files. n refers to the number of age classes in the ORCHIDEE version that simulates n diameter classes over forested areas. The number of age classes is prescribed using a def file in the parameter lists.

(3) Map-spinup.py: Draw maps that compare the nudged state of the forest against the randomized state of the forest at the end of the spinup simulation.

Requirement: All scripts require at least 100 years of an ORCHIDEE simulation with a balanced carbon content (NBP close to 0). These years are obtained after a forest clearcut.

(4) A documentation
