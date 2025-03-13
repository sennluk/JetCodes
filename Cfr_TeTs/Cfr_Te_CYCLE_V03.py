"""
Created 29/01/2025 - lsenni - 
The 'CYCLE version is equal to the MAIN version except that introduce the cycle over 
the shots of interests, creating a Database - (Python dictionary) 
with 'd', and 'vers' correspondig to every shot
Save the DB in .json
V03: inizio la pulizia e il compattamento
"""
import matplotlib.pyplot as plt
import Cfr_Te_MAIN_V03 as main
import Cfr_Te_LIB_V03 as mye
import Cfr_Te_PLOTS_V03 as graph
import time
import pickle
# %reset

plt.close('all')  

# Selection of the shots for the Database
shots =  [104520, 104521, 104522, 104523, 104524, 104525, 104526, 
            104547, 104548 ,104549, 104550, 104551, 104553, 104554,104555,104558,104559,104560,
            104574, 104575, 104990, 104991,104994] 
# 104995: PSI not present
# 
DB = {} # Empty general database 

saveDB = 0      # 0: don't save, 1: save the DB
save_figures = 0  # 0: don't save, 1: save the figures
timestr = time.strftime("%Y%m%d-%H%M%S") # Date format for the folder name generation

fig,ax=plt.subplots(1, num = 'ECE vs HRTS')
for shot in shots:
    # plt.close('all')  
    JPN = shot
    d,vars = main.pippo(JPN,save_figures)
    
    DB[shot] = {'par':d, 'var':vars} 

## Save the Dictionary in .json. If exist do not overwrite
n = len(shots)

if saveDB == 1:
    try:
        with open(f"DB_{n}_Shots_V02.pkl", "wb") as f:
            pickle.dump(DB, f)
    except FileExistsError:
        print(" 'DB.pkl' already exist!")

# if saveDB == 1:    
#     try:
#         with open(f"DB_{n}_Shots.json", "w") as f:    # 'x': non sovrascrive, 'w' sovrascrive
#             json.dump(DB, f, indent=4)     # Scrivi nel file
#     except FileExistsError:
#         print(" 'DB.json' already exist!")



## Per riaprire il dizonario
# with open("DB.json", "r") as f:
#     DB = json.load(f)

# import h5py

# # Salvataggio in formato HDF5
# with h5py.File('DB.h5', 'w') as f:
#     for key, value in DB.items():
#         f.create_dataset(key, data=value)

# # # Caricamento del file HDF5
# # with h5py.File('DB.h5', 'r') as f:
# #     loaded_db = {key: f[key][()] for key in f.keys()}


