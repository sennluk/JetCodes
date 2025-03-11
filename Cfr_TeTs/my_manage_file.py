#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 10:37:20 2024

@author: lsenni
Libreria per la gestione dei file, cartelle, 
per salvataggio plots e dettagli delle eleaborazioni
"""
import os   # per la gestione delle cartelle di salvataggio
import sys # Per salvare gli output dei 'print'
from io import StringIO

# DTE3 shot list: 104991, 992, 994,995,999   
# 104520,521,522,523,524,526
##################################################
def create_folder(d):
    path = d['mypath']
    try:
        # Create a folder at the specified path
        os.makedirs(path, exist_ok = True)
        print(f"Folder created at '{path}'")
    except OSError as e:
        # If the folder already exists or there's any other error
        print(f"Failed to create folder at '{path}': {e}")
        os.makedirs(path+'_copy')

##################################################

def save_param(d):
    shot = d['shot']
    with open(d['mypath']+f'{shot}_A0_Param.txt', 'w') as f:
        print(d, file=f)
    
##################################################


def capture_print_output():
    # Save the current stdout
    original_stdout = sys.stdout
    # Create a StringIO object to capture the output
    captured_output = StringIO()
    # Redirect stdout to the StringIO object
    sys.stdout = captured_output
    # Perform the code execution here
    # This is where your code should be executed
    # Any print statements in your code will be captured
    # Example:
    print("This is a test.")
    print("Hello, world!")
    # Restore the original stdout
    sys.stdout = original_stdout
    # Get the captured output as a string
    captured_output_str = captured_output.getvalue()
    # Close the StringIO object
    captured_output.close()
    # Now you can do whatever you want with the captured output
    # For example, save it to a file
    with open("output.txt", "w") as f:
        f.write(captured_output_str)


    
##################################################


# Funzione con controlli di esistenza cartella etc....da migliorare
# def create_folder(d):
#     path = d['mypath']
#     try:
#         # Check if the folder already exists
#         if os.path.exists(path):
#             # Ask the user if they want to overwrite the folder
#             response = input(f"The folder '{path}' already exists. Do you want to save plots in the same folder (Possible overwriting) (y/n): ")
#             if response.lower() != 'y':   # response.lower(): This part of the expression takes the string stored in the variable
#                 d['savefigs'] = 0
#                 print('You are not saving the plots')
#             # !=: This is the "not equal to" operator in Python. It checks whether the expression on its left side is not equal to
#                 # try: 
#                 #     path = (path+'_2')
#                 #     os.makedirs(path)
#                 #     print(f"Folder created at '{path}'")
#                 #     d['mypath'] = path
#                 # except OSError as e:# If there's any other error
#                 #     print(f"Failed to create folder at '{path}': {e}")
#                 return
#             # If the user confirms, remove the existing folder
#             else:
#                 print("You chosed to save/overwrite plots...")
#                 # os.rmdir(path)
#                 return
                
#         # Create the folder at the specified path
#         os.makedirs(path)
#         print(f"Folder created at '{path}'")
#     except OSError as e:
#         # If there's any other error
#         print(f"Failed to create folder at '{path}': {e}")


###############################################################   
# Save all the outputs of the 'prinmt' command executed above
# def create_text_file(file_path, content=""):
#     filepath = d['mypath']
#     try:
#         with open(file_path, 'w') as file:
#             file.write(content)
#         print(f"File '{file_path}' created successfully.")
#     except Exception as e:
#         print(f"Error occurred while creating the file: {e}")

# Define a function to redirect stdout to a file
def save_print_output_to_file(d, vars):
    path = d['mypath']
    shot = d['shot']
    filepath = (f'{path}/{shot}_details')
    sys.stdout = open(filepath, 'w')

# Define a function to restore stdout back to the default
def restore_stdout():
    sys.stdout.close()
    sys.stdout = sys.__stdout__
    
    
 