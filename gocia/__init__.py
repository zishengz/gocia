# Initialize a flag to track if the message has been printed
_initialized = False

# Define a function to print the message once
def initialize_module():
    global _initialized
    if not _initialized:
        # logo art adapted from https://patorjk.com/software/taag/#p=testall&f=Graffiti&t=GOCIA
        print("""
   ░▒▓██████▓▒░ ░▒▓██████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░░▒▓██████▓▒░  
  ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
  ░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
  ░▒▓█▓▒▒▓███▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░▒▓████████▓▒░ 
  ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░  
   ░▒▓██████▓▒░ ░▒▓██████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░  """)
        print('   Global Optimizer for Clusters, Interfaces, and Adsorbates')
        print('                Copyright © 2024 Zisheng Zhang.')
        print('     Cite: PCCP, 2025, 27, 696-706. DOI:10.1039/D4CP03801K')

        _initialized = True

# Call the function during the first import
initialize_module()

# The rest of your module initialization can go here if needed









