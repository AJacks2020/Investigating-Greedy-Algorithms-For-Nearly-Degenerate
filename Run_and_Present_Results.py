import subprocess
import matplotlib.pyplot as plt

def run_experiment(compile_command:str="ifort", if_Fortran_pre_compiled:bool= False):
    '''
        Function that runs the full experiments and plots results.
        The results of the Fortran code are not ingested directly as
        I want to keep a copy of the unchanged and full outputs of the
        Fortran so they are written to a file - in the Fortran program
        - and this function reads from there.

        INPUTS
        compile_command         : The string to use to call your compiler of choice
                                  to compile the Fortran. Deafults to using the intel
                                  OpenAPI compiler with the command "ifort". 
            
        if_Fortran_pre_compiled : if the Fortran code - "oldGreedyProgram.f90"
                                  - has already been compiled into an executable.
                                  Defaults to the Fortran NOT being already compiled.

        OUTPUTS
        NONE

        SIDE EFFECTS
        Displays a graph of ... and writes ??? to a txt file.
    '''

    # If the Fortran is not already compiled, compiles it
    if not if_Fortran_pre_compiled:
        subprocess.call([compile_command,"oldGreedyProgram.f90"])

    # Runs the Fortran code
    subprocess.call(["oldGreedyProgram.exe"])

    # Reads the stored results
    with open("GreedyOutputFile.txt") as f: Raw_Fortran_results = f.read()

    # Reformats the raw data into a list of floats
    formatted_Fortran_results = list(filter(None,
                                            Raw_Fortran_results
                                            .replace('\n',' ')# is this valid, ths splitting?
                                            .split(' ')
                                            ))
    usable_Fortran_results = [float(x) for x in formatted_Fortran_results]

    # Plots the values output by the Fortran code in a graph
    ratio_values = [0.05*y for y in range(1, 20)]
    plt.plot(ratio_values, usable_Fortran_results)

    # Adds axis labels to the graph
    plt.xlabel("Ratio of optimal to nearly optimal tour length") # this may not be correct
    plt.ylabel("greedy algorithm success probability")

    plt.show()


if __name__ == "__main__":
    run_experiment()
