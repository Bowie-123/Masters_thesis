import openpyxl
from openpyxl import Workbook
from Final_main_VSC_parallel import main, scarcityfunc
import sys
from mpi4py import MPI
import os
import numpy as np

BASE_PATH = "/scratch/leuven/347/vsc34764/python"
PARAMETERS_FILE = BASE_PATH + "/parameters.xlsx"
OUTPUT_DIRECTORY_FORMAT = BASE_PATH + "/output/Set{}___V{}__a{}_Rank{}"
COMPLETE_FILE_FORMAT = OUTPUT_DIRECTORY_FORMAT + "/complete{}.txt"
AVERAGE_RESULTS_FILE_FORMAT = OUTPUT_DIRECTORY_FORMAT + "/AVG.txt"

def read_parameters_from_excel(file_path):
    workbook = openpyxl.load_workbook(file_path)
    sheet = workbook.active

    runs_parameters = []
    for row in sheet.iter_rows(min_row=2, values_only=True):  # skip the header row
        if all(value is None for value in row): # sometimes None error even if field empty
            continue
        run_params = {
            "Vspeed": float(row[0]) if row[0] is not None else None,
            "alpha": float(row[1]) if row[1] is not None else None,
            "number_gens": int(row[2]) if row[2] is not None else None,
            "number_runs": int(row[3]) if row[3] is not None else None
        }
        # print(run_params)
        runs_parameters.append(run_params)

    return runs_parameters

def save_complete_run_to_file(run_data, file_path):
    with open(file_path, 'w') as f:
        f.write("Generation,avg_R,avg_virulence,avg_A,avg_B,scarcity,avg_daph_fitness,avg_parasite_fitness\n")
        for gen_idx, generation_result in enumerate(run_data):
            line = f"{gen_idx},"
            line += ",".join(str(generation_result[header]) for header in generation_result.keys())
            f.write(line + "\n")

def save_avg_results_to_file(avg_results, file_path):
    with open(file_path, 'w') as f:
        f.write("Generation,avg_R,avg_virulence,avg_A,avg_B,scarcity,avg_daph_fitness,avg_parasite_fitness\n")
        for gen_idx, generation_result in enumerate(avg_results):
            line = f"{gen_idx},"
            line += ",".join(str(generation_result[header]) for header in generation_result.keys())
            f.write(line + "\n")


def run_simulation(params, i, rank):
    generation_results, first_10_runs = main(params)

    output_directory = OUTPUT_DIRECTORY_FORMAT.format(
        i + 1, params["Vspeed"], params["alpha"], rank
    )
    os.makedirs(output_directory, exist_ok=True)

    # Save complete run
    complete_file = COMPLETE_FILE_FORMAT.format(
        i + 1, params["Vspeed"], params["alpha"], rank
    )
    save_complete_run_to_file(first_10_runs[i], complete_file)
    print(f"Saved complete run data to {complete_file}")

    # Save average results for every rank
    avg_results_file = AVERAGE_RESULTS_FILE_FORMAT.format(
        i + 1, params["Vspeed"], params["alpha"], rank
    )
    save_avg_results_to_file(generation_results, avg_results_file)
    print(f"Saved average results to {avg_results_file}")

if __name__ == "__main__":   
    comm = MPI.COMM_WORLD
    size = comm.Get_size() # Number of available nodes
    rank = comm.Get_rank()

    parameters_list = read_parameters_from_excel(PARAMETERS_FILE)
    
    # Calculate the number of simulations per node
    num_sims_per_node = len(parameters_list) // size
    start = rank * num_sims_per_node
    end = start + num_sims_per_node if rank != size - 1 else len(parameters_list)

    for idx in range(start, end):  # Assign simulations to each node
        run_params = parameters_list[idx]

        for i in range(run_params["number_runs"]):  # Run the simulation as many times as specified in number_runs
            params = {
                "num_daphnia": 1000,
                "num_parasites": 2000,
                "gamma": 0.5,
                "k": 0.5,
                "h": 0.2,
                "c": 0.75,
                "number_gens": run_params["number_gens"],
                "alpha": run_params["alpha"],
                "Vspeed": run_params["Vspeed"],
                "number_runs": run_params["number_runs"],
                "scarcity_function": scarcityfunc
            }

            try:
                run_simulation(params, i, rank)
            except Exception as e:
                print(f"Error occurred during simulation: {e}")