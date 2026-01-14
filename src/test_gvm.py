from neovm import tissue, database
import output_to_input as oti
import networkx as nx

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--nsim', type=int, required=True)
parser.add_argument('--tmax', type=float, required=True)
parser.add_argument('--nmax', type=int, required=True)
parser.add_argument('--sigma', type=float, required=True)
parser.add_argument('--drate', type=float, required=True)
parser.add_argument('--tens', type=float, required=True)
parser.add_argument('--lays', type=int, required=True)
parser.add_argument('--g0b', type=float, required=True)
args = parser.parse_args()

def parameters(t):
    t.Nsim = args.nsim # number of simulations
    t.Tmax = args.tmax # max termination time
    t.Nmax = args.nmax # max number of cells
    t.sig = args.sigma # edge fluctuation amplitude 'sigma'
    t.division_rate = args.drate # division rate 'k_div'
    t.g0_layer = args.tens # tension between live and necrotic layers 'Gamma_LNI'
    t.live_layers = args.lays # number of live layers 'lambda'
    t.g0_boundary = args.g0b # surface tension 'Gamma'
    t.h = 0.001 # time-integration step
    t.dg0 = 1 # baseline line tension
    t.kM = 1 # myosin turnover rate
    t.kV = 100 # volume compressibility modulus
    t.outFreq = 10.0 # output frequency

# Create the graph database
G = nx.MultiDiGraph()
# Initialize the tissue
t = tissue('/path/to/test/input.vt3d')
# Initialize the database
db = database(G, t)
# Import the model parameters
parameters(t)
# Define the output directory
outFolder = '/path/to/test_output'
# Simulate the tissue
success, last_file, new_file, last_fileCount, last_time, last_timeCount, nET, nTE, nD = t.simulate(G, db, outFolder)

# If anything goes wrong in the simulation: try again from the last stored step that was still correct
while not success:
    oti.out_to_in(outFolder, last_file, new_file)
    G = nx.MultiDiGraph()
    t=tissue('{}/{}'.format(outFolder, new_file))
    db=database(G, t)
    parameters(t)
    t.fileCount = last_fileCount
    t.Time = last_time
    t.timeCount = last_timeCount
    t.ET_count = nET
    t.TE_count = nTE
    t.D_count = nD
    success, last_file, new_file, last_fileCount, last_time, last_timeCount, nET, nTE, nD = t.simulate(G, db, outFolder)
